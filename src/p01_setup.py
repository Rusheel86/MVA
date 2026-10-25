# %% [markdown]
# ## 1 · Environment, helpers and platform probe
#
# All genotype-scale files live under an ephemeral scratch directory (`WORK`), never under `/kaggle/working`
# (which Kaggle persists with the notebook). Optional stages are wrapped so that a failure is logged and the run
# continues; the core path (VCF → ranked CSV) fails loudly.

# %%
import collections, gzip, io, itertools, json, math, os, platform, re, shutil, subprocess, sys, textwrap, threading
import time, traceback, zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

os.environ.setdefault("TORCH_DEVICE_BACKEND_AUTOLOAD", "0")    # TPU images: keep torch from auto-loading torch_xla
import numpy as np
import pandas as pd

T0 = time.time()
ON_KAGGLE = os.path.isdir("/kaggle/working")
OUT = Path("/kaggle/working") if ON_KAGGLE else Path("./mva_out").resolve()
WORK = Path("/kaggle/tmp/mva") if ON_KAGGLE else Path("./mva_work").resolve()
for _d in (OUT, OUT / "figures", WORK, WORK / "data", WORK / "ref", WORK / "tmp", WORK / "logs", WORK / "bin", WORK / "hf"):
    _d.mkdir(parents=True, exist_ok=True)

os.environ.update({
    "TMPDIR": str(WORK / "tmp"), "HF_HOME": str(WORK / "hf"), "MAMBA_ROOT_PREFIX": str(WORK / "mamba"),
    "HF_HUB_DISABLE_PROGRESS_BARS": "1", "HF_XET_FIXED_DOWNLOAD_CONCURRENCY": "4", "HF_HUB_DOWNLOAD_TIMEOUT": "60",
    "TF_FORCE_GPU_ALLOW_GROWTH": "true",
    "TF_CPP_MIN_LOG_LEVEL": "2", "TOKENIZERS_PARALLELISM": "false", "PYTHONUNBUFFERED": "1",
})
os.chdir(WORK)          # tools that drop temp files into the cwd must not leak into the saved notebook output
NCPU = os.cpu_count() or 4

RUN = collections.OrderedDict()     # stage -> dict(status, seconds, note)
S = {}                              # shared state between stages (results that later stages/report consume)


_LOG_SEQ = itertools.count()


def log(msg):
    print(f"[{time.time() - T0:7.0f}s] {msg}", flush=True)


def sh(script, name="cmd", check=True, capture=False, env=None, timeout=None):
    """Run a bash script with errexit+pipefail. stdout/stderr go to a scratch log (not to the saved notebook
    output). Returns captured stdout when capture=True. Raises RuntimeError with the log tail on failure."""
    logf = WORK / "logs" / f"{next(_LOG_SEQ):03d}_{re.sub(r'[^A-Za-z0-9_.-]+', '_', name)[:50]}.log"
    full_env = dict(os.environ)
    full_env.update(env or {})
    with open(logf, "w") as lf:
        p = subprocess.run(["bash", "-euo", "pipefail", "-c", script], stdout=subprocess.PIPE if capture else lf,
                           stderr=lf, text=True, env=full_env, timeout=timeout)
    if p.returncode != 0 and check:
        tail = "".join(open(logf, errors="replace").readlines()[-25:])
        raise RuntimeError(f"[{name}] exit {p.returncode}\n--- script ---\n{script[:600]}\n--- log tail ---\n{tail[-2500:]}")
    return p.stdout if capture else p.returncode


class stage:
    """`with stage('name'):` times a block and records the outcome. Optional stages swallow exceptions
    (logged with a short traceback) so that one failing evidence module never costs the core deliverables."""

    def __init__(self, name, required=False):
        self.name, self.required = name, required

    def __enter__(self):
        self.t = time.time()
        log(f"▶ {self.name}")
        return self

    def __exit__(self, et, ev, tb):
        dt = time.time() - self.t
        if et is None:
            RUN[self.name] = dict(status="ok", seconds=round(dt, 1), note="")
            log(f"✔ {self.name} ({dt:.0f}s)")
            return False
        note = f"{et.__name__}: {str(ev)[:300]}"
        RUN[self.name] = dict(status="FAILED" if self.required else "skipped", seconds=round(dt, 1), note=note)
        log(f"✘ {self.name} ({dt:.0f}s) — {note}")
        if self.required:
            return False
        traceback.print_exception(et, ev, tb, limit=3)
        return True


def fetch(url, dest, min_bytes=1000, tries=5):
    """Resumable HTTPS download with retries; skips files that are already complete."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size >= min_bytes and not Path(str(dest) + ".part").exists():
        return dest
    part = Path(str(dest) + ".part")
    sh(f"curl -fsSL --retry {tries} --retry-delay 5 --retry-all-errors --speed-limit 500000 --speed-time 180 -C - "
       f"-o '{part}' '{url}'", name=f"fetch_{dest.name}")
    if part.stat().st_size < min_bytes:
        raise RuntimeError(f"download too small: {url} ({part.stat().st_size} bytes)")
    part.rename(dest)
    return dest


def fetch_parallel(url, dest, n=8, min_bytes=1000, budget_s=None):
    """Large-file download over n parallel byte ranges (each resumable); falls back to `fetch` without range support.
    `budget_s` caps the wall time so a throttled mirror fails fast and the caller's fallback source is used instead."""
    dest = Path(dest)
    t_end = time.time() + budget_s if budget_s else None
    if dest.exists() and dest.stat().st_size >= min_bytes:
        return dest
    hdr = sh(f"curl -sIL --retry 3 '{url}'", name=f"head_{dest.name}", capture=True, check=False) or ""
    sizes = re.findall(r"(?im)^content-length:\s*(\d+)", hdr)
    size = int(sizes[-1]) if sizes else 0
    probe = sh(f"curl -sL -r 0-0 -o /dev/null -w '%{{http_code}}' '{url}'", name=f"probe_{dest.name}", capture=True, check=False)
    if size < 200_000_000 or (probe or "").strip() != "206":
        return fetch(url, dest, min_bytes=min_bytes)
    step = -(-size // n)
    parts = [Path(f"{dest}.r{i}") for i in range(n)]

    def get(i):
        a, b = i * step, min(size, (i + 1) * step) - 1
        expect = b - a + 1
        for _ in range(8):
            have = parts[i].stat().st_size if parts[i].exists() else 0
            if have >= expect:
                break
            left = int(t_end - time.time()) if t_end else None
            if left is not None and left <= 0:
                break
            sh(f"curl -fsSL --retry 5 --retry-all-errors --speed-limit 100000 --speed-time 120 "
               + (f"--max-time {left} " if left else "") + f"-r {a + have}-{b} '{url}' >> '{parts[i]}'",
               name=f"range{i}_{dest.name}", check=False)
        if (parts[i].stat().st_size if parts[i].exists() else 0) != expect:
            raise RuntimeError(f"range {i} of {dest.name} incomplete")

    with ThreadPoolExecutor(n) as ex:
        list(ex.map(get, range(n)))
    sh(f"cat {' '.join(repr(str(p)) for p in parts)} > '{dest}.part' && mv '{dest}.part' '{dest}'", name=f"join_{dest.name}")
    for p in parts:
        p.unlink(missing_ok=True)
    if dest.stat().st_size != size:
        raise RuntimeError(f"{dest.name}: size {dest.stat().st_size} != {size}")
    return dest


def fetch_first(urls, dest, **kw):
    """Try mirrors in order."""
    last = None
    for u in urls:
        try:
            return fetch(u, dest, **kw)
        except Exception as e:      # noqa: BLE001 - try the next mirror
            last = e
            Path(str(dest) + ".part").unlink(missing_ok=True)
    raise RuntimeError(f"all mirrors failed for {Path(dest).name}: {last}")


def free_gb(p):
    return shutil.disk_usage(p).free / 1e9


def elapsed_h():
    return (time.time() - T0) / 3600.0


# %%
# ---- platform probe (also what a secret-less API push prints, so the environment is on record) ---------------------
def _ver(mod):
    try:
        return __import__(mod).__version__
    except Exception as e:      # noqa: BLE001
        return f"missing ({type(e).__name__})"


def _first_line(cmd):
    try:
        out = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, timeout=60)
        return (out.stdout + out.stderr).strip().splitlines()[0][:200]
    except Exception as e:      # noqa: BLE001
        return f"n/a ({type(e).__name__})"


PROBE = {
    "platform": platform.platform(),
    "python": sys.version.split()[0],
    "kaggle_image_build": os.environ.get("BUILD_DATE", "n/a"),
    "cpus": NCPU,
    "ram_gb": round(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1e9, 1) if hasattr(os, "sysconf") else None,
    "scratch_free_gb": round(free_gb(WORK), 1),
    "output_free_gb": round(free_gb(OUT), 1),
    "gpus": subprocess.run(["bash", "-c", "nvidia-smi -L 2>/dev/null || true"], capture_output=True, text=True).stdout.strip().splitlines(),
    "java": _first_line("java -version 2>&1"),
}
N_GPU = len(PROBE["gpus"])
for k, v in PROBE.items():
    log(f"{k}: {v}")
if ON_KAGGLE and N_GPU < 2:
    log(f"NOTE: {N_GPU} GPU(s) visible — the deep-learning modules run on the {NCPU} CPUs instead (ESM-2 in a child process; "
        "fastest on GPU T4 x2 or a large CPU host such as the TPU v5e-8 VM). The core ranking does not need a GPU.")


# %%
# ---- access token: Kaggle Secret `HF_TOKEN` (Add-ons → Secrets), or an environment variable ------------------------
def get_hf_token():
    for name in ("HF_TOKEN", "HUGGING_FACE_HUB_TOKEN", "HUGGINGFACE_TOKEN"):
        if os.environ.get(name):
            return os.environ[name].strip()
    try:
        from kaggle_secrets import UserSecretsClient
        client = UserSecretsClient()
    except Exception:           # noqa: BLE001 - not on Kaggle
        return None
    for name in ("HF_TOKEN", "HUGGING_FACE_HUB_TOKEN", "HUGGINGFACE_TOKEN", "hf_token"):
        try:
            tok = client.get_secret(name)
            if tok:
                return tok.strip()
        except Exception:       # noqa: BLE001 - secret not attached under this label
            continue
    return None


HF_TOKEN = get_hf_token()
DEMO = HF_TOKEN is None
if DEMO:
    print("\n" + "=" * 110)
    print(" No Hugging Face token found (attach the Kaggle Secret `HF_TOKEN`: Add-ons → Secrets → tick the box).")
    print(" Running in DEMO mode on a PUBLIC genome (GIAB HG002) with a spiked-in recessive case, to verify that every")
    print(" stage of the pipeline works on this machine. No patient data is touched in DEMO mode.")
    print("=" * 110 + "\n")
else:
    os.environ["HF_TOKEN"] = HF_TOKEN
    log("HF token found — running on the gated challenge dataset.")
S["mode"] = "DEMO (public HG002 + spike-in)" if DEMO else "CHALLENGE (gated dataset)"

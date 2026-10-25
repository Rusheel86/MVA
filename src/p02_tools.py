# %% [markdown]
# ## 2 · Command-line tools and plotting style
#
# One pinned micromamba environment provides the bioinformatics binaries (bcftools / samtools / htslib /
# minimap2 / BBTools), with an `apt` fallback. Its `bin/` is *appended* to `PATH` so the image's own Python is
# never shadowed.

# %%
ENV_DIR = WORK / "env"
MICROMAMBA = WORK / "bin" / "micromamba"
BIO_PKGS = "bcftools samtools htslib minimap2 bbmap"


def _have(tool):
    return shutil.which(tool) is not None


with stage("install tools", required=True):
    if not (_have("bcftools") and _have("samtools") and _have("tabix") and _have("minimap2")):
        try:
            fetch_first(["https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64"],
                        MICROMAMBA, min_bytes=1_000_000)
            MICROMAMBA.chmod(0o755)
            if not (ENV_DIR / "bin" / "bcftools").exists():
                sh(f"'{MICROMAMBA}' create -y -q -p '{ENV_DIR}' --strict-channel-priority -c conda-forge -c bioconda {BIO_PKGS}",
                   name="micromamba_create", timeout=1800)
            os.environ["PATH"] = os.environ["PATH"] + os.pathsep + str(ENV_DIR / "bin")
        except Exception as e:      # noqa: BLE001
            log(f"micromamba route failed ({str(e)[:200]}); falling back to apt")
            sh("apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq bcftools samtools tabix minimap2 bbmap",
               name="apt_install", timeout=1800)
    for tool in ("bcftools", "samtools", "tabix", "bgzip", "minimap2"):
        if not _have(tool):
            raise RuntimeError(f"required tool missing after installation: {tool}")
    S["tool_versions"] = {
        "bcftools": _first_line("bcftools --version"), "samtools": _first_line("samtools --version"),
        "minimap2": _first_line("minimap2 --version"), "bbduk": _first_line("bbduk.sh --version 2>&1 | grep -i 'version' || echo missing"),
    }
    for k, v in S["tool_versions"].items():
        log(f"{k}: {v}")
    # small pure-Python dependencies (wheels exist for the Kaggle image's Python)
    sh(f"'{sys.executable}' -m pip install -q --no-input pysam pyfaidx 2>&1 | tail -2", name="pip_install", check=False, timeout=900)

# %%
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from IPython.display import Image, Markdown, display

# One chart style for every figure: recessive axes, thin marks, a colour-blind-validated categorical order
# (blue, orange, aqua, yellow); status colours are reserved for status.
VIZ = dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781", grid="#e1e0d9", axis="#c3c2b7",
           s1="#2a78d6", s2="#eb6834", s3="#1baf7a", s4="#eda100", neutral="#c3c2b7", critical="#d03b3b")
plt.rcParams.update({
    "figure.facecolor": VIZ["surface"], "axes.facecolor": VIZ["surface"], "savefig.facecolor": VIZ["surface"],
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": VIZ["axis"], "axes.labelcolor": VIZ["ink2"],
    "axes.titlesize": 10.5, "axes.titleweight": "bold", "axes.titlecolor": VIZ["ink"], "axes.titlelocation": "left",
    "xtick.color": VIZ["muted"], "ytick.color": VIZ["muted"], "xtick.labelcolor": VIZ["ink2"], "ytick.labelcolor": VIZ["ink2"],
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": VIZ["grid"],
    "grid.linewidth": 0.6, "grid.linestyle": "-", "axes.axisbelow": True, "legend.frameon": False,
    "legend.fontsize": 8.5, "lines.linewidth": 2.0, "figure.dpi": 110, "savefig.dpi": 160, "savefig.bbox": "tight",
})


def save_fig(fig, name, show=True):
    """Save a figure under OUT/figures and display it."""
    path = OUT / "figures" / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    S.setdefault("figures", {})[name] = str(path.relative_to(OUT))
    if show:
        display(Image(filename=str(path)))
    return path


def md(text):
    display(Markdown(text))


def df_md(df, max_rows=30, floatfmt="{:.3g}"):
    """Small DataFrame -> GitHub markdown table (no external dependency)."""
    df = df.head(max_rows)
    cols = [str(c) for c in df.columns]

    def fmt(v):
        if isinstance(v, float):
            return "" if math.isnan(v) else floatfmt.format(v)
        return str(v).replace("|", "\\|").replace("\n", " ")

    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for row in df.itertuples(index=False):
        lines.append("| " + " | ".join(fmt(v) for v in row) + " |")
    return "\n".join(lines)

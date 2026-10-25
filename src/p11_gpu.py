# %% [markdown]
# ## 11 · Deep-learning evidence on the GPUs
#
# Three orthogonal models are applied to the leading genes only (a few dozen alleles), each on its own device:
#
# * **SpliceAI** (Jaganathan 2019; five-model ensemble, 10-kb context, ±500 bp scan) — scores every rare non-coding
#   and synonymous allele in the leading genes, so that a deep-intronic or cryptic-splice second hit is not missed.
#   Scores ≥ 0.2 enter the evidence model; alleles with ≤ 0.1 keep their benign-leaning default.
# * **ESM-2** (650 M-parameter protein language model; masked-marginal log-likelihood ratios) — an in-silico saturation
#   mutagenesis of the leading protein. The candidate missense change is placed on the distribution of all 19·L
#   substitutions, and against ClinVar-pathogenic and population-observed missense variants of the same gene.
#   It is reported as orthogonal support and is not added to the ACMG points (one calibrated predictor only).
# * **AlphaFold** structural context — predicted-structure confidence, burial and 3D distance from the candidate residue
#   to known pathogenic missense positions.

# %%
def ensure_modules(mods):
    """Install missing Python modules at run time ({import name: pip arguments}); used on images without them
    (e.g. a TPU host has no CUDA builds and may lack TensorFlow or transformers)."""
    import importlib.util
    for name, spec in mods.items():
        if importlib.util.find_spec(name) is None:
            log(f"installing {spec} (module {name} not in this image)")
            sh(f"'{sys.executable}' -m pip install -q --no-input {spec}", name=f"pip_{name}", timeout=1200)
            importlib.invalidate_caches()


def load_spliceai_models():
    """The five published SpliceAI Keras models (weights installed at run time from PyPI; CC BY-NC 4.0)."""
    import importlib.util
    ensure_modules({"tensorflow": "tensorflow-cpu tf-keras"})
    if importlib.util.find_spec("spliceai") is None:
        sh(f"'{sys.executable}' -m pip install -q --no-input --no-deps spliceai", name="pip_spliceai", timeout=600)
        importlib.invalidate_caches()
    base = Path(importlib.util.find_spec("spliceai").submodule_search_locations[0]) / "models"
    paths = [base / f"spliceai{i}.h5" for i in range(1, 6)]
    import tensorflow as tf
    gpus = tf.config.list_physical_devices("GPU")
    if len(gpus) > 1:                       # leave GPU 0 to PyTorch
        tf.config.set_visible_devices(gpus[-1], "GPU")
    for g in tf.config.get_visible_devices("GPU"):
        try:
            tf.config.experimental.set_memory_growth(g, True)
        except Exception:       # noqa: BLE001
            pass
    try:
        import keras
        models = [keras.models.load_model(str(p), compile=False) for p in paths]
    except Exception as e:      # noqa: BLE001 - legacy H5 under Keras 3: fall back to tf_keras
        log(f"keras loader failed ({type(e).__name__}); using tf_keras")
        import tf_keras
        models = [tf_keras.models.load_model(str(p), compile=False) for p in paths]
    return models, len(tf.config.get_visible_devices("GPU"))


_OH = np.zeros((256, 4), dtype=np.float32)
for _i, _b in enumerate(b"ACGT"):
    _OH[_b, _i] = 1.0


def one_hot(seq):
    return _OH[np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)]


def spliceai_delta(models, fasta, chrom, pos, ref, alt, strand, dist=500):
    """Raw (unmasked) SpliceAI delta scores for one allele, following the reference implementation."""
    cov = 2 * dist + 1
    wid = 10000 + cov
    start = pos - wid // 2 - 1
    if start < 0 or len(ref) > dist or len(alt) > dist:
        return None
    seq = str(fasta[chrom][start:start + wid]).upper()
    if len(seq) != wid or seq[wid // 2:wid // 2 + len(ref)] != ref.upper():
        return None
    x_ref, x_alt = seq, seq[:wid // 2] + alt.upper() + seq[wid // 2 + len(ref):]
    xr, xa = one_hot(x_ref)[None], one_hot(x_alt)[None]
    if strand == "-":
        xr, xa = xr[:, ::-1, ::-1], xa[:, ::-1, ::-1]
    if xr.shape == xa.shape:            # SNV / MNV: reference and alternative in one batch
        y = np.mean([np.asarray(m(np.concatenate([xr, xa]), training=False)) for m in models], axis=0)
        yr, ya = y[:1], y[1:]
    else:
        yr = np.mean([np.asarray(m(xr, training=False)) for m in models], axis=0)
        ya = np.mean([np.asarray(m(xa, training=False)) for m in models], axis=0)
    if strand == "-":
        yr, ya = yr[:, ::-1], ya[:, ::-1]
    rl, al = len(ref), len(alt)
    if rl > 1 and al == 1:          # deletion: pad the alternative track
        ya = np.concatenate([ya[:, :cov // 2 + al], np.zeros((1, rl - al, 3), dtype=ya.dtype), ya[:, cov // 2 + al:]], axis=1)
    elif rl == 1 and al > 1:        # insertion: collapse the inserted positions
        ya = np.concatenate([ya[:, :cov // 2], ya[:, cov // 2:cov // 2 + al].max(axis=1, keepdims=True), ya[:, cov // 2 + al:]], axis=1)
    if yr.shape != ya.shape:
        return None
    d_acc, d_don = ya[0, :, 1] - yr[0, :, 1], ya[0, :, 2] - yr[0, :, 2]
    res = dict(DS_AG=float(d_acc.max()), DS_AL=float((-d_acc).max()), DS_DG=float(d_don.max()), DS_DL=float((-d_don).max()),
               DP_AG=int(d_acc.argmax() - cov // 2), DP_AL=int((-d_acc).argmax() - cov // 2),
               DP_DG=int(d_don.argmax() - cov // 2), DP_DL=int((-d_don).argmax() - cov // 2))
    res["max"] = max(res["DS_AG"], res["DS_AL"], res["DS_DG"], res["DS_DL"])
    return res


with stage("SpliceAI on the leading genes (GPU)"):
    if not CFG["RUN_GPU_MODELS"]:
        raise RuntimeError("disabled in CFG")
    from pyfaidx import Fasta
    todo = {}
    for v in RECS:
        if v["gene"] in TOP_GENES and v["gene"] in GENE_POS and v["csq"] not in LOF_CSQ:
            todo[(vkey(v), v["gene"])] = v
    todo = sorted(todo.values(), key=lambda v: (TOP_GENES.index(v['gene']), -v['p_path']))[:200]
    models, n_tf_gpu = load_spliceai_models()
    fa = Fasta(RES["ref"], as_raw=True, sequence_always_upper=True)
    t0, n_scored, n_new = time.time(), 0, 0
    for v in todo:
        r = spliceai_delta(models, fa, v["chrom"], int(v["pos"]), v["ref"], v["alt"], GENE_POS[v["gene"]][3])
        if r is None:
            continue
        n_scored += 1
        v["spliceai_custom"] = round(r["max"], 3)
        v["spliceai_custom_detail"] = r
        prev = _num(v.get("spliceai"))
        if prev is None or r["max"] > prev:         # wider scan window than the pre-computed scores: keep the larger delta
            v["spliceai"] = round(r["max"], 3)
            n_new += 1
    S["spliceai"] = dict(n_alleles=n_scored, seconds=round(time.time() - t0, 1), device="GPU" if n_tf_gpu else "CPU",
                         n_ge_0p2=int(sum(1 for v in todo if (v.get("spliceai_custom") or 0) >= 0.2)))
    log(f"SpliceAI: {n_scored} alleles in {TOP_GENES} scored on {S['spliceai']['device']} in {S['spliceai']['seconds']}s; "
        f"{S['spliceai']['n_ge_0p2']} with delta >= 0.2")
    del models
    rescore(RECS)
    HYPS = build_hypotheses(RECS, GENES_DX, PHENO_TOP, W)
    S["ranking_stage"] = "after deep annotation and SpliceAI"
    SUB_ROWS, parsed = write_submission(HYPS)
    if sorted(HYPS[0]["key"]) != S.get("top_key"):
        log("the leading hypothesis changed after SpliceAI: rebuilding its evidence tables")
        analyse_leader()
    S["top_gene"] = HYPS[0]["gene"]
    md("**Top hypotheses after SpliceAI**\n\n" + df_md(public_table(hyp_table(HYPS, 6))[["rank", "EPCR", "gene", "model", "variants", "evidence"]]))

# %%
AA20 = "ACDEFGHIKLMNPQRSTVWY"


def esm_landscape(seq, model_name="facebook/esm2_t33_650M_UR50D", batch=8, devices=None, budget_s=None):
    """Masked-marginal log-likelihood ratios for every substitution: array [L, 20] in AA20 order (wild type = 0).
    One model replica per GPU, positions split between them; or one float32 replica on the CPU of a large host, which
    gives up (TimeoutError) when the first batches project a run longer than `budget_s`."""
    import torch
    from transformers import AutoTokenizer, EsmForMaskedLM
    tok = AutoTokenizer.from_pretrained(model_name)
    devices = devices or [f"cuda:{i}" for i in range(torch.cuda.device_count())]
    if devices == ["cpu"]:
        torch.set_num_threads(NCPU)
    ids = tok(seq, return_tensors="pt")["input_ids"][0]         # <cls> seq <eos>
    if len(ids) != len(seq) + 2:
        raise RuntimeError("ESM tokenisation does not map one token per residue")
    aa_ids = torch.tensor([tok.convert_tokens_to_ids(a) for a in AA20])
    L = len(seq)
    out = np.zeros((L, 20), dtype=np.float32)
    replicas = [EsmForMaskedLM.from_pretrained(model_name).to(d).eval() for d in devices]
    replicas = [r.half() if d.startswith("cuda") else r for r, d in zip(replicas, devices)]
    shards = [list(range(L))[i::len(devices)] for i in range(len(devices))]
    errors = []

    def worker(m, dev, pos_list):
        try:
            t_start = time.time()
            with torch.no_grad():
                for i in range(0, len(pos_list), batch):
                    if budget_s and i == 2 * batch:
                        projected = (time.time() - t_start) / i * len(pos_list)
                        if projected > budget_s:
                            raise TimeoutError(f"projected {projected / 60:.0f} min on {dev} exceeds the {budget_s / 60:.0f}-min budget")
                    chunk = pos_list[i:i + batch]
                    x = ids.repeat(len(chunk), 1).to(dev)
                    for j, p in enumerate(chunk):
                        x[j, p + 1] = tok.mask_token_id
                    logits = m(input_ids=x).logits.float()
                    for j, p in enumerate(chunk):
                        lp = torch.log_softmax(logits[j, p + 1], dim=-1)[aa_ids.to(dev)].cpu().numpy()
                        out[p] = lp - lp[AA20.index(seq[p])] if seq[p] in AA20 else lp
        except Exception as e:      # noqa: BLE001
            errors.append(e)

    threads = [threading.Thread(target=worker, args=(m, d, sh_)) for m, d, sh_ in zip(replicas, devices, shards)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    del replicas
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    if errors:
        raise errors[0]
    return out, devices


_ESM_CPU = {}
ESM_JOB = r"""
import json, sys, time
import numpy as np
import torch
from transformers import AutoTokenizer, EsmForMaskedLM
cfg = json.load(open(sys.argv[1]))
AA20 = "ACDEFGHIKLMNPQRSTVWY"
torch.set_num_threads(cfg["threads"])
tok = AutoTokenizer.from_pretrained(cfg["model"])
m = EsmForMaskedLM.from_pretrained(cfg["model"]).eval()
seq = cfg["seq"]
ids = tok(seq, return_tensors="pt")["input_ids"][0]
assert len(ids) == len(seq) + 2, "tokenisation is not one token per residue"
aa_ids = torch.tensor([tok.convert_tokens_to_ids(a) for a in AA20])


def score(positions, budget=None):
    res, t0, B = {}, time.time(), cfg["batch"]
    with torch.no_grad():
        for i in range(0, len(positions), B):
            if budget and i == 2 * B and (time.time() - t0) / i * len(positions) > budget:
                return None, (time.time() - t0) / i * len(positions)
            chunk = positions[i:i + B]
            x = ids.repeat(len(chunk), 1)
            for j, p in enumerate(chunk):
                x[j, p + 1] = tok.mask_token_id
            lg = m(input_ids=x).logits.float()
            for j, p in enumerate(chunk):
                lp = torch.log_softmax(lg[j, p + 1], -1)[aa_ids].numpy()
                res[p] = lp - lp[AA20.index(seq[p])] if seq[p] in AA20 else lp
    return res, None


full, proj = score(list(range(len(seq))), cfg["budget_s"]) if cfg["full"] else (None, None)
if full is not None:
    np.savez(sys.argv[2], llr=np.stack([full[p] for p in range(len(seq))]), full=1)
else:
    rows, _ = score(sorted(set(cfg["positions"])))
    np.savez(sys.argv[2], pos=np.array(sorted(rows)), rows=np.stack([rows[p] for p in sorted(rows)]), full=0,
             projected_s=proj if proj else -1.0)
"""


def esm_subprocess(wseq, positions0, full, budget_s=None, model_name="facebook/esm2_t33_650M_UR50D"):
    """ESM-2 in a child process on the CPU: the full landscape (within `budget_s`, else only `positions0`) or only
    `positions0` (0-based window indices). Returns dict(llr=[L,20] or None, rows={pos0: [20]}, projected_s)."""
    job, cfg, out = WORK / "esm_job.py", WORK / "esm_job.json", WORK / "esm_job.npz"
    job.write_text(ESM_JOB)
    cfg.write_text(json.dumps(dict(model=model_name, seq=wseq, positions=[int(p) for p in positions0], full=bool(full),
                                   budget_s=budget_s, threads=NCPU, batch=8)))
    out.unlink(missing_ok=True)
    env = {**os.environ, "TORCH_DEVICE_BACKEND_AUTOLOAD": "0", "PJRT_DEVICE": "CPU", "CUDA_VISIBLE_DEVICES": ""}
    sh(f"'{sys.executable}' '{job}' '{cfg}' '{out}'", name="esm_subprocess", env=env,
       timeout=int((budget_s or 0) + 3600))
    z = np.load(out)
    if int(z["full"]):
        return dict(llr=z["llr"], rows={}, projected_s=None)
    return dict(llr=None, rows={int(p): r for p, r in zip(z["pos"], z["rows"])}, projected_s=float(z["projected_s"]))


def esm_single(seq, pos0, model_name="facebook/esm2_t33_650M_UR50D"):
    """CPU fallback: one masked position only (the model is loaded once and cached; child process without a GPU)."""
    if N_GPU == 0:
        return esm_subprocess(seq, [pos0], full=False, model_name=model_name)["rows"][pos0]
    import torch
    from transformers import AutoTokenizer, EsmForMaskedLM
    if model_name not in _ESM_CPU:
        _ESM_CPU[model_name] = (AutoTokenizer.from_pretrained(model_name), EsmForMaskedLM.from_pretrained(model_name).eval())
    tok, m = _ESM_CPU[model_name]
    ids = tok(seq, return_tensors="pt")["input_ids"]
    ids[0, pos0 + 1] = tok.mask_token_id
    with torch.no_grad():
        lp = torch.log_softmax(m(input_ids=ids).logits[0, pos0 + 1].float(), dim=-1)
    lp = lp[[tok.convert_tokens_to_ids(a) for a in AA20]].numpy()
    return lp - lp[AA20.index(seq[pos0])]


with stage("ESM-2 missense landscape of the leading protein (GPU)"):
    if not CFG["RUN_GPU_MODELS"]:
        raise RuntimeError("disabled in CFG")
    TOP = HYPS[0]
    mis = [v for v in TOP["variants"] if v["csq"] == "missense"]
    if not mis:
        raise RuntimeError("the leading hypothesis has no missense allele")
    v = mis[0]
    tx = TX[GENE_TX[v["gene"]][0]]
    seq = PROT.get(tx["protein"], "").rstrip("*")
    m = re.match(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})", v.get("hgvsp") or "")
    inv = {b: a for a, b in AA3.items()}
    if not (seq and m):
        raise RuntimeError("protein sequence or protein change unavailable")
    wt, pos1, mt = inv.get(m.group(1), "?"), int(m.group(2)), inv.get(m.group(3), "?")
    if mt not in AA20 or pos1 > len(seq) or seq[pos1 - 1] != wt:
        raise RuntimeError("protein sequence does not match the annotated protein change")
    ensure_modules({"torch": "torch --index-url https://download.pytorch.org/whl/cpu", "transformers": "transformers"})
    t0 = time.time()
    # ESM-2 was trained on <= 1,022 residues: longer proteins are scored in a 1,022-residue window around the change
    off = 0 if len(seq) <= 1022 else max(0, min(pos1 - 1 - 511, len(seq) - 1022))
    wseq = seq[off:off + 1022]
    LLR, cpu_rows = None, {}
    if N_GPU > 0:
        LLR, devs = esm_landscape(wseq)
    else:
        need = sorted({int(_num(x["aa_pos"])) - 1 - off for x in TOP["variants"] if _num(x.get("aa_pos")) is not None
                       and 0 <= int(_num(x["aa_pos"])) - 1 - off < len(wseq)} | {pos1 - 1 - off})
        r = esm_subprocess(wseq, need, full=NCPU >= 32, budget_s=CFG.get("ESM_CPU_BUDGET_MIN", 40) * 60)
        LLR, devs = r["llr"], [f"cpu x{NCPU} (child process)"]
        cpu_rows = {p0 + 1 + off: row for p0, row in r["rows"].items()}
        if LLR is None and NCPU >= 32:
            log(f"ESM-2 full landscape skipped on CPU (projected {r['projected_s'] / 60:.0f} min); scored positions {sorted(cpu_rows)}")
    if LLR is not None:
        cand = float(LLR[pos1 - 1 - off, AA20.index(mt)])
        flat = LLR.ravel()
        flat = flat[flat != 0]
        ESM_CTX = dict(gene=v["gene"], seq=seq, off=off, wseq=wseq, llr=LLR)       # reused by the allele-level triage (Track 2)
        S["esm"] = dict(model="ESM-2 650M (esm2_t33_650M_UR50D)", protein=tx["protein"], length=len(seq),
                        window=f"{off + 1}-{off + len(wseq)}", change=f"{wt}{pos1}{mt}",
                        llr=round(cand, 2), percentile_more_damaging_than=round(float((flat > cand).mean()) * 100, 1),
                        devices=devs, seconds=round(time.time() - t0, 1))
    else:
        row = cpu_rows.get(pos1)
        if row is None:
            row = esm_single(wseq, pos1 - 1 - off)
        cand = float(row[AA20.index(mt)])
        ESM_CTX = dict(gene=v["gene"], seq=seq, off=off, wseq=wseq, llr=None, rows={**cpu_rows, pos1: row})  # CPU: single positions
        S["esm"] = dict(model="ESM-2 650M (single position, CPU)", protein=tx["protein"], length=len(seq), change=f"{wt}{pos1}{mt}",
                        llr=round(cand, 2), seconds=round(time.time() - t0, 1))
    log(f"ESM-2: {S['esm']}")
    if LLR is not None:
        # reference sets from public data: ClinVar P/LP missense in this gene and population-observed missense (gnomAD)
        chrom, gs, ge = GENE_POS[v["gene"]][0], GENE_POS[v["gene"]][1], GENE_POS[v["gene"]][2]
        am_map = am_region(chrom, gs, ge)

        def llr_of(pv):
            mm = re.match(r"([A-Z])(\d+)([A-Z])$", pv)
            if not mm or not (off < int(mm.group(2)) <= off + len(wseq)) or seq[int(mm.group(2)) - 1] != mm.group(1) or mm.group(3) not in AA20:
                return None
            return float(LLR[int(mm.group(2)) - 1 - off, AA20.index(mm.group(3))])

        cvg = CLINVAR_PLP[(CLINVAR_PLP.gene == v["gene"]) & CLINVAR_PLP.is_missense]
        plp = [llr_of(am_map[(r.pos, r.ref, r.alt)][1]) for r in cvg.itertuples() if (r.pos, r.ref, r.alt) in am_map]
        plp = [x for x in plp if x is not None]
        S["esm"]["clinvar_plp_missense_n"] = len(plp)
        S["esm"]["clinvar_plp_median_llr"] = round(float(np.median(plp)), 2) if plp else None
        fig, ax = plt.subplots(figsize=(7.2, 3.0))
        ax.hist(flat, bins=60, color=VIZ["neutral"], label=f"all {len(flat):,} possible substitutions")
        ymax = ax.get_ylim()[1]
        if plp:
            ax.scatter(plp, np.full(len(plp), ymax * 0.06), s=26, color=VIZ["s2"], zorder=3, linewidths=0,
                       label=f"ClinVar P/LP missense in {v['gene']} (n={len(plp)})")
        ax.axvline(cand, color=VIZ["s1"], linewidth=2, label=f"candidate {wt}{pos1}{mt} (LLR {cand:.1f})")
        ax.set_ylim(0, ymax * 1.05)
        ax.set_xlabel("ESM-2 masked-marginal log-likelihood ratio (more negative = less tolerated)")
        ax.set_ylabel("substitutions")
        ax.set_title(f"{v['gene']} in-silico saturation mutagenesis (ESM-2): where the candidate change falls")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3)
        save_fig(fig, "esm_landscape")

# %%
with stage("AlphaFold structural context of the missense allele"):
    TOP = HYPS[0]
    mis = [v for v in TOP["variants"] if v["csq"] == "missense" and _num(v.get("aa_pos")) is not None]
    if not mis:
        raise RuntimeError("the leading hypothesis has no missense allele")
    v = mis[0]
    res_no = int(_num(v["aa_pos"]))
    u = http_json("GET", "https://rest.uniprot.org/uniprotkb/search",
                  params={"query": f"gene_exact:{v['gene']} AND organism_id:9606 AND reviewed:true", "format": "json",
                          "fields": "accession,ft_domain,ft_region,ft_variant,length"})["results"][0]
    acc = u["primaryAccession"]
    feats = [dict(type=f["type"], name=f.get("description", ""), start=f["location"]["start"]["value"], end=f["location"]["end"]["value"])
             for f in u.get("features", []) if f["type"] in ("Domain", "Region")]
    here = [f for f in feats if f["start"] and f["end"] and f["start"] <= res_no <= f["end"]]
    pdb_url = http_json("GET", f"https://alphafold.ebi.ac.uk/api/prediction/{acc}")[0]["pdbUrl"]
    pdb = HTTP.get(pdb_url, timeout=120).text
    ca, cb = {}, {}
    for line in pdb.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() in ("CA", "CB"):
            r, xyz_ = int(line[22:26]), np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if line[12:16].strip() == "CA":
                ca[r] = (xyz_, float(line[60:66]))
                cb.setdefault(r, xyz_)          # glycine: CA stands in for CB
            else:
                cb[r] = xyz_
    if res_no not in ca:
        raise RuntimeError("residue not in the AlphaFold model")
    xyz, plddt = ca[res_no]
    dists = {r: float(np.linalg.norm(c - xyz)) for r, (c, _) in ca.items()}
    chrom, gs, ge = GENE_POS[v["gene"]][0], GENE_POS[v["gene"]][1], GENE_POS[v["gene"]][2]
    am_map = am_region(chrom, gs, ge)
    cvg = CLINVAR_PLP[(CLINVAR_PLP.gene == v["gene"]) & CLINVAR_PLP.is_missense]
    plp_res = sorted({int(re.sub(r"\D", "", am_map[(r.pos, r.ref, r.alt)][1])) for r in cvg.itertuples() if (r.pos, r.ref, r.alt) in am_map} - {res_no})
    near = sorted(((r, round(dists[r], 1)) for r in plp_res if r in dists), key=lambda x: x[1])[:5]
    # disease-annotated missense variants curated by UniProt (literature-reported patient alleles, e.g. "in MVA1; ...")
    dis = []
    for f in u.get("features", []):
        alt_ = (f.get("alternativeSequence") or {})
        d_ = f.get("description", "")
        m_ = re.match(r"in ([A-Z][A-Z0-9-]{1,});?", d_)
        one_aa = len(alt_.get("alternativeSequences") or []) == 1 and f["location"]["start"]["value"] == f["location"]["end"]["value"]
        somatic = re.search(r"cancer|tumou?r|carcinoma|cell line|somatic", d_, re.I)
        if f["type"] == "Natural variant" and m_ and one_aa and not somatic:
            r = int(f["location"]["start"]["value"])
            if r != res_no and r in dists:
                dis.append(dict(residue=r, change=f"{alt_.get('originalSequence')}{r}{alt_['alternativeSequences'][0]}",
                                disease=m_.group(1), distance_A=round(dists[r], 1),
                                same_domain=any(f2["type"] == "Domain" and f2["start"] <= r <= f2["end"] for f2 in here)))
    dis.sort(key=lambda x: x["distance_A"])

    cb_ids = np.array(sorted(cb))
    cb_xyz = np.stack([cb[r] for r in cb_ids])
    n_within = ((np.linalg.norm(cb_xyz[:, None] - cb_xyz[None], axis=-1) <= 10.0).sum(1) - 1)
    BURIAL = dict(zip(cb_ids.tolist(), n_within.tolist()))      # CB neighbours within 10 Å (a burial proxy)

    def burial(r):
        return BURIAL[r]
    dom_res = [r for f in here if f["type"] == "Domain" for r in range(f["start"], f["end"] + 1) if r in cb]
    mane_len = len(PROT.get(TX[GENE_TX[v["gene"]][0]]["protein"], "").rstrip("*")) if v["gene"] in GENE_TX else 0
    AF_CTX = dict(gene=v["gene"], acc=acc, ca=ca, cb=cb, burial=BURIAL, feats=feats, disease_missense=dis,      # reused by Track 2
                  numbering_matches_mane=bool(mane_len) and max(ca) == mane_len)
    if not AF_CTX["numbering_matches_mane"]:
        log(f"AlphaFold/UniProt length {max(ca)} differs from the MANE protein ({mane_len}): residue-level reuse in Track 2 disabled")
    S["structure"] = dict(uniprot=acc, model=pdb_url.split("/")[-1], residue=res_no, plddt=round(plddt, 1),
                          neighbours_within_10A=int(sum(1 for d in dists.values() if 0 < d <= 10.0)),
                          domain="; ".join(f"{f['name']} ({f['start']}-{f['end']})" for f in here) or "none annotated",
                          nearest_clinvar_plp_missense=[dict(residue=r, distance_A=d) for r, d in near],
                          disease_missense_uniprot=dis[:6], n_disease_missense_uniprot=len(dis),
                          burial_cb_10A=burial(res_no),
                          burial_percentile_in_domain=(round(float(np.mean([burial(r) < burial(res_no) for r in dom_res])) * 100)
                                                       if dom_res else None))
    log(f"AlphaFold context: {S['structure']}")

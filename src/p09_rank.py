# %% [markdown]
# ## 9 · Ranking
#
# 1. **Phenotype scoring** — HPO-only ("blind") and HPO + stated indication.
# 2. **First-pass ranking** of every homozygous, compound-heterozygous and heterozygous hypothesis genome-wide.
# 3. **Second-hit rescue** — for the genes that already carry one credible allele, every rare non-coding variant in
#    the gene and ±50 kb is added as a potential second allele (deep-intronic, UTR, regulatory), so a hypomorphic
#    non-coding partner is not missed. They start with a small penalty until the splice model has scored them.
# 4. **Robustness** — the ranking is recomputed under random perturbation of every heuristic weight, with each HPO
#    term left out in turn, and with ClinVar or the missense predictor switched off.
# 5. **Submission file** — written immediately (and rewritten after each later evidence stage), parsed with the
#    organisers' own scoring code, and every submitted allele is checked against the source VCF.

# %%
with stage("phenotype scoring", required=True):
    HPO_TERMS = [h["term"] for h in S["hpo"] if h.get("weight", 1.0) > 0]
    HPO_WEIGHTS = {h["term"]: h.get("weight", 1.0) for h in S["hpo"]}
    GENES_BLIND, DIS_BLIND = PE.score_genes(HPO_TERMS, HPO_WEIGHTS, diagnosis=None)
    GENES_DX, _ = PE.score_genes(HPO_TERMS, HPO_WEIGHTS, diagnosis=S.get("diagnosis") or None)
    PHENO_TOP = max(v["pheno"] for v in GENES_BLIND.values()) or 1.0
    _order = sorted(GENES_BLIND, key=lambda g: -GENES_BLIND[g]["pheno"])
    GENE_RANK_BLIND = {g: i + 1 for i, g in enumerate(_order)}
    S["pheno"] = dict(
        n_diseases=PE.n_diseases, n_genes=len(GENES_BLIND),
        dx_panel=sorted(g for g, v in GENES_DX.items() if v["dx_match"]),
        top_diseases=[dict(rank=i + 1, disease=d, name=PE.d_name[d], score=round(DIS_BLIND[d][0], 3),
                           genes=",".join(sorted(PE.d2g.get(d, []))[:5]), inheritance="/".join(sorted(PE.d_inh.get(d, []))))
                      for i, d in enumerate(sorted(DIS_BLIND, key=lambda d: -DIS_BLIND[d][0])[:12])])
    md("**Top diseases by HPO similarity (blind to the diagnosis)**\n\n" + df_md(pd.DataFrame(S["pheno"]["top_diseases"])))
    if S.get("diagnosis"):
        log(f"virtual panel from the stated indication '{S['diagnosis']}': {S['pheno']['dx_panel']}")


# %%
def second_hit_rows(recs, pool_df, gene_table, pheno_top, w=W):
    """Rare non-coding variants in (and around) genes that already carry one credible heterozygous allele."""
    best = {}
    for v in recs:
        if v["zyg"] != "het" or v["p_path"] < 0.5:
            continue
        s = v["p_path"] * pheno_factor(gene_table.get(v["gene"]) or {}, pheno_top, w)
        if s > best.get(v["gene"], 0):
            best[v["gene"]] = s
    genes = [g for g in sorted(best, key=lambda g: -best[g]) if best[g] >= 0.02][:CFG["TOP_GENES_DEEP"]]
    out = []
    for g in genes:
        if g not in GENE_POS:
            continue
        chrom, start, end, strand, _ = GENE_POS[g]
        sub = pool_df[(pool_df.gene == g) & (pool_df.zyg == "het")]
        flank = pool_df[(pool_df.gene == "") & (pool_df.chrom == chrom) & (pool_df.zyg == "het")
                        & (pool_df.pos >= start - CFG["FLANK_SECOND_HIT"]) & (pool_df.pos <= end + CFG["FLANK_SECOND_HIT"])].copy()
        if len(flank):
            flank["gene"] = g
            up = (flank.pos < start) if strand == "+" else (flank.pos > end)
            flank["csq"] = np.where(up, "upstream", "downstream")
        out.append(pd.concat([sub, flank]))
    return (pd.concat(out).reset_index(drop=True) if out else pool_df.iloc[0:0]), genes


def rank_case(cand_df, pool_df, gene_table, pheno_top, w=W, rescue=True):
    recs = score_variants(cand_df, gene_table, LOF_GENES)
    rescue_genes = []
    if rescue and pool_df is not None and len(pool_df):
        extra, rescue_genes = second_hit_rows(recs, pool_df, gene_table, pheno_top, w)
        if len(extra):
            have = {(vkey(v), v["gene"]) for v in recs}
            recs += [v for v in score_variants(extra, gene_table, LOF_GENES) if (vkey(v), v["gene"]) not in have]
    return recs, build_hypotheses(recs, gene_table, pheno_top, w), rescue_genes


with stage("first-pass ranking", required=True):
    RECS, HYPS, RESCUE_GENES = rank_case(CAND, POOL, GENES_DX, PHENO_TOP)
    _, HYPS_BLIND, _ = rank_case(CAND, POOL, GENES_BLIND, PHENO_TOP)
    S["ranking_stage"] = "first pass (offline annotations)"
    S["n_hypotheses"] = len(HYPS)
    top = HYPS[0]
    S["top_gene"] = top["gene"]
    md(f"**{len(HYPS):,} genotype hypotheses** from {len(RECS):,} scored variant rows. Top hypotheses (with the stated indication):\n\n"
       + df_md(public_table(hyp_table(HYPS, 12))))
    md("**Blind to the diagnosis (HPO terms only):**\n\n"
       + df_md(public_table(hyp_table(HYPS_BLIND, 5))[["rank", "EPCR", "gene", "model", "variants", "phenotype_fit"]]))
    S["top_gene_pheno_rank_blind"] = GENE_RANK_BLIND.get(top["gene"])
    S["top_same_blind"] = bool(HYPS_BLIND and HYPS_BLIND[0]["key"] == top["key"])
    log(f"top hypothesis: {top['gene']} {top['model']} EPCR={top['epcr']:.3f}; phenotype rank of the gene (HPO only): "
        f"{S['top_gene_pheno_rank_blind']} of {len(GENES_BLIND):,}; same top hypothesis without the diagnosis: {S['top_same_blind']}")


# %%
def robustness(recs, hyps, n_perturb=300, seed=11):
    """How stable is the top hypothesis? (a) every heuristic weight multiplied by a random factor in [0.5, 2];
    (b) each HPO term left out; (c) ClinVar ignored; (d) missense predictor ignored; (e) diagnosis ignored."""
    rng = np.random.default_rng(seed)
    top_key = hyps[0]["key"]
    out = {}
    same, ranks = 0, []
    for _ in range(n_perturb):
        w = {}
        for k, v in W.items():
            f = float(np.exp(rng.uniform(-math.log(2), math.log(2))))
            w[k] = v * f if k == "pheno_power" else min(1.0, v * f)
        h2 = build_hypotheses(recs, GENES_DX, PHENO_TOP, w)
        r = next((i + 1 for i, h in enumerate(h2) if h["key"] == top_key), None)
        ranks.append(r or 999)
        same += int(r == 1)
    out["weight_perturbation"] = dict(n=n_perturb, top1_unchanged=round(same / n_perturb, 3), worst_rank=int(max(ranks)))
    loo = []
    for h in S["hpo"]:
        if h.get("weight", 1.0) <= 0:
            continue
        terms = [t for t in HPO_TERMS if t != h["term"]]
        gt, _ = PE.score_genes(terms, HPO_WEIGHTS, diagnosis=None)        # blind, so that the term itself matters
        ptop = max(v["pheno"] for v in gt.values()) or 1.0
        _, h2, _ = rank_case(CAND, POOL, gt, ptop)
        loo.append(dict(term_left_out=h["term"], label=h.get("label", ""), rank_of_top=next((i + 1 for i, x in enumerate(h2) if x["key"] == top_key), None)))
    out["leave_one_term_out_blind"] = loo

    def ablate(mutator, table):
        c2 = CAND.copy()
        mutator(c2)
        _, h2, _ = rank_case(c2, POOL, table, PHENO_TOP)
        return next((i + 1 for i, x in enumerate(h2) if x["key"] == top_key), None)

    out["rank_without_clinvar"] = ablate(lambda c: c.__setitem__("clnsig", ""), GENES_DX)
    out["rank_without_missense_predictor"] = ablate(lambda c: c.__setitem__("am", np.nan), GENES_DX)
    out["rank_without_diagnosis"] = ablate(lambda c: None, GENES_BLIND)
    out["rank_without_clinvar_and_diagnosis"] = ablate(lambda c: c.__setitem__("clnsig", ""), GENES_BLIND)
    return out


with stage("robustness analysis"):
    S["robustness"] = robustness(RECS, HYPS)
    rb = S["robustness"]
    log(f"weight perturbation: top-1 unchanged in {rb['weight_perturbation']['top1_unchanged']:.1%} of "
        f"{rb['weight_perturbation']['n']} runs (worst rank {rb['weight_perturbation']['worst_rank']})")
    log(f"rank of the top hypothesis without ClinVar: {rb['rank_without_clinvar']}; without the missense predictor: "
        f"{rb['rank_without_missense_predictor']}; without the diagnosis: {rb['rank_without_diagnosis']}; "
        f"without both ClinVar and the diagnosis: {rb['rank_without_clinvar_and_diagnosis']}")
    md("**Leave-one-HPO-term-out (blind mode)**\n\n" + df_md(pd.DataFrame(rb["leave_one_term_out_blind"])))


# %%
# ---- the organisers' public Track 1 scorer (evaluation.py of the challenge Space), embedded verbatim in logic so the
# ---- file we write is parsed and scored exactly as it will be on submission
RANK_POINT_TIERS = [(1, 100), (3, 50), (5, 25), (10, 10)]


def official_load(csv_path):
    import csv as _csv
    rows = []
    with open(csv_path, newline="") as f:
        for row in _csv.DictReader(f):
            def pv(c, p, r, a):
                return None if (not c or not p) else (c.strip(), int(p), r.strip().upper(), a.strip().upper())
            v1 = pv(row["chrom_1"], row["pos_1"], row["ref_1"], row["alt_1"])
            v2 = pv(row.get("chrom_2", ""), row.get("pos_2", ""), row.get("ref_2", ""), row.get("alt_2", ""))
            if v1 is None:
                raise ValueError("row missing primary variant")
            epcr = float(row["epcr"])
            if not (0 < epcr <= 1):
                raise ValueError(f"EPCR {epcr} out of range (0,1]")
            ft = (row.get("finding_type") or "primary").strip().lower()
            if ft not in ("primary", "secondary"):
                raise ValueError(f"bad finding_type {ft}")
            rows.append((row["proband_id"].strip(), frozenset([v1, v2]) if v2 else frozenset([v1]), epcr, ft))
    if len(rows) > 10:
        raise ValueError(f"{len(rows)} rows; max is 10")
    order = sorted(range(len(rows)), key=lambda i: (-rows[i][2], i))
    return [(rows[i][0], rows[i][1], rows[i][2], rank + 1) for rank, i in enumerate(order)]


def official_score(sub_rows, true_variants):
    full = next((r for _, vs, _, r in sub_rows if vs == true_variants), None)
    partial = None
    if len(true_variants) == 2 and full is None:
        partial = next((r for _, vs, _, r in sub_rows if vs & true_variants), None)
    pts = lambda r: next((p for m, p in RANK_POINT_TIERS if r <= m), 0)
    rank_points = float(pts(full)) if full else (0.5 * pts(partial) if partial else 0.0)
    best_f = 0.0
    for t in sorted({e for _, _, e, _ in sub_rows}, reverse=True):
        pred = set()
        for _, vs, e, _ in sub_rows:
            if e >= t:
                pred |= vs
        tp, fp, fn = len(pred & true_variants), len(pred - true_variants), len(true_variants - pred)
        pr = tp / (tp + fp) if tp + fp else 0.0
        rc = tp / (tp + fn) if tp + fn else 0.0
        best_f = max(best_f, 2 * pr * rc / (pr + rc) if pr + rc else 0.0)
    return dict(rank_points=rank_points, f_max=round(best_f, 3), full_match_rank=full, partial_match_rank=partial)


SUBMISSION_CSV = OUT / f"{CFG['TEAM']}_{CFG['MODEL_NAME']}.csv"


def select_rows(hyps, max_rows=10, per_gene=5, min_epcr=0.01, min_epcr_other_gene=0.25):
    """Rows for the CSV: the top hypothesis, alternatives within its gene (other pairings, single alleles: partial
    credit if only one allele is right), and hypotheses in other genes only when they are serious competitors."""
    rows, seen, n_gene = [], set(), collections.Counter()
    for h in hyps:
        if rows and h["epcr"] < min_epcr:
            break
        if rows and h["gene"] != rows[0]["gene"] and h["epcr"] < min_epcr_other_gene:
            continue
        if h["key"] in seen or n_gene[h["gene"]] >= per_gene:
            continue
        if rows and len(h["key"]) == 1 and h["key"] <= frozenset().union(*(r["key"] for r in rows)):
            continue        # a single allele already listed in a higher pair adds nothing to rank points or F-max
        if rows and any(v.get("points", 0) < 0 for v in h["variants"]):
            continue        # a variant with net benign ACMG evidence is not put forward as a causal candidate
        seen.add(h["key"])
        n_gene[h["gene"]] += 1
        rows.append(h)
        if len(rows) >= max_rows:
            break
    return rows


def write_submission(hyps, secondary=None):
    """Top hypotheses -> the organisers' CSV format (chr-prefixed GRCh38, integer positions, strictly decreasing EPCR)."""
    import csv as _csv
    secondary = (secondary or [])[:2] if CFG["INCLUDE_SECONDARY_ROWS"] else []
    primary = select_rows(hyps, max_rows=10 - len(secondary))
    out_rows, prev = [], None
    for h, ftype in [(h, "primary") for h in primary] + [(h, "secondary") for h in secondary]:
        e = round(min(h["epcr"], 0.99), 4)
        if ftype == "secondary":
            e = min(e, 0.009)
        if prev is not None and e >= prev:
            e = round(prev - 0.0005, 4)
        e = max(e, 0.0001)
        prev = e
        vs = sorted(h["variants"], key=lambda v: (PRIMARY.index(v["chrom"]) if v["chrom"] in PRIMARY else 99, int(v["pos"])))
        v1, v2 = vs[0], (vs[1] if len(vs) > 1 else None)
        note = h.get("note") or f"{h['gene']} {h['model'].replace('_', '-')}" + (f", phase {h['phase']}" if v2 else "")
        out_rows.append([CFG["PROBAND_ID"], v1["chrom"], int(v1["pos"]), v1["ref"], v1["alt"],
                         v2["chrom"] if v2 else "", int(v2["pos"]) if v2 else "", v2["ref"] if v2 else "", v2["alt"] if v2 else "",
                         f"{e:.4f}", ftype, note.replace(",", ";")])
    with open(SUBMISSION_CSV, "w", newline="") as f:
        wr = _csv.writer(f)
        wr.writerow(["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"])
        wr.writerows(out_rows)
    parsed = official_load(SUBMISSION_CSV)             # raises if the organisers' parser would reject the file
    assert parsed[0][1] == primary[0]["key"], "top row of the CSV is not the top hypothesis"
    assert len({e for _, _, e, _ in parsed}) == len(parsed), "EPCR values must be strictly decreasing"
    return primary, parsed


def verify_in_source(variants):
    """Every submitted allele must be present in the annotated copy of the source VCF (and must not be synthetic)."""
    res = []
    for v in variants:
        out = sh(f"bcftools view -H -r {v['chrom']}:{int(v['pos'])}-{int(v['pos'])} '{ANN_BCF}' | cut -f1-5", name="verify", capture=True)
        hit = [l.split("\t") for l in out.splitlines() if l.strip()]
        ok = [h for h in hit if int(h[1]) == int(v["pos"]) and h[3] == v["ref"] and h[4] == v["alt"]]
        res.append(dict(variant=f"{v['chrom']}:{int(v['pos'])} {v['ref']}>{v['alt']}", found=bool(ok),
                        synthetic=bool(ok and ok[0][2].startswith("SPIKE_"))))
    return res


with stage("write + validate submission (first pass)", required=True):
    SUB_ROWS, parsed = write_submission(HYPS)
    chk = verify_in_source([v for h in SUB_ROWS for v in h["variants"]])
    S["source_check"] = chk
    if not DEMO:
        assert all(c["found"] and not c["synthetic"] for c in chk), f"submitted allele not found in the source VCF: {chk}"
    if DEMO:
        truth = frozenset(SPIKE_CASES[0]["variants"])
        S["demo_score"] = official_score(parsed, truth)
        log(f"DEMO self-score against the spiked truth (organisers' scorer): {S['demo_score']}")
    log(f"submission written: {SUBMISSION_CSV.name} ({len(parsed)} rows)")
    md("**Submission rows (first pass)**\n\n" + df_md(pd.read_csv(SUBMISSION_CSV, dtype=str, keep_default_na=False)))

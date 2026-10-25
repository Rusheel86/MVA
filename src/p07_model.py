# %% [markdown]
# ## 7 · Evidence model and genotype hypotheses
#
# **Variant level.** Each rare variant collects ACMG/AMP evidence in Bayesian points (Tavtigian et al. 2018, 2020:
# supporting = 1, moderate = 2, strong = 4, very strong = 8; prior 0.10; odds of pathogenicity 350 for very strong), so
# `P = O·0.1 / ((O−1)·0.1 + 1)` with `O = 350^(points/8)`. All thresholds are taken from the literature and were
# fixed before the pipeline was run; none is specific to a gene or a case.
#
# | Evidence | Automated rule | Points |
# |---|---|---|
# | PVS1 | stop-gained / frameshift / canonical ±1,2 splice on the MANE transcript, in a gene with an established loss-of-function mechanism; nonsense-mediated decay predicted (Abou Tayoun 2018) | +8; NMD-escaping: +4 if > 10 % of the protein is lost, else +2; mechanism not established: +2 |
# | PM2 | absent from gnomAD, or filtering AF ≤ 1e-4 with no homozygotes (ClinGen SVI PM2 v1.0: supporting) | +1 |
# | BS1 / BA1 / BS2 | AF > 0.1 % (possible hypomorph: −1), ≥ 1 % (−4), ≥ 5 % (−8); ≥ 2 homozygotes in gnomAD (−4). Waived for ClinVar P/LP alleles below BA1 | −1 … −8 |
# | PP3 / BP4 missense | AlphaMissense in ClinGen-SVI calibrated bands (Pejaver 2022; Bergquist 2025): ≥ 0.990 / ≥ 0.972 / ≥ 0.906 / ≥ 0.792; ≤ 0.169 / ≤ 0.099 / ≤ 0.070. REVEL bands are used only when AlphaMissense is missing | +4 / +3 / +2 / +1; −1 / −2 / −3 |
# | PP3 / BP4 splice | SpliceAI max delta ≥ 0.2 (Walker 2023; graded ≥ 0.5, ≥ 0.8); ≤ 0.1 | +1 / +2 / +4; −1 (synonymous: −2, BP7) |
# | PM4 | in-frame indel, stop-loss | +2 |
# | ClinVar | used as a documented prior, never as a filter: ≥ 2-star Pathogenic sets P ≥ 0.99, P/LP ≥ 0.95, LP ≥ 0.90; 1-star P/LP +2; benign assertions mirrored | floor / cap |
# | PM3 | applied once, to the weaker allele of a pair (ClinGen SVI PM3 v1.0, single proband): partner Pathogenic — *trans* +2, phase unknown +1; partner Likely pathogenic — *trans* +1, unknown 0; homozygous +1 | genotype level |
#
# **Genotype level.** Hypotheses are enumerated per gene — homozygous / hemizygous, compound-heterozygous pairs, single
# heterozygous — and scored as
#
# `EPCR = P(genotype pathogenic) × phenotype fit × inheritance-mode fit × call quality`
#
# A compound-het pair multiplies the two variant posteriors and a phase factor (read-backed *trans* 1.0; unknown 0.9 —
# if the gene is causal and both alleles are pathogenic they must be in trans, bar a cis pair plus a hidden third
# allele; read-backed *cis* 0.01). The phenotype enters once, through the gene-level phenotype fit (PP4 is therefore
# not added to variant points). The heuristic weights live in one dictionary and are stress-tested in section 9.

# %%
PRIOR_P, ODDS_VS = 0.10, 350.0
LOF_CSQ = {"stop_gained", "frameshift", "splice_donor", "splice_acceptor"}
NONCODING_CSQ = {"synonymous", "intron", "splice_region", "utr5", "utr3", "upstream", "downstream", "other"}

W = dict(
    phase_trans=1.0, phase_unknown=0.9, phase_cis=0.01,
    mode_biallelic_ok=1.0, mode_biallelic_unknown=0.6, mode_biallelic_dominant_gene=0.4,
    mode_het_dominant=0.5, mode_het_unknown=0.1, mode_het_single_hit=0.1, mode_het_recessive_gene=0.02,
    q_nonpass=0.3, q_allele_balance=0.5, q_low_gq=0.5,
    pheno_floor=0.02, pheno_power=3.0, dx_floor=0.95,
)
AM_BANDS = [(0.990, 4), (0.972, 3), (0.906, 2), (0.792, 1)]
AM_BENIGN = [(0.070, -3), (0.099, -2), (0.169, -1)]
REVEL_BANDS = [(0.932, 4), (0.879, 3), (0.773, 2), (0.644, 1)]
REVEL_BENIGN = [(0.016, -4), (0.052, -3), (0.183, -2), (0.290, -1)]


def points_to_posterior(points):
    odds = ODDS_VS ** (points / 8.0)
    return odds * PRIOR_P / ((odds - 1.0) * PRIOR_P + 1.0)


def clinvar_stars(revstat):
    r = (revstat or "").lower()
    if "practice_guideline" in r:
        return 4
    if "expert_panel" in r:
        return 3
    if r.startswith("no_") or "no_assertion" in r or "no_classification" in r:
        return 0
    if "multiple_submitters" in r and "no_conflict" in r:
        return 2
    return 1 if "criteria_provided" in r else 0


def clinvar_class(sig):
    s = (sig or "").lower()
    if not s or s == ".":
        return None
    if s.startswith("conflicting"):
        return "conflicting"
    if "pathogenic" in s and "benign" not in s:
        return "P/LP"
    if "benign" in s and "pathogenic" not in s:
        return "B/LB"
    if "uncertain" in s:
        return "VUS"
    return "other"


def _num(x):
    """None for missing values (NaN / None / '.' / '')."""
    if x is None:
        return None
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(x) else x


def _band(value, pathogenic, benign):
    for thr, pts in pathogenic:
        if value >= thr:
            return pts
    for thr, pts in benign:
        if value <= thr:
            return pts
    return 0


def variant_evidence(v, lof_ok=True):
    """(points, floor, cap, [(code, points, explanation)]) for one annotated variant (dict)."""
    ev, floor, cap = [], 0.0, 1.0
    af, nhom = _num(v.get("af")), _num(v.get("nhomalt"))
    cls, stars = clinvar_class(v.get("clnsig")), clinvar_stars(v.get("clnrevstat"))
    known_path = cls == "P/LP" and stars >= 1
    # ---- population data --------------------------------------------------------------------------------------------
    if af is None:
        ev.append(("PM2_supporting", 1, "absent from gnomAD"))
    elif af >= 0.05:
        ev.append(("BA1", -8, f"gnomAD AF {af:.2g} >= 5%"))
    elif af >= 0.01:
        if not known_path:
            ev.append(("BS1", -4, f"gnomAD AF {af:.2g} >= 1%"))
    elif af > 1e-3:
        if not known_path:
            ev.append(("BS1_supporting", -1, f"gnomAD AF {af:.2g} > 0.1% (too common for a fully penetrant allele; possible hypomorph)"))
    elif af <= 1e-4 and not nhom:
        ev.append(("PM2_supporting", 1, f"gnomAD AF or FAF95 (the higher) {af:.2g} <= 1e-4, no homozygotes"))
    if nhom is not None and nhom >= 2 and not known_path:
        ev.append(("BS2", -4, f"{int(nhom)} homozygotes in gnomAD"))
    ov = _num(v.get("common_overlap"))
    if ov is not None and ov >= 0.01 and not known_path:
        ev.append(("BS1_overlap", -4, f"overlaps a common gnomAD indel (AF {ov:.2g}): representation artefact likely"))
    # ---- ClinVar: documented prior, never a requirement ------------------------------------------------------------
    sig = (v.get("clnsig") or "").lower()
    if cls == "P/LP":
        if stars >= 2:
            floor = 0.90 if sig.startswith("likely_pathogenic") else (0.95 if "likely_pathogenic" in sig else 0.99)
            ev.append((f"ClinVar_{stars}star", 0, f"ClinVar {v.get('clnsig')} ({stars}-star): posterior >= {floor}"))
        else:
            ev.append((f"ClinVar_{stars}star", 2 if stars == 1 else 1, f"ClinVar {v.get('clnsig')} ({stars}-star)"))
    elif cls == "B/LB":
        if stars >= 2:
            cap = 0.01
            ev.append((f"ClinVar_benign_{stars}star", 0, f"ClinVar {v.get('clnsig')} ({stars}-star): posterior <= 0.01"))
        else:
            ev.append((f"ClinVar_benign_{stars}star", -2 if stars == 1 else -1, f"ClinVar {v.get('clnsig')}"))
    # ---- consequence -------------------------------------------------------------------------------------------------
    c = v.get("csq")
    sp = _num(v.get("spliceai"))
    if c in LOF_CSQ:
        frac = _num(v.get("frac_lost"))
        if not lof_ok:
            ev.append(("PVS1_moderate", 2, f"{c}; loss of function is not an established mechanism for this gene"))
        elif v.get("nmd_escape") is True:
            if frac is not None and frac > 0.10:
                ev.append(("PVS1_strong", 4, f"{c}; escapes NMD, {frac:.0%} of the protein lost"))
            else:
                ev.append(("PVS1_moderate", 2, f"{c}; escapes NMD, < 10% of the protein lost"))
        else:
            ev.append(("PVS1", 8, f"{c}; nonsense-mediated decay predicted"))
    elif c == "start_lost":
        ev.append(("PVS1_moderate", 2, "start codon lost"))
    elif c in ("stop_lost", "inframe_indel"):
        ev.append(("PM4", 2, c.replace("_", " ")))
    elif c == "missense":
        am, revel = _num(v.get("am")), _num(v.get("revel"))
        tool, val, pts = (("AlphaMissense", am, _band(am, AM_BANDS, AM_BENIGN)) if am is not None else
                          (("REVEL", revel, _band(revel, REVEL_BANDS, REVEL_BENIGN)) if revel is not None else (None, None, 0)))
        if pts:
            name = {4: "strong", 3: "moderate+", 2: "moderate", 1: "supporting"}[abs(pts)]
            ev.append((f"{'PP3' if pts > 0 else 'BP4'}_{name}", pts, f"{tool} {val:.3f}"))
    if c not in LOF_CSQ:
        if sp is not None and sp >= 0.2:
            code, pts = ("PP3_splice_strong", 4) if sp >= 0.8 else (("PP3_splice_moderate", 2) if sp >= 0.5 else ("PP3_splice_supporting", 1))
            ev.append((code, pts, f"SpliceAI max delta {sp:.2f}"))
        elif c in NONCODING_CSQ and not (cls == "P/LP"):
            if sp is not None and sp <= 0.1:
                ev.append(("BP4_BP7" if c == "synonymous" else "BP4", -2 if c == "synonymous" else -1, f"{c}; SpliceAI {sp:.2f}"))
            elif sp is None:
                ev.append(("BP4_default", -1, f"{c}; no predicted functional effect (splice model not run)"))
    return sum(p for _, p, _ in ev), floor, cap, ev


def posterior(v, extra=0.0):
    """Variant posterior with optional genotype-level points (PM3), bounded by the ClinVar floor / cap."""
    return float(min(v["cap"], max(v["floor"], points_to_posterior(v["points"] + extra))))


def quality_factor(v, w=W):
    q = 1.0
    if str(v.get("filt", "PASS")) not in ("PASS", "."):
        q *= w["q_nonpass"]
    a, b = _num(v.get("ad_ref")), _num(v.get("ad_alt"))
    if a is not None and b is not None and (a + b) > 0 and v.get("zyg") == "het":
        if not 0.25 <= b / (a + b) <= 0.75:
            q *= w["q_allele_balance"]
    gq = _num(v.get("gq"))
    if gq is not None and gq < 20:
        q *= w["q_low_gq"]
    return q


def phase_of(a, b):
    """Read-backed phase from GATK/Sentieon PID/PGT tags (same PID: same PGT = cis, opposite = trans)."""
    pa, pb = str(a.get("pid") or "."), str(b.get("pid") or ".")
    ga, gb = str(a.get("pgt") or "."), str(b.get("pgt") or ".")
    if pa not in (".", "", "nan") and pa == pb and ga in ("0|1", "1|0") and gb in ("0|1", "1|0"):
        return "cis" if ga == gb else "trans"
    return "unknown"


def pheno_factor(info, top, w=W):
    """Gene-level phenotype fit in (0, 1]: cubic in the HPO similarity relative to the best-matching gene; genes
    without a Mendelian phenotype annotation get a small floor; genes on the indication's panel are lifted."""
    if not info:
        return w["pheno_floor"]
    f = max(w["pheno_floor"], (info["pheno"] / top) ** w["pheno_power"]) if top > 0 else w["pheno_floor"]
    if info.get("dx_match"):
        f = max(f, w["dx_floor"])
    return min(f, 1.0)


def vkey(v):
    return (v["chrom"], int(v["pos"]), v["ref"], v["alt"])


def score_variants(cand_df, gene_table, lof_genes):
    recs = cand_df.to_dict("records")
    for v in recs:
        inh = set((gene_table.get(v["gene"]) or {}).get("inheritance_all") or [])
        lof_ok = v["gene"] in lof_genes or bool(inh & {"AR", "XLR", "XL"})
        v["points"], v["floor"], v["cap"], v["evidence"] = variant_evidence(v, lof_ok=lof_ok)
        v["p_path"] = posterior(v)
    return recs


def build_hypotheses(recs, gene_table, pheno_top, w=W, max_pair_variants=14):
    by_gene = collections.defaultdict(list)
    for v in recs:
        by_gene[v["gene"]].append(v)
    hyps = []
    for g, vs in by_gene.items():
        info = gene_table.get(g) or {}
        inh = set(info.get("inheritance") or [])
        pf = pheno_factor(info, pheno_top, w)
        recessive, dominant = bool(inh & {"AR", "XLR", "XL"}), bool(inh & {"AD", "XLD", "XL"})
        bi_mode = w["mode_biallelic_ok"] if recessive else (w["mode_biallelic_dominant_gene"] if dominant else w["mode_biallelic_unknown"])
        for v in vs:
            q = quality_factor(v, w)
            af = _num(v.get("af"))
            if v["zyg"] in ("hom", "hemi"):
                extra = 1.0 if (af is None or af <= 1e-4) else (-3.0 if af > 1e-3 else 0.0)     # PM3 (homozygous) / BS1 in a homozygous model
                hyps.append(dict(gene=g, model="homozygous" if v["zyg"] == "hom" else "hemizygous", variants=[v], phase="n/a",
                                 p_geno=posterior(v, extra), pheno=pf, mode=bi_mode, quality=q))
            else:
                if dominant:
                    mode = w["mode_het_dominant"] * (1.0 if (af is None or af <= 1e-4) else 0.1)        # dominant rarity gate
                elif recessive:
                    mode = w["mode_het_single_hit"] if v["p_path"] >= 0.9 else w["mode_het_recessive_gene"]
                else:
                    mode = w["mode_het_unknown"]
                hyps.append(dict(gene=g, model="heterozygous", variants=[v], phase="n/a", p_geno=v["p_path"], pheno=pf,
                                 mode=mode, quality=q))
        hets = sorted((v for v in vs if v["zyg"] == "het"), key=lambda v: -v["p_path"])[:max_pair_variants]
        for a, b in itertools.combinations(hets, 2):       # a is the stronger allele
            if vkey(a) == vkey(b):
                continue
            ph = phase_of(a, b)
            b_af = _num(b.get("af"))
            pm3 = 0.0
            if ph != "cis" and (b_af is None or b_af <= 1e-4):         # PM3 once, on the weaker allele, if it meets PM2
                if a["p_path"] >= 0.99:                                 # partner Pathogenic: trans moderate, unknown supporting
                    pm3 = 2.0 if ph == "trans" else 1.0
                elif a["p_path"] >= 0.90:                               # partner Likely pathogenic: trans supporting, unknown none
                    pm3 = 1.0 if ph == "trans" else 0.0
            pa, pb = a["p_path"], posterior(b, pm3)
            hyps.append(dict(gene=g, model="compound_het", variants=[a, b], phase=ph, p_geno=pa * pb * w[f"phase_{ph}"],
                             pheno=pf, mode=bi_mode, quality=quality_factor(a, w) * quality_factor(b, w),
                             p_each=(round(pa, 4), round(pb, 4)), pm3=pm3))
    for h in hyps:
        h["epcr"] = float(min(0.99, max(1e-4, h["p_geno"] * h["pheno"] * h["mode"] * h["quality"])))
        h["key"] = frozenset(vkey(v) for v in h["variants"])
    hyps.sort(key=lambda h: (-h["epcr"], -len(h["variants"])))
    return hyps


def describe_variant(v):
    hg = v.get("hgvsp") or ""
    return f"{v['chrom']}:{int(v['pos'])} {v['ref']}>{v['alt']} {v['csq']}" + (f" {hg}" if hg else "")


def hyp_table(hyps, n=15):
    rows = []
    for i, h in enumerate(hyps[:n], 1):
        rows.append(dict(rank=i, EPCR=round(h["epcr"], 4), gene=h["gene"], model=h["model"], phase=h["phase"],
                         variants=" + ".join(describe_variant(v) for v in h["variants"]),
                         evidence=" + ".join("[" + ", ".join(c for c, _, _ in v["evidence"]) + "]" for v in h["variants"]),
                         P_genotype=round(h["p_geno"], 3), phenotype_fit=round(h["pheno"], 2), mode_fit=h["mode"],
                         quality=round(h["quality"], 2)))
    return pd.DataFrame(rows)


def reportable_gene(g):
    """Gene names outside the stated indication's panel and the leading hypothesis are withheld from displayed tables
    and reports (organisers' guidance on incidental findings is pending); set CFG['MASK_INCIDENTAL']=False to show them."""
    if DEMO or not CFG.get("MASK_INCIDENTAL", True):
        return True
    return g == S.get("top_gene") or g in set((S.get("pheno") or {}).get("dx_panel", []))


def public_table(df):
    if "gene" not in df.columns or not len(df):
        return df
    df = df.copy()
    hide = ~df["gene"].map(reportable_gene)
    for c in ("variants", "evidence"):
        if c in df.columns:
            df.loc[hide, c] = "(withheld)"
    df.loc[hide, "gene"] = "other gene (withheld)"
    return df

# %% [markdown]
# ## 13 · Does the method generalise? A spike-in benchmark on a public genome
#
# One solved case says little about a method. The same pipeline is therefore run on synthetic recessive cases built
# on the public GIAB HG002 genome: for each case, one ClinVar pathogenic truncating allele and one ClinVar pathogenic
# missense allele of a single-gene autosomal-recessive disease are spiked in as an unphased heterozygous pair, and
# the "patient" is given six HPO terms sampled from that disease — some replaced by a less specific parent term, plus
# one unrelated noise term. The ranking is **blind to the diagnosis** (HPO terms only). Two settings are reported:
#
# * **ClinVar-aware** — the realistic setting for known alleles;
# * **ClinVar-blind** — the spiked alleles' ClinVar records are hidden, so they must be found from consequence,
#   rarity and the missense predictor alone (a proxy for novel alleles).
#
# The benchmark uses only public data and can be reproduced by anyone without access to the patient genome.

# %%
def run_benchmark(cand_all, pool_all, cases):
    bg = cand_all[~cand_all.rsid.str.startswith("SPIKE_")].reset_index(drop=True)
    pool_bg = pool_all[~pool_all.rsid.str.startswith("SPIKE_")].reset_index(drop=True)
    bg_recs = score_variants(bg, GENES_BLIND, LOF_GENES)
    rows = []
    for cs in cases:
        ids = [f"SPIKE_{cs['case']}_0", f"SPIKE_{cs['case']}_1"]
        sp = cand_all[cand_all.rsid.isin(ids) & (cand_all.gene == cs["gene"])]
        truth = frozenset(cs["variants"])
        found = frozenset(zip(sp.chrom, sp.pos.astype(int), sp.ref, sp.alt))
        if found != truth:      # an allele did not survive normalisation / annotation as a candidate in the intended gene
            rows.append(dict(case=cs["case"], gene=cs["gene"], disease=cs["disease_name"][:60], rank_aware=None, rank_blind=None,
                             epcr_aware=None, epcr_blind=None, top_gene_blind=None,
                             note="spiked allele not annotated as a candidate in the intended gene"))
            continue
        gt, _ = PE.score_genes(cs["hpo"], None, diagnosis=None)
        ptop = max(v["pheno"] for v in gt.values()) or 1.0
        res = {}
        for mode in ("aware", "blind"):
            sp2 = sp.copy()
            if mode == "blind":
                sp2["clnsig"], sp2["clnrevstat"] = "", ""
            recs = bg_recs + score_variants(sp2, gt, LOF_GENES)
            hyps = build_hypotheses(recs, gt, ptop, W)
            r = next((i + 1 for i, h in enumerate(hyps) if h["key"] == truth), None)
            res[mode] = (r, hyps[r - 1]["epcr"] if r else None, hyps[0]["gene"])
        rows.append(dict(case=cs["case"], gene=cs["gene"], disease=cs["disease_name"][:60], rank_aware=res["aware"][0],
                         rank_blind=res["blind"][0], epcr_aware=res["aware"][1], epcr_blind=res["blind"][1],
                         top_gene_blind=res["blind"][2], note=""))
    return pd.DataFrame(rows)


with stage("spike-in benchmark (public genome)"):
    if not CFG["RUN_BENCHMARK"]:
        raise RuntimeError("disabled in CFG")
    if DEMO:
        bench_cand, bench_pool, bench_cases = CAND_ALL, POOL_ALL, SPIKE_CASES
    else:
        hg = fetch(HG002_VCF, WORK / "data" / "HG002.vcf.gz", min_bytes=100_000_000)
        fetch(HG002_VCF + ".tbi", WORK / "data" / "HG002.vcf.gz.tbi", min_bytes=100_000)
        bench_cases = make_spike_cases(60)
        bench_cand, bench_pool, _, _ = annotate_vcf(str(hg), "bench", spikes=bench_cases, sex="unknown")
    BENCH = run_benchmark(bench_cand, bench_pool, bench_cases)
    ok = BENCH[BENCH.rank_aware.notna() | BENCH.rank_blind.notna()]
    n = len(ok)

    def rate(col, k):
        return float((ok[col].fillna(10 ** 6) <= k).mean()) if n else float("nan")

    S["benchmark"] = dict(n_cases=int(n), n_not_evaluable=int(len(BENCH) - n),
                          aware=dict(top1=round(rate("rank_aware", 1), 3), top3=round(rate("rank_aware", 3), 3), top10=round(rate("rank_aware", 10), 3)),
                          blind=dict(top1=round(rate("rank_blind", 1), 3), top3=round(rate("rank_blind", 3), 3), top10=round(rate("rank_blind", 10), 3)),
                          median_epcr_when_top1=round(float(ok.loc[ok.rank_aware == 1, "epcr_aware"].median()), 3) if (ok.rank_aware == 1).any() else None)
    log(f"benchmark on {n} synthetic recessive cases (HPO only, no diagnosis): {S['benchmark']}")
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    ks, x = [1, 3, 10], np.arange(3)
    for j, (mode, label, col) in enumerate((("aware", "ClinVar-aware", VIZ["s1"]), ("blind", "ClinVar-blind (alleles treated as novel)", VIZ["s2"]))):
        vals = [S["benchmark"][mode][f"top{k}"] * 100 for k in ks]
        bars = ax.bar(x + (j - 0.5) * 0.36, vals, width=0.32, color=col, label=label)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 1.5, f"{v:.0f}%", ha="center", va="bottom", fontsize=8, color=VIZ["ink2"])
    ax.set_xticks(x)
    ax.set_xticklabels([f"top {k}" for k in ks])
    ax.set_ylim(0, 112)
    ax.set_ylabel("cases with the true pair recovered (%)")
    ax.set_title(f"Spike-in benchmark: {n} synthetic recessive cases on GIAB HG002")
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
    save_fig(fig, "benchmark")
    miss = ok[(ok.rank_blind.fillna(10 ** 6) > 1)][["case", "gene", "disease", "rank_aware", "rank_blind", "top_gene_blind"]]
    if len(miss):
        md("**Cases where the true pair was not ranked first in the ClinVar-blind setting**\n\n" + df_md(miss, max_rows=15))

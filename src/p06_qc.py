# %% [markdown]
# ## 6 · Genome QC, sex, autozygosity and a chromosome-level aneuploidy screen
#
# The disorder in question is defined cytogenetically (mosaic aneuploidies that differ from cell to cell), so before
# ranking variants we ask what the bulk genome itself can and cannot show:
#
# * **QC** – Ti/Tv, het/hom ratio, depth: is this a clean germline call set?
# * **Sex and autozygosity** – X heterozygosity / depth, and long runs of homozygosity (parental relatedness changes the
#   prior between homozygous and compound-heterozygous models).
# * **Aneuploidy screen** – per-chromosome depth ratio and B-allele-frequency dispersion at heterozygous SNVs. A *clonal*
#   gain or loss shifts both; *variegated* aneuploidy (a different chromosome in every cell) averages out in bulk DNA and
#   is expected to be invisible. The screen therefore reports a detection limit rather than a diagnosis.

# %%
def load_site_table(vcf_path):
    """Compact per-site numeric table for QC (biallelic records only). Stays on scratch."""
    tsv = WORK / "tmp" / "qc_sites.tsv"
    awk = r'''BEGIN{OFS="\t"} {
        n=split($7,a,","); gt=$6; gsub(/\|/,"/",gt);
        z=(gt=="0/1"||gt=="1/0")?1:((gt=="1/1")?2:0);
        snv=(length($3)==1&&length($4)==1)?1:0; ti=0;
        if(snv){p=$3 $4; if(p=="AG"||p=="GA"||p=="CT"||p=="TC") ti=1}
        print $1,$2,snv,ti,z,a[1]+0,a[2]+0,$8+0,$9+0,(($5=="PASS"||$5==".")?1:0),((snv&&($3=="G"||$3=="C"))?1:0) }'''
    sh(f"bcftools query -i 'N_ALT=1' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER[\\t%GT\\t%AD\\t%DP\\t%GQ]\\n' '{vcf_path}' "
       f"| awk -F'\\t' '{awk}' > '{tsv}'", name="qc_query")
    t = pd.read_csv(tsv, sep="\t", header=None,
                    names=["chrom", "pos", "snv", "ti", "zyg", "ad_ref", "ad_alt", "dp", "gq", "is_pass", "gc"],
                    dtype={"chrom": "category", "pos": np.int32, "snv": np.int8, "ti": np.int8, "zyg": np.int8,
                           "ad_ref": np.int32, "ad_alt": np.int32, "dp": np.int32, "gq": np.int32, "is_pass": np.int8,
                           "gc": np.int8})
    tsv.unlink(missing_ok=True)
    t["chrom"] = t["chrom"].astype(str).str.replace("^chr", "", regex=True)
    return t


AUTOSOMES = [str(i) for i in range(1, 23)]
PAR_X = [(10_001, 2_781_479), (155_701_383, 156_030_895)]        # GRCh38 pseudo-autosomal regions


def qc_summary(t):
    q = {}
    p = t[t.is_pass == 1]
    q["biallelic_records"] = int(len(t))
    q["pass_fraction"] = round(float(len(p) / max(len(t), 1)), 4)
    snv = p[p.snv == 1]
    q["pass_snvs"], q["pass_indels"] = int(len(snv)), int(len(p) - len(snv))
    q["titv"] = round(float(snv.ti.sum() / max((snv.ti == 0).sum(), 1)), 3)
    auto = snv[snv.chrom.isin(AUTOSOMES)]
    q["het_hom_ratio"] = round(float((auto.zyg == 1).sum() / max((auto.zyg == 2).sum(), 1)), 3)
    hq = auto[(auto.zyg == 1) & (auto.gq >= 60)]
    q["median_depth_het_sites"] = float((hq.ad_ref + hq.ad_alt).median()) if len(hq) else float("nan")
    # sex from chrX: heterozygosity outside the PARs and depth relative to autosomes
    x = snv[(snv.chrom == "X") & (snv.gq >= 60)]
    x = x[~(((x.pos >= PAR_X[0][0]) & (x.pos <= PAR_X[0][1])) | ((x.pos >= PAR_X[1][0]) & (x.pos <= PAR_X[1][1])))]
    if len(x) > 500:
        q["x_het_fraction"] = round(float((x.zyg == 1).mean()), 4)
        xd = float((x[x.zyg == 2].ad_alt).median())
        ad = float((auto[auto.zyg == 2].ad_alt).median())
        q["x_to_autosome_depth"] = round(xd / ad, 3) if ad else float("nan")
        q["inferred_sex"] = "male (XY)" if (q["x_het_fraction"] < 0.15 and q["x_to_autosome_depth"] < 0.75) else "female (XX)"
    else:
        q["inferred_sex"] = "not determined (no chrX records)"
    return q


def runs_of_homozygosity(t, min_len=1_000_000, min_hom=100):
    """Het-gap method: a stretch >= min_len between consecutive confident het SNVs that still contains >= min_hom
    hom-alt SNVs is called a run of homozygosity (centromeres/gaps have no variants at all and are not counted)."""
    rows = []
    g = t[(t.is_pass == 1) & (t.snv == 1) & (t.gq >= 30) & t.chrom.isin(AUTOSOMES)]
    for ch, sub in g.groupby("chrom", observed=True):
        het = np.sort(sub.pos[sub.zyg == 1].to_numpy())
        hom = np.sort(sub.pos[sub.zyg == 2].to_numpy())
        if len(het) < 2:
            continue
        gaps = np.diff(het)
        for i in np.where(gaps >= min_len)[0]:
            a, b = het[i], het[i + 1]
            n_hom = int(np.searchsorted(hom, b) - np.searchsorted(hom, a))
            if n_hom >= min_hom and n_hom / ((b - a) / 1e6) >= 50:
                rows.append(dict(chrom=ch, start=int(a), end=int(b), mb=round((b - a) / 1e6, 2), n_hom=n_hom))
    roh = pd.DataFrame(rows, columns=["chrom", "start", "end", "mb", "n_hom"])
    return roh, float(roh.mb.sum() / 2875.0) if len(roh) else 0.0        # ~2,875 Mb of autosome


def _trimmed_mean(x):
    lo, hi = x.quantile(0.05), x.quantile(0.95)
    return float(x[(x >= lo) & (x <= hi)].mean())


def aneuploidy_screen(t, bin_mb=5, min_sites=300):
    """Depth and allele balance at confident heterozygous SNVs, per 5-Mb bin and per autosome.
    Depth is corrected for base composition with a within-sample proxy (GC fraction of the reference bases at
    SNV sites in the bin); a chromosome is flagged only if depth AND allele-balance dispersion deviate together."""
    p = t[(t.is_pass == 1) & (t.snv == 1) & t.chrom.isin(AUTOSOMES)].copy()
    p["bin"] = (p.pos // (bin_mb * 1_000_000)).astype(int)
    gc = p.groupby(["chrom", "bin"], observed=True)["gc"].mean().rename("gc_proxy")
    h = p[(p.zyg == 1) & (p.gq >= 60) & ((p.ad_ref + p.ad_alt) >= 15)].copy()
    h["depth"] = h.ad_ref + h.ad_alt
    h["ab"] = h.ad_alt / h.depth
    bins = h.groupby(["chrom", "bin"], observed=True).agg(
        n=("ab", "size"), depth=("depth", _trimmed_mean), baf_dev=("ab", lambda s: float(np.mean(np.abs(s - 0.5)))),
        skew_frac=("ab", lambda s: float(((s < 0.35) | (s > 0.65)).mean()))).join(gc).reset_index()
    bins = bins[bins["n"] >= min_sites].copy()
    bins["raw_ratio"] = bins["depth"] / float(bins["depth"].median())
    # collapsed-repeat artefact bins: excess of strongly skewed "hets" and/or inflated depth
    bins["artefact"] = (bins["skew_frac"] > 0.15) | (bins["raw_ratio"] > 1.3)
    ok = ~bins["artefact"]
    coef = np.polyfit(bins.loc[ok, "gc_proxy"], bins.loc[ok, "raw_ratio"], 2) if ok.sum() >= 5 else np.array([0.0, 0.0, 1.0])
    bins["ratio"] = bins["raw_ratio"] - np.polyval(coef, bins["gc_proxy"]) + 1.0
    order = {c: i for i, c in enumerate(AUTOSOMES)}
    bins["ord"] = bins["chrom"].map(order)
    bins = bins.sort_values(["ord", "bin"]).reset_index(drop=True)
    ok = ~bins["artefact"]
    good = bins[ok]
    chrom = good.groupby("chrom", observed=True).agg(bins=("ratio", "size"), depth_ratio=("ratio", "median"),
                                                     raw_depth_ratio=("raw_ratio", "median"),
                                                     baf_dev=("baf_dev", "median")).reset_index()
    chrom["ord"] = chrom["chrom"].map(order)
    chrom = chrom.sort_values("ord").drop(columns="ord").reset_index(drop=True)
    sd_ratio = max(float(chrom["depth_ratio"].std()), 0.01)        # includes residual systematic chromosome effects
    sd_baf = max(float(chrom["baf_dev"].std()), 0.0006)
    chrom["z_depth"] = (chrom["depth_ratio"] - chrom["depth_ratio"].median()) / sd_ratio
    chrom["z_baf"] = (chrom["baf_dev"] - chrom["baf_dev"].median()) / sd_baf
    chrom["flag"] = (chrom["z_depth"].abs() > 4) & (chrom["z_baf"] > 4)
    # a clonal trisomy (monosomy) in a fraction f of cells changes the depth ratio by +f/2 (-f/2)
    limit = dict(sd_depth_ratio=round(sd_ratio, 4), sd_baf_dev=round(sd_baf, 5),
                 gc_correlation=round(float(np.corrcoef(bins.loc[ok, "gc_proxy"], bins.loc[ok, "raw_ratio"])[0, 1]), 3),
                 min_detectable_clone_fraction=round(2 * 4 * sd_ratio, 3))
    return bins, chrom, limit


with stage("QC / sex / autozygosity / aneuploidy screen"):
    sites = load_site_table(S["vcf_raw"])
    S["qc"] = qc_summary(sites)
    roh, f_roh = runs_of_homozygosity(sites)
    S["qc"]["roh_segments_ge_1Mb"] = int(len(roh))
    S["qc"]["f_roh"] = round(f_roh, 4)
    S["qc"]["autozygosity_note"] = ("elevated — parental relatedness likely; homozygous models up-weighted" if f_roh > 0.02
                                    else "no evidence of parental relatedness")
    bins, chrom_tab, limit = aneuploidy_screen(sites)
    S["aneuploidy"] = dict(limit=limit, flagged=chrom_tab[chrom_tab.flag].chrom.tolist(),
                           table=chrom_tab.round(4).to_dict("records"),
                           artefact_bins=int(bins.artefact.sum()), bins=int(len(bins)))
    del sites
    for k, v in S["qc"].items():
        log(f"QC {k}: {v}")
    log(f"aneuploidy screen: flagged chromosomes = {S['aneuploidy']['flagged'] or 'none'}; detection limit = {limit}")

    # figure: two small multiples sharing the genome axis (one measure per panel, never a dual axis)
    off, ticks, x0 = {}, [], 0.0
    for c in AUTOSOMES:
        sub = bins[bins.chrom == c]
        if not len(sub):
            continue
        off[c] = x0
        width = float(sub.bin.max() + 1)
        ticks.append((x0 + width / 2, c))
        x0 += width + 1
    bins["x"] = bins.apply(lambda r: off[r.chrom] + r.bin + 0.5, axis=1)
    fig, axes = plt.subplots(2, 1, figsize=(11, 5.2), sharex=True, gridspec_kw=dict(hspace=0.32))
    for ax, col, title, base in ((axes[0], "ratio", "Depth ratio per 5-Mb bin (heterozygous SNV sites, GC-corrected; 1.0 = autosomal median)", 1.0),
                                 (axes[1], "baf_dev", "Allele-balance dispersion per 5-Mb bin (mean |alt fraction − 0.5|)", None)):
        ok, bad = bins[~bins.artefact], bins[bins.artefact]
        ax.scatter(ok.x, ok[col], s=9, color=VIZ["s1"], linewidths=0, label="bins analysed")
        ax.scatter(bad.x, bad[col].clip(upper=ok[col].max() * 1.6), s=9, color=VIZ["neutral"], linewidths=0,
                   label="excluded: collapsed-repeat artefact (clipped)")
        if base is not None:
            ax.axhline(base, color=VIZ["axis"], linewidth=0.8)
        for c, o in off.items():
            ax.axvline(o - 0.5, color=VIZ["grid"], linewidth=0.6)
        ax.set_title(title)
        ax.grid(axis="x", visible=False)
    axes[0].legend(loc="upper left", ncol=2, handletextpad=0.2, columnspacing=1.2)
    axes[1].set_xticks([t for t, _ in ticks])
    axes[1].set_xticklabels([c for _, c in ticks])
    axes[1].set_xlabel("chromosome")
    save_fig(fig, "aneuploidy_screen")
    md("**Chromosome-level table (artefact bins excluded)**\n\n" + df_md(chrom_tab.round(4), max_rows=25))

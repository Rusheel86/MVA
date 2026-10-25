# %% [markdown]
# ## 12 · Phase of the leading pair: what the data can and cannot say
#
# A compound-heterozygous diagnosis needs the two alleles on opposite haplotypes. With one short-read genome and no
# parents, phase is only observable when both sites sit on the same read pair, or can be chained through heterozygous
# sites that do. This section measures that directly instead of assuming it:
#
# 1. read-backed phase tags from the caller (`PID` / `PGT`);
# 2. the heterozygous sites of the proband around and between the two alleles, and the largest gap that reads would
#    have to bridge;
# 3. if the alleles cannot be linked, an explicit statement of what was inferred instead (a model-based probability,
#    not a measurement) and of the laboratory test that would settle it.

# %%
def phase_dossier(hyp, max_link=600):
    a, b = sorted(hyp["variants"], key=lambda v: int(v["pos"]))
    chrom, lo, hi = a["chrom"], int(a["pos"]), int(b["pos"])
    dist = hi - lo
    raw_chrom = chrom if S.get("vcf_chr_style", True) else chrom[3:]
    out = dict(distance_bp=dist, caller_phase=hyp["phase"], pid_a=a.get("pid"), pid_b=b.get("pid"))
    txt = sh(f"bcftools query -r {raw_chrom}:{max(1, lo - 2000)}-{hi + 2000} -i 'GT=\"het\"' -f '%POS\\t%REF\\t%ALT\\t%FILTER[\\t%GT]\\n' '{S['vcf_raw']}'",
             name="phase_hets", capture=True)
    hets = sorted({int(l.split("\t")[0]) for l in txt.splitlines() if l.strip()} | {lo, hi})
    between = [p for p in hets if lo <= p <= hi]
    gaps = np.diff(between) if len(between) > 1 else np.array([dist])
    out.update(n_het_sites_between=len(between) - 2, largest_gap_bp=int(gaps.max()), max_read_link_bp=max_link,
               nearest_het_to_a_bp=int(min([abs(p - lo) for p in hets if p != lo], default=10 ** 9)),
               nearest_het_to_b_bp=int(min([abs(p - hi) for p in hets if p != hi], default=10 ** 9)))
    if hyp["phase"] in ("cis", "trans"):
        out["verdict"] = f"read-backed phase from the caller: {hyp['phase']}"
        out["phaseable_by_reads"] = True
    elif gaps.max() <= max_link:
        out["verdict"] = ("the two alleles can be chained through intervening heterozygous sites; read-level chaining is attempted "
                          "in the read-level module")
        out["phaseable_by_reads"] = True
    else:
        out["verdict"] = (f"not phaseable with short reads: the alleles are {dist:,} bp apart and the chain of heterozygous sites is "
                          f"broken by a gap of {int(gaps.max()):,} bp (read pairs link at most ~{max_link} bp)")
        out["phaseable_by_reads"] = False
    amplicon = dist + 1000
    cpos = [re.search(r"c\.(\d+)[ACGT]>", v.get("hgvsc") or "") for v in (a, b)]
    cdna = abs(int(cpos[0].group(1)) - int(cpos[1].group(1))) if all(cpos) else None
    out["cdna_distance_nt"] = cdna
    out["confirmation"] = ("testing either parent for the two alleles (the standard clinical route)"
                           + (f"; or a single long-range PCR amplicon of ~{amplicon / 1000:.1f} kb spanning both sites, read on a long-read platform"
                              if amplicon <= 20000 else "; or long-read sequencing of the locus")
                           + (f"; or, because the alleles are only {cdna:,} nt apart in the mRNA, one RT-PCR amplicon from the child's cultured "
                              "cells with an NMD inhibitor (see the Track 2 report)" if cdna and cdna < dist / 3 else ""))
    return out


with stage("phase analysis"):
    TOP = HYPS[0]
    if TOP["model"] != "compound_het":
        S["phase"] = dict(verdict=f"not applicable: the leading hypothesis is {TOP['model']}")
        log(S["phase"]["verdict"])
    else:
        hdr_names = sh(f"bcftools index -s '{S['vcf_raw']}'", name="phase_idx", capture=True)
        S["vcf_chr_style"] = "chr1\t" in hdr_names or hdr_names.startswith("chr")
        S["phase"] = phase_dossier(TOP)
        b = S.get("bystander") or {}
        if b:
            S["phase"]["model_based_probability"] = dict(
                statement=("If the gene is the cause and the stronger allele is pathogenic, the weaker allele is either the second "
                           "causal allele (necessarily in trans) or a coincidental bystander while the true second hit is hidden."),
                bystander_probability=round(b.get("carrier_probability_strict", b["carrier_probability"]), 5),
                posterior_second_allele=b["posterior_second_allele"])
        for k, v in S["phase"].items():
            log(f"phase {k}: {v}")

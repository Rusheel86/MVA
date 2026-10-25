# %% [markdown]
# ## 8 · Normalise and annotate
#
# 1. primary chromosomes, UCSC contig names, multi-allelic sites split, indels left-aligned against the reference
#    (`bcftools norm`); spanning-deletion `*` alleles dropped;
# 2. population frequency (gnomAD v4.1 joint exomes+genomes from a compact `echtvar` archive: AF, grpmax filtering AF,
#    homozygote count; two independent fallbacks);
# 3. ClinVar assertions — looked up genome-wide and *before* any frequency filter, so that known pathogenic
#    regulatory or founder alleles are not lost;
# 4. consequences on MANE Select / MANE Plus Clinical transcripts (`bcftools csq`), an automated PVS1 decision
#    (nonsense-mediated-decay rule) and AlphaMissense for every rare missense variant.
#
# The VCF's own `INFO/AF` is the within-sample allele fraction and is never used as a population frequency.

# %%
CSQ_CLASS = {
    "stop_gained": "stop_gained", "frameshift": "frameshift", "splice_acceptor": "splice_acceptor",
    "splice_donor": "splice_donor", "start_lost": "start_lost", "stop_lost": "stop_lost",
    "inframe_deletion": "inframe_indel", "inframe_insertion": "inframe_indel", "inframe_altering": "inframe_indel",
    "missense": "missense", "splice_region": "splice_region", "synonymous": "synonymous",
    "stop_retained": "synonymous", "start_retained": "synonymous", "coding_sequence": "inframe_indel", "feature_truncation": "frameshift",
    "5_prime_utr": "utr5", "3_prime_utr": "utr3", "intron": "intron", "non_coding": "other",
    "intergenic": "other", "feature_elongation": "inframe_indel",
}
CSQ_RANK = ["stop_gained", "frameshift", "splice_acceptor", "splice_donor", "start_lost", "stop_lost", "inframe_indel",
            "missense", "splice_region", "synonymous", "utr5", "utr3", "intron", "upstream", "downstream", "other"]
CODING_CLASSES = set(CSQ_RANK[:10])
AA3 = dict(A="Ala", R="Arg", N="Asn", D="Asp", C="Cys", Q="Gln", E="Glu", G="Gly", H="His", I="Ile", L="Leu", K="Lys",
           M="Met", F="Phe", P="Pro", S="Ser", T="Thr", W="Trp", Y="Tyr", V="Val")
AA3["*"] = "Ter"


def parse_bcsq(bcsq):
    """INFO/BCSQ -> {gene: (class, transcript, aa_change)} keeping the most severe consequence per gene."""
    out = {}
    if not bcsq or bcsq == ".":
        return out
    for entry in bcsq.split(","):
        f = entry.split("|")
        if len(f) < 3 or entry.startswith("@"):
            continue
        terms = [CSQ_CLASS.get(t.lstrip("*"), "other") for t in f[0].split("&")]
        cls = min(terms, key=CSQ_RANK.index)
        gene, tid = f[1], f[2]
        aa = f[5] if len(f) > 5 else ""
        if gene not in out or CSQ_RANK.index(cls) < CSQ_RANK.index(out[gene][0]):
            out[gene] = (cls, tid, aa)
    return out


def protein_hgvs(aa):
    """bcftools csq amino-acid string ('737L>737*', '1002N>1002K') -> ('p.Leu737Ter', 737)."""
    m = re.match(r"^(\d+)([A-Z*]+)>(\d+)([A-Z*]+)$", aa or "")
    if not m:
        m2 = re.match(r"^(\d+)", aa or "")
        return "", (int(m2.group(1)) if m2 else None)
    pos, ref, alt = int(m.group(1)), m.group(2), m.group(4)
    if len(ref) == 1 and len(alt) == 1:
        return f"p.{AA3.get(ref, ref)}{pos}{AA3.get(alt, alt)}", pos
    return f"p.{AA3.get(ref[0], ref[0])}{pos}delins", pos


def nmd_status(gene, pos):
    """(escapes_nmd, fraction_of_protein_lost) for a truncating variant, from the MANE CDS structure:
    NMD is not predicted when the stop falls in the last coding exon or within 50 nt of the end of the penultimate one."""
    tids = GENE_TX.get(gene) or []
    t = TX.get(tids[0]) if tids else None
    if not t or not t["cds"]:
        return None, None
    plus, cum = t["strand"] == "+", 0
    for i, (s, e) in enumerate(t["cds"]):
        if s - 3 <= pos <= e + 3:
            p = min(max(pos, s), e)
            off, dist3 = (p - s, e - p) if plus else (e - p, p - s)
            last = len(t["cds"]) - 1
            escape = last == 0 or i == last or (i == last - 1 and dist3 < 50)
            return escape, round(1.0 - (cum + off) / max(t["cds_len"], 1), 3)
        cum += e - s + 1
    return None, None


def _has(hdr, kind, tag):
    return f"##{kind}=<ID={tag}," in hdr


def annotate_vcf(vcf_in, tag, spikes=None, sex="unknown"):
    """Raw single-sample VCF -> (candidate table, non-coding pool, annotated BCF path)."""
    d = WORK / "data"
    hdr = sh(f"bcftools view -h '{vcf_in}'", name=f"{tag}_hdr", capture=True)
    idx_names = [l.split("\t")[0] for l in sh(f"bcftools index -s '{vcf_in}'", name=f"{tag}_idx", capture=True).splitlines() if l.strip()]
    idx_set = set(idx_names)
    if "chr1" in idx_set or "chrX" in idx_set:
        chr_style = True
    elif "1" in idx_set or "X" in idx_set:
        chr_style = False
    else:
        raise RuntimeError(f"[{tag}] unrecognised contig naming: {idx_names[:8]}")
    primary_in = [c for c in (PRIMARY if chr_style else [p[3:] for p in PRIMARY]) if c in idx_set]
    if not primary_in:
        raise RuntimeError(f"[{tag}] no primary chromosomes with records")
    rename = ""
    if not chr_style:
        names = set(re.findall(r"##contig=<ID=([^,>]+)", hdr)) | set(idx_names)
        mapf = d / f"{tag}.chrmap.txt"
        mapf.write_text("".join(f"{c}\tchr{'M' if c == 'MT' else c}\n" for c in sorted(names)))
        rename = f"| bcftools annotate --rename-chrs '{mapf}' -Ou"
    norm = d / f"{tag}.norm.bcf"
    sh(f"bcftools view -r {','.join(primary_in)} -Ou '{vcf_in}' {rename} "
       f"| bcftools norm -m -any -f '{RES['ref']}' -c w -Ou "
       f"| bcftools view -e 'ALT=\"*\"' --threads {NCPU} -Ob -o '{norm}' && bcftools index -f '{norm}'", name=f"{tag}_norm")
    n_norm = int(sh(f"bcftools index -n '{norm}'", name=f"{tag}_n", capture=True).strip())
    log(f"[{tag}] normalised biallelic records on primary chromosomes: {n_norm:,}")

    if spikes:      # synthetic benchmark alleles (public ClinVar variants), tagged in the ID column
        order = {c: i for i, c in enumerate(PRIMARY)}
        rows = sorted(((c, p, r, a, f"SPIKE_{cs['case']}_{k}") for cs in spikes for k, (c, p, r, a) in enumerate(cs["variants"])),
                      key=lambda x: (order[x[0]], x[1]))
        spk = d / f"{tag}.spikes.vcf"
        with open(spk, "w") as f:
            f.write(sh(f"bcftools view -h '{norm}'", name=f"{tag}_hdr2", capture=True))
            for c, p, r, a, vid in rows:
                f.write(f"{c}\t{p}\t{vid}\t{r}\t{a}\t500\tPASS\t.\tGT:AD:DP:GQ\t0/1:20,20:40:99\n")
        merged = d / f"{tag}.spiked.bcf"
        sh(f"bgzip -f '{spk}' && bcftools index -f '{spk}.gz' && bcftools concat -a -Ou '{norm}' '{spk}.gz' "
           f"| bcftools sort -T '{WORK / 'tmp'}' -m 2G -Ob -o '{merged}' && bcftools index -f '{merged}'", name=f"{tag}_spike")
        norm.unlink()
        norm = merged
        n_norm += len(rows)

    # ---- population frequency: first source that works (fallback archives are downloaded lazily) --------------------
    gn, af_source, rare_expr = d / f"{tag}.gn.bcf", None, None
    for src in AF_SOURCES:
        try:
            if src["name"] == "echtvar":
                if "af_tool" not in RES or "af_data" not in RES:
                    raise RuntimeError("echtvar archive not downloaded")
                raw_gn = d / f"{tag}.gn.raw.bcf"
                sh(f"'{RES['af_tool']}' anno -e '{RES['af_data']}' '{norm}' '{raw_gn}' && bcftools view --threads {NCPU} -Ob -o '{gn}' '{raw_gn}' "
                   f"&& rm -f '{raw_gn}' && bcftools index -f '{gn}'", name=f"{tag}_echtvar")
                need, tagname = "gnomad_faf95max", "gnomad_faf95max"
            elif src["name"] == "slivar":
                tool = fetch(src["tool"], WORK / "bin" / "slivar", min_bytes=3_000_000)
                Path(tool).chmod(0o755)
                data = fetch(src["data"], WORK / "ref" / "af_slivar.zip", min_bytes=src["min_bytes"])
                sh(f"'{tool}' expr --vcf '{norm}' -g '{data}' --info 'true' | bcftools view --threads {NCPU} -Ob -o '{gn}' "
                   f"&& bcftools index -f '{gn}'", name=f"{tag}_slivar", env={"SLIVAR_QUIET": "1"})
                need, tagname = "gnomad_popmax_af", "gnomad_popmax_af"
            else:
                data = fetch(src["data"], WORK / "ref" / "af-only-gnomad.hg38.vcf.gz", min_bytes=src["min_bytes"])
                fetch(src["data"] + ".tbi", WORK / "ref" / "af-only-gnomad.hg38.vcf.gz.tbi", min_bytes=1_000_000)
                sh(f"bcftools annotate -a '{data}' -c 'INFO/gnomad_popmax_af:=INFO/AF' --threads {NCPU} -Ob -o '{gn}' '{norm}' "
                   f"&& bcftools index -f '{gn}'", name=f"{tag}_afonly")
                need, tagname = "gnomad_popmax_af", "gnomad_popmax_af"
            n_gn = int(sh(f"bcftools index -n '{gn}'", name=f"{tag}_n2", capture=True).strip())
            if n_gn < 0.98 * n_norm or not _has(sh(f"bcftools view -h '{gn}'", name=f"{tag}_h3", capture=True), "INFO", need):
                raise RuntimeError(f"annotation incomplete ({n_gn:,} of {n_norm:,} records)")
            rare_expr = f'(INFO/{tagname}<{2 * CFG["MAX_AF"]} || INFO/{tagname}=".")'
            af_source = src["label"]
            break
        except Exception as e:      # noqa: BLE001
            log(f"[{tag}] population-frequency source '{src['name']}' failed: {str(e)[:220]}")
    if af_source is None:
        raise RuntimeError("no population-frequency resource could be applied")
    norm.unlink(missing_ok=True)

    # ---- ClinVar (before any frequency filter), rarity filter, consequences ---------------------------------------
    keep = rare_expr + ' || INFO/CLNSIG~"athogenic" || ID~"^SPIKE"'
    ann = d / f"{tag}.ann.bcf"
    sh(f"bcftools annotate -a '{RES['clinvar_chr']}' -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/GENEINFO -Ou '{gn}' "
       f"| bcftools view -i '{keep}' -Ou "
       f"| bcftools csq -f '{RES['ref']}' -g '{RES['csq_gff']}' --local-csq --force --ncsq 64 --phase a -Ob -o '{ann}' "
       f"&& bcftools index -f '{ann}'", name=f"{tag}_annot")
    gn.unlink(missing_ok=True)
    Path(str(gn) + ".csi").unlink(missing_ok=True)

    h = sh(f"bcftools view -h '{ann}'", name=f"{tag}_h4", capture=True)
    info_tags = [t for t in ("gnomad_faf95max", "gnomad_AF", "gnomad_AC", "gnomad_HomAlt", "gnomad_popmax_af", "gnomad_nhomalt",
                             "CLNSIG", "CLNREVSTAT", "CLNDN", "GENEINFO", "BCSQ") if _has(h, "INFO", t)]
    fmt_tags = [t for t in ("GT", "AD", "DP", "GQ", "PGT", "PID") if _has(h, "FORMAT", t)]
    fmt = ("%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\t%QUAL\\t%FILTER" + "".join(f"\\t%INFO/{t}" for t in info_tags)
           + "[" + "".join(f"\\t%{t}" for t in fmt_tags) + "]\\n")
    tsv = d / f"{tag}.ann.tsv"
    sh(f"bcftools query -f '{fmt}' '{ann}' > '{tsv}'", name=f"{tag}_query")
    cols = ["chrom", "pos", "ref", "alt", "rsid", "qual", "filt"] + [t.lower() for t in info_tags] + [t.lower() for t in fmt_tags]
    raw = pd.read_csv(tsv, sep="\t", header=None, names=cols, dtype=str, keep_default_na=False, quoting=3)
    tsv.unlink(missing_ok=True)
    for c in ("gnomad_faf95max", "gnomad_af", "gnomad_ac", "gnomad_homalt", "gnomad_popmax_af", "gnomad_nhomalt", "clnsig",
              "clnrevstat", "clndn", "geneinfo", "bcsq", "gt", "ad", "dp", "gq", "pgt", "pid"):
        if c not in raw.columns:
            raw[c] = "."
    log(f"[{tag}] rare or ClinVar-pathogenic records: {len(raw):,} ({af_source})")

    # ---- one row per (variant, gene) ---------------------------------------------------------------------------------
    recs = []
    for r in raw.itertuples(index=False):
        alleles = re.split(r"[/|]", r.gt)
        n_alt = sum(1 for x in alleles if x == "1")
        if n_alt == 0:
            continue
        zyg = "hom" if n_alt == len(alleles) else "het"
        male = sex.startswith("male")
        x_nonpar = r.chrom == "chrY" or (r.chrom == "chrX" and not any(a <= int(r.pos) <= b for a, b in PAR_X))
        if zyg == "hom" and (len(alleles) == 1 or (male and x_nonpar)):
            zyg = "hemi"
        elif zyg == "het" and male and x_nonpar:
            zyg = "het_xy"          # implausible in a male: single-allele model only, never paired
        ad = (r.ad.split(",") + [".", "."])[:2]
        # unified population frequency: the larger of the global AF and the grpmax filtering AF (echtvar), or the
        # max population AF (slivar; "absent" is encoded as -1), or the global AF (af-only). None = not in gnomAD.
        f95, gaf, pmx = _num(r.gnomad_faf95max), _num(r.gnomad_af), _num(r.gnomad_popmax_af)
        af = max([x for x in (f95, gaf, pmx) if x is not None and x > 0], default=None)
        nhom = _num(r.gnomad_homalt) if _num(r.gnomad_homalt) is not None else _num(r.gnomad_nhomalt)
        base = dict(chrom=r.chrom, pos=int(r.pos), ref=r.ref, alt=r.alt, rsid=r.rsid, qual=_num(r.qual), filt=r.filt,
                    zyg=zyg, ad_ref=_num(ad[0]), ad_alt=_num(ad[1]), dp=_num(r.dp), gq=_num(r.gq),
                    pgt=r.pgt, pid=r.pid, af=af, nhomalt=nhom, gnomad_ac=_num(r.gnomad_ac),
                    clnsig=r.clnsig if r.clnsig != "." else "", clnrevstat=r.clnrevstat if r.clnrevstat != "." else "",
                    clndn=(r.clndn if r.clndn != "." else "")[:120])
        genes = parse_bcsq(r.bcsq)
        if not genes and base["clnsig"] and r.geneinfo != ".":      # e.g. a known pathogenic regulatory allele
            genes = {r.geneinfo.split(":")[0]: ("other", "", "")}
        for g, (cls, tid, aa) in genes.items():
            hg, aa_pos = protein_hgvs(aa)
            recs.append(dict(base, gene=g, transcript=tid, csq=cls, aa_change=aa, hgvsp=hg, aa_pos=aa_pos))
        if not genes:
            recs.append(dict(base, gene="", transcript="", csq="other", aa_change="", hgvsp="", aa_pos=None))
    allv = pd.DataFrame(recs)
    allv["nmd_escape"], allv["frac_lost"], allv["am"], allv["spliceai"], allv["revel"], allv["cadd"] = None, np.nan, np.nan, np.nan, np.nan, np.nan
    lof_mask = allv.csq.isin(LOF_CSQ)
    for i in allv.index[lof_mask]:
        esc, frac = nmd_status(allv.at[i, "gene"], allv.at[i, "pos"])
        allv.at[i, "nmd_escape"], allv.at[i, "frac_lost"] = esc, (np.nan if frac is None else frac)

    # ---- AlphaMissense for every rare missense variant (tabix lookups in the public BGZF table) ------------------
    mis = allv[allv.csq == "missense"]
    if len(mis) and "alphamissense" in RES:
        try:
            reg = d / f"{tag}.am.regions.tsv"
            mis[["chrom", "pos"]].drop_duplicates().assign(p2=lambda x: x.pos).to_csv(reg, sep="\t", header=False, index=False)
            out = sh(f"tabix -R '{reg}' '{RES['alphamissense']}'", name=f"{tag}_am", capture=True)
            am = {}
            for line in out.splitlines():
                c = line.split("\t")
                if len(c) >= 9:
                    am[(c[0], int(c[1]), c[2], c[3])] = float(c[8])
            keys = list(zip(allv.chrom, allv.pos, allv.ref, allv.alt))
            allv.loc[mis.index, "am"] = [am.get(keys[i], np.nan) for i in mis.index]
            log(f"[{tag}] AlphaMissense scores attached to {int(allv.loc[mis.index, 'am'].notna().sum()):,} of {len(mis):,} rare missense variants")
        except Exception as e:      # noqa: BLE001
            log(f"[{tag}] AlphaMissense lookup failed: {str(e)[:200]}")

    is_plp = allv.clnsig.map(lambda s: clinvar_class(s) == "P/LP")
    rare = allv.af.isna() | (allv.af < CFG["MAX_AF"]) | is_plp
    cand = allv[rare & (allv.gene != "") & (allv.csq.isin(CODING_CLASSES) | is_plp)].reset_index(drop=True)
    pool = allv[rare & ~(allv.csq.isin(CODING_CLASSES) | is_plp)].reset_index(drop=True)      # non-coding, for second-hit rescue
    log(f"[{tag}] candidate table: {len(cand):,} coding/splice or ClinVar-P/LP rows in {cand.gene.nunique():,} genes; "
        f"non-coding pool: {len(pool):,} rows")
    return cand, pool, str(ann), af_source


with stage("normalise + annotate", required=True):
    _sex = (S.get("qc") or {}).get("inferred_sex", "unknown")
    CAND_ALL, POOL_ALL, ANN_BCF, S["af_source"] = annotate_vcf(S["vcf_raw"], "case", spikes=SPIKE_CASES if DEMO else None, sex=_sex)
    _is_spike = CAND_ALL.rsid.str.startswith("SPIKE_")
    if DEMO:        # the demo proband = HG002 background + the two alleles of spike case 0
        CAND = CAND_ALL[~_is_spike | CAND_ALL.rsid.isin(["SPIKE_0_0", "SPIKE_0_1"])].reset_index(drop=True)
    else:
        CAND = CAND_ALL[~_is_spike].reset_index(drop=True)
        assert not _is_spike.any(), "synthetic alleles must never enter the challenge candidate table"
    POOL = POOL_ALL[~POOL_ALL.rsid.str.startswith("SPIKE_")].reset_index(drop=True)
    S["n_candidates"], S["n_candidate_genes"] = int(len(CAND)), int(CAND.gene.nunique())
    md("**Candidate rows by consequence and zygosity**\n\n" + df_md(CAND.groupby(["csq", "zyg"]).size().unstack(fill_value=0).reset_index()))

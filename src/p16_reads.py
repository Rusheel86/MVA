# %% [markdown]
# ## 16 · Read-level evidence from the raw FASTQ (optional; the long stage)
#
# Aligning 85 GB of reads genome-wide is not possible on this machine, and is not needed. Instead:
#
# 1. **Bait** — the leading gene ± 50 kb, the next genes ± 5 kb, and ~500 random 2-kb control windows (autosomes and, as a
#    built-in single-copy control in a male, chromosome X).
# 2. **Unique-minimizer recruitment** — a `deacon` minimizer index of the bait, minus every minimizer that also occurs
#    elsewhere in GRCh38, so reads from repeats (Alu, L1) are not pulled in. Each lane is downloaded, streamed through the
#    recruiter (a pair is kept if either mate hits) and deleted.
# 3. **Re-alignment of the recruited reads to the whole genome** (`minimap2 -ax sr`, duplicates marked), so reads that
#    belong elsewhere receive low mapping quality instead of piling up on the bait.
# 4. **Analyses** — independent allele counts (per strand, per lane) at the leading alleles; re-calling with a second
#    caller and concordance with the provided VCF; read-pair links between consecutive heterozygous sites between the two
#    alleles (empirical test of phaseability); depth along the gene relative to autosomal controls (chrX as positive
#    control for a one-copy change); clipped-read and discordant-pair clusters (structural-variant screen).
#
# The module stops starting new lanes after `CFG["READ_LEVEL_MAX_HOURS"]` of total run time; one lane (~10×) is already
# enough for allele-level validation. In demo mode it runs on the first 2 M read pairs of a public HG002 run, as a
# toolchain test only.

# %%
READS = WORK / "reads"
HG002_FQ = "https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x/HG002.novaseq.pcr-free.30x.R{r}.fastq.gz"


def ngs_env():
    """deacon (+ samtools / minimap2 already on PATH) in its own micromamba environment."""
    env = WORK / "env_reads"
    if not (env / "bin" / "deacon").exists():
        if not MICROMAMBA.exists():
            fetch("https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64", MICROMAMBA, min_bytes=1_000_000)
            MICROMAMBA.chmod(0o755)
        sh(f"'{MICROMAMBA}' create -y -q -p '{env}' -c conda-forge -c bioconda deacon", name="mamba_deacon", timeout=1800)
    return env / "bin" / "deacon"


def make_bait(genes, fa, n_auto=400, n_x=100, seed=3):
    regs = []
    for i, g in enumerate(genes):
        if g in GENE_POS:
            c, s, e, _, _ = GENE_POS[g]
            f = 50_000 if i == 0 else 5_000
            regs.append((c, max(1, s - f), min(CHROM_LEN[c], e + f), f"target_{g}"))
    rng = np.random.default_rng(seed)
    autos = [c for c in PRIMARY[:22]]
    w = np.array([CHROM_LEN[c] for c in autos], dtype=float)
    ctrl = []
    while len(ctrl) < n_auto + n_x:
        if len(ctrl) < n_auto:
            c = autos[int(rng.choice(len(autos), p=w / w.sum()))]
            s = int(rng.integers(5_000_000, CHROM_LEN[c] - 5_000_000))
            tag = "ctrl_auto"
        else:
            c, s, tag = "chrX", int(rng.integers(3_000_000, 150_000_000)), "ctrl_X"
        if any(c == r[0] and s < r[2] and s + 2000 > r[1] for r in regs):
            continue
        seq = str(fa[c][s - 1:s + 1999])
        if seq.count("N") > 20 or len(seq) < 2000:
            continue
        ctrl.append((c, s, s + 1999, tag))
    return regs + ctrl


def write_masked_genome_script(bait, path):
    """Stream GRCh38 with the bait intervals replaced by N (input for subtracting non-unique minimizers)."""
    by = collections.defaultdict(list)
    for c, s, e, _ in bait:
        by[c].append((s, e))
    path.write_text(textwrap.dedent(f'''
        import sys
        from pyfaidx import Fasta
        bait = {dict(by)!r}
        fa = Fasta({RES["ref"]!r}, as_raw=True, sequence_always_upper=True)
        out = sys.stdout
        for name in fa.keys():
            seq = fa[name][:]
            for s, e in sorted(bait.get(name, [])):
                seq = seq[:s - 1] + "N" * (e - s + 1) + seq[e:]
            out.write(">" + name + "\\n")
            for i in range(0, len(seq), 1000000):
                out.write(seq[i:i + 1000000] + "\\n")
    '''))


def fastq_sizes():
    from huggingface_hub import HfApi
    files = [CFG["FASTQ_PATTERN"].format(lane=l, read=r) for l in CFG["FASTQ_LANES"] for r in (1, 2)]
    info = HfApi(token=HF_TOKEN).get_paths_info(CFG["HF_DATASET"], files, repo_type="dataset")
    return {f.path: f.size for f in info}


def recruit_lane(deacon, idx, lane, sizes, fq_dir, out_dir):
    """Download both mates of one lane (or stream them when scratch space is short), recruit, delete."""
    fns = [CFG["FASTQ_PATTERN"].format(lane=lane, read=r) for r in (1, 2)]
    o1, o2, js = out_dir / f"L{lane}.R1.fq.gz", out_dir / f"L{lane}.R2.fq.gz", out_dir / f"L{lane}.json"
    need = sum(sizes.get(f, 11e9) for f in fns) / 1e9
    stream = free_gb(WORK) < need + 8
    t0 = time.time()
    if not stream:
        from huggingface_hub import hf_hub_download
        with ThreadPoolExecutor(2) as ex:
            paths = list(ex.map(lambda f: hf_hub_download(CFG["HF_DATASET"], f, repo_type="dataset", local_dir=str(fq_dir), token=HF_TOKEN), fns))
        for p, f in zip(paths, fns):
            if f in sizes and Path(p).stat().st_size != sizes[f]:
                raise RuntimeError(f"size mismatch for {f}")
        t_dl = time.time() - t0
        sh(f"'{deacon}' filter --threads {NCPU} -a 2 -s '{js}' '{idx}' '{paths[0]}' '{paths[1]}' -o '{o1}' -O '{o2}'", name=f"deacon_L{lane}", timeout=4 * 3600)
        for p in paths:
            Path(p).unlink(missing_ok=True)
        shutil.rmtree(fq_dir / ".cache", ignore_errors=True)
    else:
        base = f"https://huggingface.co/datasets/{CFG['HF_DATASET']}/resolve/main/"
        p1, p2 = fq_dir / f"L{lane}_R1.fq.gz", fq_dir / f"L{lane}_R2.fq.gz"
        sh(f"""rm -f '{p1}' '{p2}'; mkfifo '{p1}' '{p2}'
               curl -sSfL --retry 3 -H "Authorization: Bearer $HF_TOKEN" '{base}{fns[0]}' > '{p1}' & c1=$!
               curl -sSfL --retry 3 -H "Authorization: Bearer $HF_TOKEN" '{base}{fns[1]}' > '{p2}' & c2=$!
               '{deacon}' filter --threads {NCPU} -a 2 -s '{js}' '{idx}' '{p1}' '{p2}' -o '{o1}' -O '{o2}'
               wait $c1; wait $c2; rm -f '{p1}' '{p2}'""", name=f"deacon_stream_L{lane}", timeout=4 * 3600)
        t_dl = float("nan")
    summ = json.loads(js.read_text()) if js.exists() else {}
    return dict(lane=lane, mode="stream" if stream else "download", minutes=round((time.time() - t0) / 60, 1),
                download_minutes=round(t_dl / 60, 1) if t_dl == t_dl else None,
                **{k: v for k, v in summ.items() if isinstance(v, (int, float)) and k.startswith(("seqs", "bp", "seconds", "time"))})


def demo_reads(deacon, idx, out_dir, n_pairs=2_000_000):
    fq = {}
    for r in (1, 2):
        fq[r] = READS / f"demo.R{r}.fq"
        sh(f"set +o pipefail; curl -sSfL '{HG002_FQ.format(r=r)}' | gzip -dc 2>/dev/null | head -n {4 * n_pairs} > '{fq[r]}'; true",
           name=f"demo_fq_R{r}", timeout=1800)
    js = out_dir / "L1.json"
    sh(f"'{deacon}' filter --threads {NCPU} -a 2 -s '{js}' '{idx}' '{fq[1]}' '{fq[2]}' -o '{out_dir / 'L1.R1.fq.gz'}' -O '{out_dir / 'L1.R2.fq.gz'}'",
       name="deacon_demo")
    for p in fq.values():
        p.unlink(missing_ok=True)
    summ = json.loads(js.read_text()) if js.exists() else {}
    return [dict(lane="demo (HG002, 2M pairs)", mode="stream",
                 **{k: v for k, v in summ.items() if isinstance(v, (int, float)) and k.startswith(("seqs", "bp", "seconds", "time"))})]


def allele_counts(bam, v, min_mapq=20, min_bq=20):
    """Per-strand and per-lane read counts supporting REF / ALT at an SNV (fragment-level: one count per read pair)."""
    import pysam
    out = dict(ref_fwd=0, ref_rev=0, alt_fwd=0, alt_rev=0, other=0, lanes=collections.Counter())
    seen = set()
    pos0 = int(v["pos"]) - 1
    with pysam.AlignmentFile(str(bam)) as f:
        for col in f.pileup(v["chrom"], pos0, pos0 + 1, truncate=True, min_mapping_quality=min_mapq, min_base_quality=min_bq,
                            stepper="nofilter", ignore_orphans=False):
            for pr in col.pileups:
                a = pr.alignment
                if pr.is_del or pr.is_refskip or a.is_duplicate or a.is_secondary or a.is_supplementary or a.query_name in seen:
                    continue
                seen.add(a.query_name)
                b = a.query_sequence[pr.query_position]
                strand = "rev" if a.is_reverse else "fwd"
                if b == v["ref"]:
                    out[f"ref_{strand}"] += 1
                elif b == v["alt"]:
                    out[f"alt_{strand}"] += 1
                    parts = a.query_name.split(":")
                    out["lanes"][parts[3] if len(parts) > 4 else "?"] += 1
                else:
                    out["other"] += 1
    out["lanes"] = dict(out["lanes"])
    return out


def fragment_links(bam, sites, min_mapq=20):
    """Number of read pairs that cover both members of each consecutive pair of heterozygous SNV sites."""
    import pysam
    cover = collections.defaultdict(set)
    with pysam.AlignmentFile(str(bam)) as f:
        for p in sites:
            for a in f.fetch(sites[0][0], p[1] - 1, p[1]):
                if a.mapping_quality >= min_mapq and not a.is_duplicate and not a.is_secondary and not a.is_supplementary:
                    if a.reference_start <= p[1] - 1 < a.reference_end:
                        cover[p[1]].add(a.query_name)
    return [dict(from_pos=a[1], to_pos=b[1], gap_bp=b[1] - a[1], linking_pairs=len(cover[a[1]] & cover[b[1]]))
            for a, b in zip(sites, sites[1:])]


def sv_screen(bam, chrom, start, end, min_mapq=20):
    import pysam
    clips, disc = collections.Counter(), collections.Counter()
    with pysam.AlignmentFile(str(bam)) as f:
        for a in f.fetch(chrom, start, end):
            if a.mapping_quality < min_mapq or a.is_duplicate or a.is_secondary or a.is_supplementary:
                continue
            cig = a.cigartuples or []
            if cig and cig[0][0] == 4 and cig[0][1] >= 20:
                clips[a.reference_start // 10 * 10] += 1
            if cig and cig[-1][0] == 4 and cig[-1][1] >= 20:
                clips[a.reference_end // 10 * 10] += 1
            if a.is_paired and not a.mate_is_unmapped and (a.reference_name != a.next_reference_name or abs(a.template_length) > 2000):
                disc[a.reference_start // 500 * 500] += 1
    return ([dict(pos=p, clipped_reads=n) for p, n in sorted(clips.items()) if n >= 5],
            [dict(bin_start=p, discordant_pairs=n) for p, n in sorted(disc.items()) if n >= 4])


with stage("read-level evidence from the raw FASTQ"):
    if not CFG["RUN_READ_LEVEL"]:
        raise RuntimeError("disabled in CFG")
    import pysam
    from pyfaidx import Fasta
    for d_ in (READS, READS / "fq", READS / "rec"):
        d_.mkdir(parents=True, exist_ok=True)
    # free scratch space: the genome-wide frequency archive is no longer needed
    for k in ("af_data",):
        if k in RES and Path(RES[k]).exists():
            Path(RES[k]).unlink()
    TOPH = HYPS[0]
    genes = list(dict.fromkeys([TOPH["gene"]] + [h["gene"] for h in HYPS[:20]]))[:3]
    deacon = ngs_env()
    fa = Fasta(RES["ref"], as_raw=True, sequence_always_upper=True)
    bait = make_bait(genes, fa)
    with open(READS / "bait.fa", "w") as f:
        for c, s, e, tag in bait:
            f.write(f">{tag}_{c}_{s}\n{str(fa[c][s - 1:e])}\n")
    pd.DataFrame(bait, columns=["chrom", "start", "end", "tag"]).to_csv(READS / "bait.tsv", sep="\t", index=False)
    write_masked_genome_script(bait, READS / "mask.py")
    sh(f"'{deacon}' index build '{READS / 'bait.fa'}' > '{READS / 'bait.idx'}'", name="deacon_index")
    sh(f"'{sys.executable}' '{READS / 'mask.py'}' | '{deacon}' index diff '{READS / 'bait.idx'}' - > '{READS / 'bait.uniq.idx'}'",
       name="deacon_diff", timeout=3600)
    info = sh(f"'{deacon}' index info '{READS / 'bait.uniq.idx'}' 2>&1 || true", name="deacon_info", capture=True)
    log(f"bait: {len(bait)} regions, {sum(e - s + 1 for _, s, e, _ in bait) / 1e6:.2f} Mb; unique-minimizer index: {info.strip()[:200]}")
    idx = READS / "bait.uniq.idx"
    lanes_done = []
    if DEMO:
        lanes_done = demo_reads(deacon, idx, READS / "rec")
    else:
        sizes = fastq_sizes()
        for lane in CFG["FASTQ_LANES"]:
            if elapsed_h() > CFG["READ_LEVEL_MAX_HOURS"]:
                log(f"time budget reached ({elapsed_h():.1f} h): no further lanes")
                break
            try:
                r = recruit_lane(deacon, idx, lane, sizes, READS / "fq", READS / "rec")
                lanes_done.append(r)
                log(f"lane {lane}: {r}")
            except Exception as e:      # noqa: BLE001
                log(f"lane {lane} failed ({str(e)[:200]}); continuing")
                for p in (READS / "fq").glob("*"):
                    p.unlink() if p.is_file() else shutil.rmtree(p, ignore_errors=True)
    if not list((READS / "rec").glob("L*.R1.fq.gz")):
        raise RuntimeError("no lane could be recruited")
    # ---- align the recruited reads to the whole genome --------------------------------------------------------------
    bam = READS / "mini.bam"
    sh(f"cat {READS / 'rec'}/L*.R1.fq.gz > {READS / 'all.R1.fq.gz'} && cat {READS / 'rec'}/L*.R2.fq.gz > {READS / 'all.R2.fq.gz'} && "
       f"minimap2 -ax sr -t {NCPU} -R '@RG\\tID:proband\\tSM:proband\\tPL:ILLUMINA' '{RES['ref']}' {READS / 'all.R1.fq.gz'} {READS / 'all.R2.fq.gz'} "
       f"| samtools fixmate -m -u - - | samtools sort -m 1500M -@ 2 -T '{WORK / 'tmp' / 'srt'}' -u - "
       f"| samtools markdup -@ 2 -d 2500 - '{bam}' && samtools index '{bam}'", name="align", timeout=3 * 3600)
    flag = sh(f"samtools flagstat '{bam}'", name="flagstat", capture=True)
    n_reads = int(flag.split()[0])
    # ---- depth: bait targets vs control windows ---------------------------------------------------------------------
    bed = READS / "bins.bed"
    rows = []
    for c, s, e, tag in bait:
        if tag.startswith("target"):
            for b in range(s, e, 500):
                rows.append((c, b - 1, min(b + 499, e), tag))
        else:
            rows.append((c, s - 1, e, tag))
    pd.DataFrame(rows).to_csv(bed, sep="\t", header=False, index=False)
    cov = pd.read_csv(io.StringIO(sh(f"samtools bedcov -Q 20 -j '{bed}' '{bam}'", name="bedcov", capture=True)), sep="\t", header=None,
                      names=["chrom", "start", "end", "tag", "sum"])
    cov["depth"] = cov["sum"] / (cov.end - cov.start)
    auto_med = float(cov[cov.tag == "ctrl_auto"].depth.median()) or float("nan")
    x_ratio = float(cov[cov.tag == "ctrl_X"].depth.median()) / auto_med if auto_med == auto_med and auto_med > 0 else float("nan")
    tg = cov[cov.tag == f"target_{TOPH['gene']}"].copy()
    tg["ratio"] = tg.depth / auto_med
    gs, ge = GENE_POS[TOPH["gene"]][1], GENE_POS[TOPH["gene"]][2]
    in_gene = tg[(tg.end > gs) & (tg.start < ge)]
    low = in_gene[in_gene.ratio < 0.65]
    runs = (low.start.diff().fillna(500) > 500).cumsum()
    del_runs = [dict(start=int(g.start.min()), end=int(g.end.max()), bins=len(g), mean_ratio=round(float(g.ratio.mean()), 2))
                for _, g in low.groupby(runs) if len(g) >= 3]
    clip_cl, disc_cl = sv_screen(bam, GENE_POS[TOPH["gene"]][0], gs, ge)
    # ---- leading alleles: independent counts and re-calling -----------------------------------------------------------
    allele_rows = []
    for v in TOPH["variants"]:
        if len(v["ref"]) == 1 and len(v["alt"]) == 1:
            ac = allele_counts(bam, v)
            n_ref, n_alt = ac["ref_fwd"] + ac["ref_rev"], ac["alt_fwd"] + ac["alt_rev"]
            from scipy.stats import binomtest
            allele_rows.append(dict(variant=describe_variant(v), ref_reads=n_ref, alt_reads=n_alt,
                                    alt_fraction=round(n_alt / max(n_ref + n_alt, 1), 3),
                                    alt_strands=f"+{ac['alt_fwd']}/-{ac['alt_rev']}", alt_by_lane=ac["lanes"],
                                    het_binom_p=round(binomtest(n_alt, max(n_ref + n_alt, 1), 0.5).pvalue, 3) if n_ref + n_alt else None,
                                    vcf_AD="/".join(fmt(int(x)) if isinstance(x, (int, float)) and x == x else fmt(x) for x in (v.get('ad_ref'), v.get('ad_alt')))))
    tgt_bed = READS / "targets.bed"
    pd.DataFrame([(c, s - 1, e) for c, s, e, t in bait if t.startswith("target")]).to_csv(tgt_bed, sep="\t", header=False, index=False)
    sh(f"bcftools mpileup -f '{RES['ref']}' -R '{tgt_bed}' -q 20 -Q 20 -a AD,DP -d 10000 -Ou '{bam}' "
       f"| bcftools call -mv -Oz -o '{READS / 'recall.vcf.gz'}' && bcftools index -f '{READS / 'recall.vcf.gz'}'", name="recall")
    recall = {(c[0], int(c[1]), c[2], c[3]): c[4] for c in
              (l.split("\t") for l in sh(f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' '{READS / 'recall.vcf.gz'}'",
                                          name="recall_q", capture=True).splitlines() if l)}
    _names = {l.split("\t")[0] for l in sh(f"bcftools index -s '{S['vcf_raw']}'", name="raw_idx", capture=True).splitlines() if l.strip()}
    raw_chr = "chr1" in _names or "chrX" in _names
    raw_bed = READS / "targets.raw.bed"
    pd.DataFrame([(c if raw_chr else c[3:], s - 1, e) for c, s, e, t in bait if t.startswith("target")]).to_csv(
        raw_bed, sep="\t", header=False, index=False)
    orig = {(("" if raw_chr else "chr") + c[0], int(c[1]), c[2], c[3]): c[4] for c in
            (l.split("\t") for l in sh(f"bcftools query -R '{raw_bed}' -i 'GT=\"het\" && FILTER=\"PASS\" && TYPE=\"snp\"' "
                                        f"-f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' '{S['vcf_raw']}'", name="orig_q",
                                        capture=True).splitlines() if l)}
    conc = sum(1 for k in orig if k in recall) / max(len(orig), 1)
    lead_recalled = {describe_variant(v): (vkey(v) in recall) for v in TOPH["variants"]}
    # ---- empirical phaseability between the two leading alleles ------------------------------------------------------
    links = []
    if TOPH["model"] == "compound_het":
        a, b = sorted(TOPH["variants"], key=lambda v: int(v["pos"]))
        between = sorted({(k[0], k[1]) for k, gt in recall.items() if k[0] == a["chrom"] and int(a["pos"]) <= k[1] <= int(b["pos"])
                          and gt.replace("|", "/") in ("0/1", "1/0") and len(k[2]) == 1 and len(k[3]) == 1}
                         | {(a["chrom"], int(a["pos"])), (b["chrom"], int(b["pos"]))}, key=lambda x: x[1])
        links = fragment_links(bam, between)
    S["reads"] = dict(lanes=lanes_done, mapped_reads=n_reads, autosomal_control_depth=round(auto_med, 1), chrX_to_autosome_depth=round(x_ratio, 3),
                      leading_alleles=allele_rows, recalled_leading=lead_recalled, vcf_het_concordance=round(conc, 3), n_vcf_hets_in_bait=len(orig),
                      low_depth_runs_in_gene=del_runs, clip_clusters=clip_cl[:10], discordant_clusters=disc_cl[:10], phase_links=links,
                      demo=DEMO)
    broken = [l for l in links if l["linking_pairs"] == 0]
    S["reads"]["summary"] = (
        ("*Demo mode: toolchain test on 2 M public read pairs; coverage too low for interpretation.* " if DEMO else "")
        + f"{len(lanes_done)} lane(s) streamed through unique-minimizer recruitment; {n_reads:,} recruited reads re-aligned genome-wide. "
        f"Control-window depth {auto_med:.1f}× on autosomes; chrX/autosome depth {x_ratio:.2f} (expected 0.5 in a male, which validates "
        f"the depth method for one-copy changes). "
        + " ".join(f"{r['variant']}: {r['alt_reads']}/{r['ref_reads'] + r['alt_reads']} alt reads (fraction {r['alt_fraction']}; strands "
                   f"{r['alt_strands']}; lanes {r['alt_by_lane']}; original VCF AD {r['vcf_AD']})." for r in allele_rows)
        + f" Independent re-calling (bcftools) recovered {100 * conc:.0f} % of the {len(orig)} PASS heterozygous SNVs of the provided VCF in the "
          f"targets, and the leading alleles: {lead_recalled}. Depth along {TOPH['gene']}: "
        + (f"{len(del_runs)} run(s) of ≥ 3 bins below 0.65× ({del_runs})" if del_runs else "no run of low-depth bins (no exon-scale deletion)")
        + f"; clipped-read clusters: {len(clip_cl)}; discordant-pair clusters: {len(disc_cl)}. "
        + (f"Between the two alleles the independent re-call finds {len(links) - 1} heterozygous SNV site(s) (the provided VCF lists "
           f"{fmt((S.get('phase') or {}).get('n_het_sites_between'))}); {len(broken)} of the {len(links)} consecutive site pairs are linked by no read pair "
           f"(largest gap {max((l['gap_bp'] for l in links), default=0):,} bp)"
           + (": the reads cannot phase the alleles." if broken else ": every consecutive pair is linked, so read-level chaining is possible.")
           if links else ""))
    log(S["reads"]["summary"])
    if allele_rows:
        md("**Leading alleles in the re-aligned reads**\n\n" + df_md(pd.DataFrame(allele_rows)))
    fig, ax = plt.subplots(figsize=(9, 2.8))
    ax.scatter((tg.start + 250) / 1e6, tg.ratio, s=8, color=VIZ["s1"], linewidths=0, label=f"{TOPH['gene']} locus, 500-bp bins")
    ax.axhline(1.0, color=VIZ["axis"], linewidth=0.8)
    ax.axhline(0.5, color=VIZ["s2"], linewidth=1.2, label=f"one-copy level (chrX control measured {x_ratio:.2f})")
    for v in TOPH["variants"]:
        ax.axvline(int(v["pos"]) / 1e6, color=VIZ["ink2"], linewidth=0.8)
    ax.axvspan(gs / 1e6, ge / 1e6, color=VIZ["grid"], alpha=0.5, linewidth=0)
    ax.set_ylim(0, max(2.0, float(tg.ratio.quantile(0.99)) if len(tg) else 2.0))
    ax.set_xlabel(f"{GENE_POS[TOPH['gene']][0]} position (Mb); shaded = gene body; vertical lines = leading alleles")
    ax.set_ylabel("depth / autosomal control")
    ax.set_title(f"Read depth across the {TOPH['gene']} locus from recruited reads")
    ax.legend(loc="upper right", ncol=2)
    save_fig(fig, "locus_depth")

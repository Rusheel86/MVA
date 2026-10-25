# %% [markdown]
# ## 3 · Public reference resources
#
# Genome-wide steps run offline against downloaded public files; only a short list of top candidates is later sent
# to public annotation services.
#
# | Resource | Use | Source |
# |---|---|---|
# | GRCh38 no-alt analysis set | normalisation, consequence calling, read alignment (same primary assembly and coordinates as the VCF's reference) | NCBI (GCS mirror) |
# | MANE Select / Plus Clinical gene models and proteins | consequences on the clinically agreed transcript | NCBI RefSeq / Ensembl |
# | gnomAD v4.1 joint exomes+genomes, compact `echtvar` archive (AF, grpmax FAF95, homozygotes) | first-pass rarity, genome-wide and offline | gnomAD via Zenodo 22297079 (fallbacks: slivar gnotate v3, GATK af-only) |
# | gnomAD v4.1 joint sites (remote tabix, short list only) | exact AF, filtering AF, homozygotes | gnomAD |
# | ClinVar GRCh38 | prior clinical assertions, loss-of-function mechanism per gene | NCBI |
# | AlphaMissense hg38 | calibrated missense predictor (PP3/BP4) | Google DeepMind (CC BY-NC-SA) |
# | HPO ontology and disease annotations | phenotype engine | HPO / Monarch |

# %%
REF_GCS = "https://storage.googleapis.com/deepvariant/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
REF_NCBI = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/"
            "seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")
MANE_DIR = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
HPO_BASE = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/"
GNOMAD_V4_JOINT = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.{chrom}.vcf.bgz"
# first-pass population-frequency sources, tried in this order
AF_SOURCES = [
    dict(name="echtvar", label="gnomAD v4.1 joint exomes+genomes (echtvar archive, Zenodo 22297079): AF, grpmax FAF95, homozygotes",
         tool="https://github.com/brentp/echtvar/releases/download/v0.2.4/echtvar-v0.2.4-x86_64-unknown-linux-musl",
         data="https://zenodo.org/records/22297079/files/gnomad_4.1_whole_genome.zip?download=1", min_bytes=3_500_000_000),
    dict(name="slivar", label="gnomAD v3.0 genomes (slivar gnotate): max population AF, homozygotes",
         tool="https://github.com/brentp/slivar/releases/download/v0.3.4/slivar",
         data="https://slivar.s3.amazonaws.com/gnomad.hg38.genomes.v3.fix.zip", min_bytes=3_500_000_000),
    dict(name="afonly", label="gnomAD v2 lifted to GRCh38 (GATK af-only resource): global AF", tool=None,
         data="https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz", min_bytes=2_500_000_000),
]
R = WORK / "ref"
RES = {}
PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def _mane_files():
    listing = sh(f"curl -fsSL --retry 5 '{MANE_DIR}'", name="mane_listing", capture=True)
    out = {}
    for key, pat in (("gff", r"MANE\.GRCh38\.v[\d.]+\.ensembl_genomic\.gff\.gz"),
                     ("summary", r"MANE\.GRCh38\.v[\d.]+\.summary\.txt\.gz"),
                     ("protein", r"MANE\.GRCh38\.v[\d.]+\.ensembl_protein\.faa\.gz")):
        m = re.search(pat, listing)
        if not m:
            raise RuntimeError(f"MANE file not found in listing: {key}")
        out[key] = m.group(0)
    return out


with stage("download public resources", required=True):
    mane = _mane_files()
    S["mane_version"] = re.search(r"v[\d.]+\d", mane["gff"]).group(0)
    jobs = [
        ("ref", [REF_GCS], R / "ref.fa", 3_000_000_000),
        ("ref_fai", [REF_GCS + ".fai"], R / "ref.fa.fai", 5_000),
        ("mane_gff", [MANE_DIR + mane["gff"]], R / "mane.gff.gz", 5_000_000),
        ("mane_summary", [MANE_DIR + mane["summary"]], R / "mane.summary.txt.gz", 500_000),
        ("mane_protein", [MANE_DIR + mane["protein"]], R / "mane.protein.faa.gz", 5_000_000),
        ("clinvar", ["https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"], R / "clinvar.vcf.gz", 50_000_000),
        ("alphamissense", ["https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz",
                           "https://zenodo.org/records/10813168/files/AlphaMissense_hg38.tsv.gz?download=1"],
         R / "AlphaMissense_hg38.tsv.gz", 300_000_000),
        ("hp_obo", [HPO_BASE + "hp.obo", "http://purl.obolibrary.org/obo/hp.obo"], R / "hp.obo", 1_000_000),
        ("hpoa", [HPO_BASE + "phenotype.hpoa", "http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa"], R / "phenotype.hpoa", 1_000_000),
        ("genes_to_disease", [HPO_BASE + "genes_to_disease.txt", "http://purl.obolibrary.org/obo/hp/hpoa/genes_to_disease.txt"],
         R / "genes_to_disease.txt", 100_000),
        ("af_tool", [AF_SOURCES[0]["tool"]], WORK / "bin" / "echtvar", 3_000_000),
        ("af_data", [AF_SOURCES[0]["data"]], R / "af_echtvar.zip", AF_SOURCES[0]["min_bytes"]),
    ]

    def _job(j):
        key, urls, dest, min_bytes = j
        t = time.time()
        try:
            if min_bytes >= 2_000_000_000 and len(urls) == 1:        # multi-GB archives: parallel byte ranges
                # the gnomAD archive has fallbacks (AF_SOURCES), so a throttled mirror is abandoned after 15 min
                fetch_parallel(urls[0], dest, n=8, min_bytes=min_bytes, budget_s=900 if key == "af_data" else None)
            else:
                fetch_first(urls, dest, min_bytes=min_bytes)
            return key, str(dest), None, time.time() - t
        except Exception as e:      # noqa: BLE001
            return key, None, str(e)[:300], time.time() - t

    with ThreadPoolExecutor(max_workers=5) as ex:
        for key, path, err, dt in ex.map(_job, jobs):
            if err:
                log(f"  {key}: FAILED ({err})")
            else:
                RES[key] = path
                log(f"  {key}: {Path(path).stat().st_size / 1e6:,.0f} MB in {dt:.0f}s")
    for must in ("mane_gff", "mane_protein", "clinvar", "hp_obo", "hpoa", "genes_to_disease"):
        if must not in RES:
            raise RuntimeError(f"required resource missing: {must}")
    if "af_tool" in RES:
        Path(RES["af_tool"]).chmod(0o755)
    S["resource_dates"] = {"downloaded_utc": time.strftime("%Y-%m-%d %H:%M", time.gmtime())}

# %%
PRIMARY_SET = set(PRIMARY)


def load_mane(gff_gz, csq_gff_path):
    """MANE gene models as Python structures, plus a minimal Ensembl-style GFF3 that `bcftools csq` understands."""
    tx, gene_tx, gene_pos = {}, collections.defaultdict(list), {}
    with gzip.open(gff_gz, "rt") as fin, open(csq_gff_path, "w") as fout:
        fout.write("##gff-version 3\n")
        for line in fin:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if c[0] not in PRIMARY_SET:         # MANE models on fix / alt patch contigs: not in the reference or the VCF
                continue
            a = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            ftype = c[2]
            gid = a.get("gene_id", "").split(".")[0]
            tid = a.get("transcript_id", "").split(".")[0]
            if ftype == "gene":
                gene_pos[a["gene_name"]] = (c[0], int(c[3]), int(c[4]), c[6], gid)
                fout.write("\t".join(c[:8] + [f"ID=gene:{gid};biotype=protein_coding;Name={a['gene_name']}"]) + "\n")
            elif ftype == "transcript":
                tx[tid] = dict(gene=a["gene_name"], chrom=c[0], strand=c[6], start=int(c[3]), end=int(c[4]), cds=[], exons=[],
                               protein=a.get("protein_id", ""), refseq=a.get("Dbxref", "").replace("RefSeq:", ""),
                               tag=a.get("tag", ""), version=a.get("transcript_id", ""))
                gene_tx[a["gene_name"]].append(tid)
                fout.write("\t".join(c[:8] + [f"ID=transcript:{tid};Parent=gene:{gid};biotype=protein_coding"]) + "\n")
            elif ftype in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR") and tid in tx:
                if ftype == "CDS":
                    tx[tid]["cds"].append((int(c[3]), int(c[4])))
                elif ftype == "exon":
                    tx[tid]["exons"].append((int(c[3]), int(c[4])))
                fout.write("\t".join(c[:8] + [f"Parent=transcript:{tid}"]) + "\n")
    for t in tx.values():       # transcription order
        rev = t["strand"] == "-"
        t["cds"].sort(reverse=rev)
        t["exons"].sort(reverse=rev)
        t["cds_len"] = sum(e - s + 1 for s, e in t["cds"])
    # prefer the MANE Select transcript when a gene also has a MANE Plus Clinical one
    for g in gene_tx:
        gene_tx[g].sort(key=lambda t: 0 if "MANE_Select" in tx[t]["tag"] else 1)
    return tx, gene_tx, gene_pos


def load_proteins(faa_gz):
    prot, name = {}, None
    with gzip.open(faa_gz, "rt") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                prot[name] = []
            elif name:
                prot[name].append(line.strip())
    return {k: "".join(v) for k, v in prot.items()}


def load_clinvar_plp(tsv):
    """P/LP (non-conflicting) ClinVar records; loss of function counts as an established disease mechanism for a
    gene when ClinVar holds >= 2 P/LP truncating alleles in it."""
    cv = pd.read_csv(tsv, sep="\t", header=None, dtype=str, keep_default_na=False, quoting=3,
                     names=["chrom", "pos", "ref", "alt", "clnsig", "clnrevstat", "geneinfo", "mc", "clndn"])
    cv = cv[~cv.clnsig.str.startswith("Conflicting")].copy()
    cv["pos"] = cv["pos"].astype(int)
    cv["gene"] = cv.geneinfo.str.split(":").str[0]
    cv["is_lof"] = cv.mc.str.contains("nonsense|frameshift|splice_donor|splice_acceptor", regex=True)
    cv["is_missense"] = cv.mc.str.contains("missense")
    n_lof = cv[cv.is_lof].groupby("gene").size()
    return cv, set(n_lof[n_lof >= 2].index)


with stage("prepare reference, gene models, ClinVar", required=True):
    # ---- reference FASTA (GRCh38 no-alt analysis set; same primary sequence and coordinates as the VCF's reference) ----
    REF = R / "ref.fa"
    if "ref" not in RES or "ref_fai" not in RES:        # fall back to the NCBI copy (plain gzip)
        gz = fetch(REF_NCBI + ".gz", R / "ref.fa.gz", min_bytes=500_000_000)
        sh(f"(pigz -dc '{gz}' || gzip -dc '{gz}') > '{REF}' && samtools faidx '{REF}'", name="ref_unpack")
        Path(gz).unlink(missing_ok=True)
    RES["ref"] = str(REF)
    CHROM_LEN = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(str(REF) + ".fai")}
    assert all(c in CHROM_LEN for c in PRIMARY), "reference is missing primary chromosomes"

    # ---- MANE ----------------------------------------------------------------------------------------------------------
    RES["csq_gff"] = str(R / "mane.csq.gff3")
    TX, GENE_TX, GENE_POS = load_mane(RES["mane_gff"], RES["csq_gff"])
    PROT = load_proteins(RES["mane_protein"])
    log(f"MANE {S['mane_version']}: {len(GENE_POS):,} genes, {len(TX):,} transcripts, {len(PROT):,} proteins")

    # ---- ClinVar: chr-prefixed, primary chromosomes only, with contig lines; plus a P/LP table ----------------------
    cv = R / "clinvar.chr.vcf.gz"
    if not Path(str(cv) + ".tbi").exists():
        keep_chr = {c[3:] for c in PRIMARY}
        plain = R / "clinvar.chr.vcf"
        with gzip.open(RES["clinvar"], "rt") as fin, open(plain, "w") as fout:
            for line in fin:
                if line.startswith("##"):
                    fout.write(line)
                elif line.startswith("#CHROM"):
                    fout.write("".join(f"##contig=<ID={c},length={CHROM_LEN[c]}>\n" for c in PRIMARY))
                    fout.write(line)
                else:
                    c, rest = line.split("\t", 1)
                    if c in keep_chr:
                        fout.write("chr" + c + "\t" + rest)
        sh(f"bgzip -f -@ {NCPU} '{plain}' && tabix -f -p vcf '{cv}'", name="clinvar_chr")        # writes clinvar.chr.vcf.gz
    RES["clinvar_chr"] = str(cv)
    S["clinvar_date"] = (sh(f"(gzip -dc '{RES['clinvar']}' | head -40 | grep -m1 '^##fileDate') || true", name="clinvar_date",
                            capture=True).strip().replace("##fileDate=", "") or "n/a")
    plp_tsv = R / "clinvar.plp.tsv"
    sh(f"bcftools query -i 'INFO/CLNSIG~\"athogenic\"' -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/GENEINFO\t%INFO/MC\t%INFO/CLNDN\n' '{cv}' > '{plp_tsv}'",
       name="clinvar_plp")
    CLINVAR_PLP, LOF_GENES = load_clinvar_plp(plp_tsv)
    log(f"ClinVar {S['clinvar_date']}: {len(CLINVAR_PLP):,} P/LP records; LoF mechanism established in {len(LOF_GENES):,} genes")

    # ---- AlphaMissense: the public file is BGZF, so it only needs a tabix index --------------------------------------
    if "alphamissense" in RES and not Path(RES["alphamissense"] + ".tbi").exists():
        try:
            sh(f"tabix -f -s 1 -b 2 -e 2 -c '#' '{RES['alphamissense']}'", name="am_index")
        except Exception as e:      # noqa: BLE001
            log(f"AlphaMissense index failed ({str(e)[:120]}); missense scores will come from the short-list annotation only")
            RES.pop("alphamissense", None)

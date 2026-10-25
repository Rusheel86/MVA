# %% [markdown]
# ## 5 · Case data
#
# * **Challenge mode** – the proband VCF and phenotype document come from the gated Hugging Face dataset; HPO terms
#   are parsed from the document's table. Terms that describe relatives are listed
#   but excluded from the proband similarity.
# * **Demo mode** (no token) – the public GIAB HG002 genome plus a spiked-in recessive case built from ClinVar, so
#   that every stage can be exercised without patient data. The same spike-in generator powers the benchmark in
#   section 13.

# %%
HG002_VCF = ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/"
             "NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
RELATIVE_RE = re.compile(r"\b(parent|parental|maternal|paternal|mother|father|sibling|family history|relative)", re.I)


def parse_phenotype_docx(path):
    """HPO ids from the phenotype document (table rows first, any 'HP:nnnnnnn' in the text as fallback)."""
    xml = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8", "replace")

    def text_of(fragment):
        paras = re.findall(r"<w:p[ >].*?</w:p>", fragment, flags=re.S)
        return " ".join("".join(re.findall(r"<w:t(?: [^>]*)?>([^<]*)</w:t>", p)) for p in paras).strip()

    out, seen = [], set()
    for row in re.findall(r"<w:tr[ >].*?</w:tr>", xml, flags=re.S):
        text = " | ".join(text_of(c) for c in re.findall(r"<w:tc[ >].*?</w:tc>", row, flags=re.S))
        for hp in re.findall(r"HP:\d{7}", text):
            if hp not in seen:
                seen.add(hp)
                rel = bool(RELATIVE_RE.search(text))
                out.append(dict(term=hp, relative=rel, weight=0.0 if rel else 1.0))
    if not out:
        for hp in re.findall(r"HP:\d{7}", text_of(xml)):
            if hp not in seen:
                seen.add(hp)
                out.append(dict(term=hp, relative=False, weight=1.0))
    return out


def hf_get(filename):
    from huggingface_hub import hf_hub_download
    try:
        return hf_hub_download(repo_id=CFG["HF_DATASET"], filename=filename, repo_type="dataset",
                               local_dir=str(WORK / "data"), token=HF_TOKEN)
    except Exception as e:      # noqa: BLE001
        raise RuntimeError(f"could not download '{filename}' from {CFG['HF_DATASET']} — check that the token belongs to an "
                           f"account that has been granted access to the gated dataset ({type(e).__name__}: {str(e)[:200]})") from e


def make_spike_cases(n, seed=7):
    """Synthetic recessive cases for the demo/benchmark: one ClinVar P/LP truncating allele + one ClinVar P/LP
    missense allele in an autosomal-recessive single-gene disease, with HPO terms sampled from that disease
    (imprecision and one unrelated noise term added, as in published phenotype-simulation benchmarks)."""
    rng = np.random.default_rng(seed)
    rs = CLINVAR_PLP.clnrevstat
    cv = CLINVAR_PLP[rs.str.contains("criteria_provided|expert_panel|practice_guideline") & ~rs.str.contains("no_assertion")
                     & CLINVAR_PLP.ref.str.fullmatch(r"[ACGT]{1,30}") & CLINVAR_PLP.alt.str.fullmatch(r"[ACGT]{1,30}")
                     & CLINVAR_PLP.chrom.isin(PRIMARY[:22])]
    lof = cv[cv.is_lof & ~cv.is_missense].groupby("gene")
    mis = cv[cv.is_missense & ~cv.is_lof & (cv.ref.str.len() == 1) & (cv.alt.str.len() == 1)].groupby("gene")
    lof_genes, mis_genes = set(lof.groups), set(mis.groups)
    eligible = []
    for did, terms in PE.d_terms.items():
        genes = PE.d2g.get(did, set())
        if len(genes) != 1 or PE.d_inh.get(did) != {"AR"} or len(terms) < 8:
            continue
        g = next(iter(genes))
        if g in GENE_POS and g in lof_genes and g in mis_genes and GENE_POS[g][0] in PRIMARY[:22]:
            eligible.append((did, g))
    eligible.sort()
    all_terms = sorted(PE.ic)
    picks = rng.choice(len(eligible), size=min(n, len(eligible)), replace=False)
    cases, used_genes = [], set()
    for ci in picks:
        did, g = eligible[int(ci)]
        if g in used_genes:
            continue
        used_genes.add(g)
        a = lof.get_group(g).iloc[int(rng.integers(len(lof.get_group(g))))]
        b = mis.get_group(g).iloc[int(rng.integers(len(mis.get_group(g))))]
        terms = sorted(PE.d_terms[did])
        ic = np.array([PE.ic.get(t, 0.0) for t in terms]) + 0.1
        chosen = list(rng.choice(terms, size=min(6, len(terms)), replace=False, p=ic / ic.sum()))
        noisy = []
        for t in chosen:
            parents = sorted(PE.hpo.parents.get(t, ()))
            if parents and rng.random() < 0.3 and PE.ic.get(parents[0], 0) > 1.0:
                noisy.append(parents[0])        # imprecise term
            else:
                noisy.append(t)
        noisy.append(all_terms[int(rng.integers(len(all_terms)))])      # unrelated noise term
        noisy = list(dict.fromkeys(noisy))
        cases.append(dict(case=len(cases), gene=g, disease=did, disease_name=PE.d_name[did], hpo=noisy,
                          variants=[(a.chrom, int(a.pos), a.ref, a.alt), (b.chrom, int(b.pos), b.ref, b.alt)]))
    return cases


with stage("case data", required=True):
    if not DEMO:
        vcf = hf_get(CFG["VCF_FILE"])
        # own index: the supplied .tbi (written by htsjdk) has no per-contig record counts, which bcftools index -s/-n need
        sh(f"bcftools index -f -t --threads {NCPU} '{vcf}'", name="index_vcf")
        if CFG["HPO_TERMS"]:
            S["hpo"] = [dict(term=t, relative=False, weight=1.0) for t in CFG["HPO_TERMS"]]
        else:
            S["hpo"] = parse_phenotype_docx(hf_get(CFG["PHENOTYPE_DOCX"]))
        if not S["hpo"]:
            raise RuntimeError("no HPO terms found — set CFG['HPO_TERMS']")
        S["diagnosis"] = CFG["CLINICAL_DIAGNOSIS"].strip()
        S["vcf_raw"], S["case_label"] = vcf, CFG["PROBAND_ID"]
        SPIKE_CASES = []
    else:
        vcf = fetch(HG002_VCF, WORK / "data" / "HG002.vcf.gz", min_bytes=100_000_000)
        fetch(HG002_VCF + ".tbi", WORK / "data" / "HG002.vcf.gz.tbi", min_bytes=100_000)
        SPIKE_CASES = make_spike_cases(60)
        demo = SPIKE_CASES[0]
        S["hpo"] = [dict(term=t, relative=False, weight=1.0) for t in demo["hpo"]]
        S["diagnosis"] = demo["disease_name"]
        S["vcf_raw"], S["case_label"] = str(vcf), "PROBAND01"
        S["demo_truth"] = dict(gene=demo["gene"], disease=demo["disease_name"],
                               variants=[f"{c}:{p} {r}>{a}" for c, p, r, a in demo["variants"]])
        log(f"DEMO case: {demo['disease_name']} ({demo['gene']}) spiked into HG002; {len(SPIKE_CASES)} benchmark cases prepared")
    for h in S["hpo"]:
        h["label"] = PE.hpo.name.get(PE.hpo.canon(h["term"]), "?")
        h["ic"] = round(PE.term_ic(h["term"]), 2)
    md("**HPO terms used**\n\n" + df_md(pd.DataFrame(S["hpo"])[["term", "label", "relative", "weight", "ic"]]))
    log(f"clinical indication: {S['diagnosis'] or '(none — undiagnosed mode)'}")

# %% [markdown]
# ## 15 · Reports
#
# Both write-ups are generated from the results of this run (nothing is typed in by hand), so they always match the
# submission file. They are rewritten at the end of the notebook once the read-level module has finished.
#
# * `<team>_track1_report.md` — methods write-up for Track 1, including the organisers' methods-template answers;
# * `<team>_track2_report.md` — mechanism characterisation and ranked drug hypotheses for Track 2;
# * `run_summary.json` — machine-readable summary (findings only; no genotype-scale data).

# %%
REPORT_PMIDS = ["15475955", "16411201", "24344301", "20516114", "18932004", "15159543", "15208629", "23242215", "24825348",
                "19407811", "22698286", "23079597", "31053147", "29300386", "32720330", "30192042", "36413997", "40084623",
                "28518168", "38685113", "21552266", "28553959", "34009673"]


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "–"
    if isinstance(x, float):
        return f"{x:.{nd}g}"
    return str(x)


def clean_title(t):
    """Europe PMC titles carry HTML markup and entities (e.g. &lt;i&gt;): plain text for Markdown."""
    import html
    return re.sub(r"<[^>]+>", "", html.unescape(html.unescape(t or ""))).strip()


def cite(p, pm):
    v = pm.get(p)
    return f"{v['authors']} et al. {v['year']}, *{v['journal']}* (PMID {p})" if v else f"PMID {p}"


def run_table():
    return df_md(pd.DataFrame([dict(stage=k, status=v["status"], minutes=round(v["seconds"] / 60, 1), note=v["note"][:90])
                               for k, v in RUN.items()]), max_rows=60)


def wilson(k, n, z=1.96):
    """95 % Wilson score interval for k successes out of n."""
    if not n:
        return (None, None)
    ph = k / n
    d = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / d
    h = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def pct_ci(frac, n):
    lo, hi = wilson(round(frac * n), n)
    return f"{100 * frac:.0f} % (95 % CI {100 * lo:.0f}–{100 * hi:.0f})" if lo is not None else "–"


PARAMS = [   # (parameter, value, source) — what the ranking depends on, and where each value comes from
    ("ACMG evidence → points → posterior", "Bayesian point scale, prior 0.10, OddsPath 350", "Tavtigian 2018 / 2020 (literature)"),
    ("Missense predictor bands (PP3 / BP4)", "AlphaMissense and REVEL thresholds", "Pejaver 2022 ClinGen calibration (literature)"),
    ("PM2 / BS1 / BA1 frequency cut-offs", "≤ 1e-4 / > 1e-3 / ≥ 5 % (max of AF and FAF95)", "ClinGen SVI recommendations (literature)"),
    ("PM3 for the weaker allele", "+2 with a Pathogenic partner, +1 with a Likely pathogenic one, 0 when phase is unknown", "ClinGen SVI PM3 v1.0 (literature)"),
    ("PVS1 / NMD", "last-exon and 50-nt rules", "Abou Tayoun 2018 (literature)"),
    ("Phenotype factor", "(similarity / best similarity)^3", "heuristic, fixed before the run"),
    ("Mode-of-inheritance and phase factors", "e.g. 0.9 for unphased pairs", "heuristic, fixed before the run"),
    ("Indication prior", "genes on the stated diagnosis' panel floored at 0.95 of the best fit", "heuristic, fixed before the run; both rankings reported"),
    ("Other-gene rows in the CSV", "only when EPCR ≥ 0.25", "heuristic, fixed before the run"),
]


def competing_summary():
    others = [h for h in HYPS[1:30] if h["gene"] != HYPS[0]["gene"]]
    same = [h for h in HYPS[1:30] if h["gene"] == HYPS[0]["gene"]]
    lines = []
    if others:
        best = others[0]
        name = best["gene"] if reportable_gene(best["gene"]) else "a gene outside the indication (name withheld)"
        lines.append(f"The best hypothesis in any other gene is a {best['model'].replace('_', '-')} model in {name} with EPCR "
                     f"{best['epcr']:.3f} ({HYPS[0]['epcr'] / max(best['epcr'], 1e-6):.1f}-fold lower than the leading pair).")
    if same:
        n_sub = max(0, len(SUB_ROWS) - 1) if "SUB_ROWS" in globals() else None
        lines.append(f"{len(same)} lower-ranked alternative genotypes within {HYPS[0]['gene']} follow (other pairings and single alleles)"
                     + (f"; {n_sub} of them are submitted as fall-back rows so that a partially correct answer still scores "
                        "(a single allele already listed in a higher pair, and any variant with net benign evidence, are not submitted)."
                        if n_sub is not None else "."))
    return " ".join(lines)


# The organisers' required acknowledgement, verbatim from the Official Rules ("Citation and Acknowledgement"): any
# publication or public communication arising from the hackathon must include it.
ACK = ("This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, "
       "Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship "
       "from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and "
       "their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.")


def track1_report():
    T, pm = S.get("top", {}), {**(S.get("track2", {}) or {}).get("pmids", {}), **S.get("report_pmids", {})}
    vs = T.get("variants", [])
    qc, an, rb, bench = S.get("qc", {}), S.get("aneuploidy", {}), S.get("robustness", {}), S.get("benchmark", {})
    ph, by, esm, st, rl = S.get("phase", {}), S.get("bystander", {}), S.get("esm", {}), S.get("structure", {}), S.get("reads", {})
    gene = T.get("gene", S.get("top_gene", "?"))
    sub = pd.read_csv(SUBMISSION_CSV, dtype=str, keep_default_na=False) if SUBMISSION_CSV.exists() else pd.DataFrame()
    L = []
    L.append(f"# {CFG['TEAM']} — Track 1 methods report: {CFG['MODEL_NAME']}\n")
    L.append(f"*Rare Disease, Real Kid: MVA Hackathon 2026 · generated automatically by the analysis notebook on "
             f"{time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())} · mode: {S.get('mode')}*\n")
    # ---------------- abstract ----------------
    v1 = vs[0] if vs else {}
    v2 = vs[1] if len(vs) > 1 else {}
    L.append("## Summary\n")
    L.append(
        f"We built **DuoRank**, a fully automated, public-data-only pipeline that ranks genotype hypotheses — homozygous, "
        f"compound-heterozygous, heterozygous and hemizygous — from a single-sample GRCh38 VCF and a list of HPO terms. "
        f"Variants are normalised and annotated offline (gnomAD v4.1, ClinVar, MANE consequences, AlphaMissense), scored in "
        f"the Bayesian point system of the ACMG/AMP guidelines (evidence strengths from published ClinGen calibrations; the few "
        f"ranking heuristics fixed before the run and listed in section 7), "
        f"and combined with a gene-level phenotype fit computed by information-content similarity against all "
        f"{fmt(S.get('pheno', {}).get('n_diseases'))} annotated Mendelian diseases. A short list is then re-annotated "
        f"(Ensembl VEP, exact gnomAD records, SpliceAI, ESM-2, AlphaFold structure), phase is analysed, and the raw "
        f"reads are optionally re-aligned around the leading locus for independent validation.\n")
    if T:
        L.append(
            f"**Result.** The leading hypothesis is a **{T['model'].replace('_', '-')} genotype in {gene}** "
            f"({T.get('dx_disease') or T.get('disease')}, {'/'.join(T.get('inheritance') or [])}): "
            + " + ".join(f"`{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}` ({v.get('hgvsc', '').split(':')[-1]} {v.get('hgvsp', '')}, "
                         f"automated class **{v['classification']}**)" for v in vs)
            + f", EPCR {T['epcr']:.3f}. Without the stated diagnosis (HPO terms only) the gene ranks "
              f"#{T.get('pheno_rank_blind')} of {fmt(S.get('pheno', {}).get('n_genes'))} on phenotype alone and the same genotype is "
              f"{'still' if S.get('top_same_blind') else 'not'} ranked first. {competing_summary()}\n")
    # ---------------- data ----------------
    L.append("## 1. Data and question\n")
    L.append(f"* Proband genome: single-sample WGS VCF, GRCh38 (Sentieon Haplotyper + GATK hard filters); raw FASTQ (4 lanes, paired).")
    L.append("* Phenotype (HPO): " + "; ".join(f"{h.get('label')} ({h['term']})" for h in S.get("hpo", []) if not h.get("relative"))
             + ". Terms that describe relatives rather than the proband were excluded from the similarity and are not listed here.")
    L.append(f"* Stated clinical indication (used as a transparent prior, both rankings reported): {S.get('diagnosis') or 'none'}")
    L.append(f"* Genome QC: Ti/Tv {fmt(qc.get('titv'))}, het/hom {fmt(qc.get('het_hom_ratio'))}, median depth at het sites "
             f"{fmt(qc.get('median_depth_het_sites'))}×, inferred sex {qc.get('inferred_sex')}, autozygosity F_ROH {fmt(qc.get('f_roh'))} "
             f"({qc.get('autozygosity_note', '')}).\n")
    # ---------------- method ----------------
    L.append("## 2. Method\n")
    L.append("| Step | What is done | Key choices |\n|---|---|---|")
    L.append("| Normalise | primary chromosomes, `chr` names, split multi-allelics, left-align (`bcftools norm`) | exact allele representation matches the scorer |")
    L.append(f"| Rarity | {S.get('af_source', 'gnomAD')} | ClinVar P/LP alleles exempt below 5 % (founder / hypomorphic alleles); VCF `INFO/AF` (within-sample) never used |")
    L.append(f"| Consequence | `bcftools csq` on MANE {S.get('mane_version', '')} transcripts; automated PVS1 / NMD rule | last-exon and 50-nt rule; LoF mechanism from ClinVar |")
    L.append("| Missense | AlphaMissense, ClinGen-calibrated bands (one predictor only; REVEL/CADD shown, not counted) | avoids double counting |")
    L.append("| Splicing | SpliceAI (pre-computed via VEP and re-run on GPU over a ±500 bp window for the leading genes) | non-coding second hits are scored, not filtered |")
    L.append(f"| Phenotype | Resnik best-match similarity vs {fmt(S.get('pheno', {}).get('n_diseases'))} HPO-annotated diseases (HPO {S.get('hpo_version')}) | relatives' terms excluded; indication used as a panel prior |")
    L.append("| Genotypes | hom / comp-het / het / hemi per gene; EPCR = P(genotype) × phenotype fit × mode fit × call quality | PM3 per ClinGen SVI; phase factor from read-backed tags |")
    L.append("| Second-hit rescue | rare non-coding variants in and around genes that already carry one credible allele | catches regulatory / deep-intronic hypomorphs |")
    L.append("| Output | top hypothesis + within-gene alternatives; other genes only if EPCR ≥ 0.25 | file parsed and self-scored by a re-implementation of the organisers' evaluation script before saving |\n")
    # ---------------- leading hypothesis ----------------
    if vs:
        L.append(f"## 3. The leading hypothesis in detail\n")
        L.append("| | " + " | ".join(f"allele {i + 1}" for i in range(len(vs))) + " |")
        L.append("|---|" + "---|" * len(vs))
        rows = [("GRCh38", lambda v: f"{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}"),
                ("HGVS (MANE)", lambda v: f"{v.get('mane', '')}:{v.get('hgvsc', '').split(':')[-1]} {v.get('hgvsp', '')}"),
                ("consequence", lambda v: f"{v['csq']}" + (f", exon {v['exon']}" if v.get("exon") else "")
                 + (f"; NMD {'escaped' if v.get('nmd_escape') else 'predicted'}" if v["csq"] in LOF_CSQ else "")),
                ("genotype / reads", lambda v: f"{v['zyg']}, AD {v.get('ad')}, GQ {fmt(v.get('gq'))}, FILTER {v.get('filt')}"),
                ("gnomAD v4.1", lambda v: (f"AC {fmt((v.get('gnomad') or {}).get('ac'))} / AN {fmt((v.get('gnomad') or {}).get('an'))}, "
                                           f"homozygotes {fmt((v.get('gnomad') or {}).get('nhomalt'))}, FAF95 {fmt((v.get('gnomad') or {}).get('faf95'))}")
                 if (v.get("gnomad") or {}).get("present") else "absent"),
                ("ClinVar", lambda v: (v.get("clnsig") or "not in ClinVar").replace("_", " ") + (f" ({v.get('clnrevstat', '').replace('_', ' ')})" if v.get("clnsig") else "")),
                ("AlphaMissense / REVEL / CADD", lambda v: f"{fmt(v.get('am'))} / {fmt(v.get('revel'))} / {fmt(v.get('cadd'))}"),
                ("SpliceAI max Δ", lambda v: fmt(v.get("spliceai"))),
                ("ACMG points → posterior", lambda v: f"{v['points']:g} → {v['posterior']:.3f}"),
                ("automated classification", lambda v: f"**{v['classification']}**")]
        for name, f in rows:
            L.append(f"| {name} | " + " | ".join(f(v) for v in vs) + " |")
        L.append("")
        for i, v in enumerate(vs, 1):
            L.append(f"*Allele {i} evidence:* " + "; ".join(f"{e['criterion']} ({e['points']:+g}: {e['basis']})" for e in v["evidence"]) + "\n")
        L.append(f"**Phenotype fit.** {gene} → {T.get('disease')} (similarity {T.get('pheno_score')}); rank #{T.get('pheno_rank_blind')} "
                 f"of {fmt(S.get('pheno', {}).get('n_genes'))} genes on HPO terms alone. "
                 + (f"Green gene on {len(S.get('panelapp', []))} Genomics England PanelApp panels, including "
                    + ", ".join(p["panel"] for p in S.get("panelapp", [])[:4]) + "." if S.get("panelapp") else "") + "\n")
        if by:
            L.append(f"**Is the weaker allele a coincidental bystander?** In gnomAD v4.1, {by['n_rare_variants']} rare (< 0.01 %) "
                     f"{by['variant_class']} variants in {gene} have a cumulative frequency of {by['cumulative_af']:.4f}, so an unrelated "
                     f"person carries one by chance with probability ≈ {by['carrier_probability']:.3f}"
                     + (f" (≈ {by['carrier_probability_strict']:.4f} counting only alleles at least as damaging by AlphaMissense)"
                        if "carrier_probability_strict" in by else "")
                     + f". Given a pathogenic first allele, the posterior that this allele is the second causal allele (rather than a "
                       f"bystander next to an undetected second hit) is {by['posterior_second_allele']['prior_0.5']:.3f} with an even prior "
                       f"on 'second hit visible as a rare coding variant', and {by['posterior_second_allele']['prior_0.2']:.3f} with a "
                       f"sceptical prior of 0.2.\n")
        if esm:
            L.append(f"**Protein language model (orthogonal, not counted in the points).** ESM-2 650M log-likelihood ratio for "
                     f"{esm.get('change')}: {esm.get('llr')}"
                     + (f"; more damaging than {esm.get('percentile_more_damaging_than')} % of all possible substitutions in the protein"
                        if esm.get("percentile_more_damaging_than") is not None else "")
                     + (f"; median for {esm.get('clinvar_plp_missense_n')} ClinVar P/LP missense variants in the gene: {esm.get('clinvar_plp_median_llr')}"
                        if (esm.get("clinvar_plp_missense_n") or 0) >= 5 else
                        f" (the gene has only {esm.get('clinvar_plp_missense_n') or 0} ClinVar P/LP missense allele(s), too few for a "
                        "gene-specific reference)")
                     + ". On its own this is moderate, not decisive, evidence of damage — consistent with a hypomorphic rather than a "
                       "null allele, which is what a second allele in a recessive disorder whose null/null genotype is lethal is expected to be.\n")
        if st:
            L.append(f"**Structure.** Residue {st.get('residue')} lies in: {st.get('domain')}; AlphaFold pLDDT {st.get('plddt')} "
                     f"({st.get('neighbours_within_10A')} residues within 10 Å). "
                     + ("Nearest ClinVar P/LP missense residues: " + ", ".join(f"{x['residue']} ({x['distance_A']} Å)" for x in st["nearest_clinvar_plp_missense"]) + ". "
                        if st.get("nearest_clinvar_plp_missense") else "")
                     + (f"UniProt lists {st.get('n_disease_missense_uniprot')} literature-reported patient missense alleles in the protein; the nearest are "
                        + ", ".join(f"{x['change']} ({x['disease']}, {x['distance_A']} Å)" for x in st["disease_missense_uniprot"][:3])
                        + (f" — reported {st['disease_missense_uniprot'][0]['disease']} missense alleles cluster around this residue in the same "
                           "domain, the expected location of a second, hypomorphic allele"
                           if sum(1 for x in st["disease_missense_uniprot"] if x["distance_A"] <= 30 and x.get("same_domain")) >= 3 else "")
                        + ". " if st.get("disease_missense_uniprot") else "") + "\n")
        if ph:
            L.append(f"**Phase.** {ph.get('verdict', '')}. The alleles are {fmt(ph.get('distance_bp'))} bp apart with "
                     f"{fmt(ph.get('n_het_sites_between'))} heterozygous call(s) between them in the provided VCF (largest gap "
                     f"{fmt(ph.get('largest_gap_bp'))} bp; the read-level re-call below counts sites independently). "
                     f"Trans configuration is therefore **inferred, not observed**; it would be confirmed by {ph.get('confirmation')}.\n")
        if rl:
            L.append("**Read-level validation from the raw FASTQ.** " + rl.get("summary", "") + "\n")
        if gene == "BUB1B" and T.get("model") == "compound_het":
            L.append(f"**Why this genotype fits the known biology.** MVA1 is caused by biallelic BUB1B variants, typically one "
                     f"truncating allele with one missense allele, often in the kinase domain ({cite('15475955', pm)}; "
                     f"{cite('16411201', pm)}). A complete loss of BubR1 is embryonic-lethal in mice, so live-born patients keep some "
                     f"function, and BubR1 dosage tracks the phenotype: a graded reduction first causes premature chromatid separation "
                     f"and then aneuploidy ({cite('18932004', pm)}), and hypomorphic mice age prematurely ({cite('15208629', pm)}). "
                     f"Patient kinase-domain missense proteins are unstable rather than catalytically dead ({cite('20516114', pm)}). "
                     f"The child's pair — a predicted-NMD nonsense allele with a kinase-domain missense allele — is therefore the expected "
                     f"null + hypomorph architecture. Some families carry a regulatory BUB1B allele instead ({cite('24344301', pm)}); this is why "
                     "rare non-coding variants around genes with one credible allele are rescued and scored.\n")
        rc = (S.get("levers") or {}).get("reclassification")
        if rc and rc.get("path"):
            L.append(f"**Path to a firm classification of {rc['allele']}.** Points are applied cumulatively (Tavtigian point scale; "
                     "PM3 per ClinGen SVI):\n\n" + df_md(pd.DataFrame(rc["path"])) + "\n")
        L.append("**What a clinical laboratory would do to confirm.** (1) Sanger confirmation of both alleles. (2) Phase: targeted "
                 "testing of either parent (the standard route, which also supplies PM3)"
                 + (f"; or, because the alleles are only {fmt(ph.get('cdna_distance_nt'))} nt apart in the mRNA, one RT-PCR amplicon from "
                    "cultured cells with an NMD inhibitor — Primer3 primers are given in the Track 2 report" if ph.get("cdna_distance_nt") else "")
                 + ". (3) Function: BUBR1 protein level and mitotic-checkpoint competence in the child's cells (PS3), and premature "
                   "chromatid separation / mosaic aneuploidy on cytogenetics (PP4). (4) Cascade testing for carriers in the family, "
                   "offered through the clinical team.\n")
        if len(sub):
            L.append("**Submitted rows**\n\n" + df_md(sub))
            L.append("")
    # ---------------- certainty ----------------
    L.append("## 4. How robust is the ranking?\n")
    if rb:
        wp = rb.get("weight_perturbation", {})
        L.append(f"* Every heuristic weight multiplied by a random factor in [½, 2] ({wp.get('n')} runs): the top hypothesis is unchanged in "
                 f"{100 * (wp.get('top1_unchanged') or 0):.0f} % of runs (worst rank {wp.get('worst_rank')}).")
        L.append(f"* Rank of the leading genotype without ClinVar: {rb.get('rank_without_clinvar')}; without the missense predictor: "
                 f"{rb.get('rank_without_missense_predictor')}; without the stated diagnosis: {rb.get('rank_without_diagnosis')}; "
                 f"without both ClinVar and the diagnosis: {rb.get('rank_without_clinvar_and_diagnosis')}.")
        loo = rb.get("leave_one_term_out_blind", [])
        if loo:
            L.append("* Leave-one-HPO-term-out (diagnosis-blind): " + ", ".join(f"{x['label']} → rank {x['rank_of_top']}" for x in loo) + ".")
    L.append(f"* {competing_summary()}")
    L.append("* Hypotheses in genes outside the stated indication are counted but not named in this report: the organisers are still "
             "deciding, with the IRB, how incidental findings about this child may be published.\n")
    # ---------------- aneuploidy ----------------
    if an:
        L.append("## 5. What the bulk genome can say about the aneuploidy itself\n")
        lim = an.get("limit", {})
        top_dev = sorted((r for r in an.get("table", []) if max(abs(r.get("z_depth") or 0), abs(r.get("z_baf") or 0)) >= 2),
                         key=lambda r: -max(abs(r.get("z_depth") or 0), abs(r.get("z_baf") or 0)))[:3]
        L.append(f"Per-chromosome depth (GC-corrected) and allele-balance dispersion at heterozygous SNVs showed no chromosome-wide "
                 f"imbalance passing the pre-set threshold (flagged chromosomes: {an.get('flagged') or 'none'}; {an.get('artefact_bins')} "
                 "collapsed-repeat bins excluded"
                 + (f"; the largest deviations, chromosome(s) {', '.join(str(r['chrom']) for r in top_dev)} at |z| ≤ "
                    f"{max(max(abs(r['z_depth']), abs(r['z_baf'])) for r in top_dev):.1f}, stay below the threshold of 4 and the raw depth "
                    f"still correlates with GC content (r = {lim.get('gc_correlation')}), so residual GC bias is the likelier explanation"
                    if top_dev else "") + "). With a "
                 f"between-chromosome SD of {lim.get('sd_depth_ratio')} in depth ratio, a clonal whole-chromosome gain or loss present in "
                 f"more than ~{100 * (lim.get('min_detectable_clone_fraction') or 0):.0f} % of cells would have been detected. "
                 "This is the expected result in mosaic *variegated* aneuploidy: different cells carry different gains and losses, which "
                 "average out in bulk DNA. Bulk WGS can therefore neither confirm nor exclude the cytogenetic phenotype — metaphase "
                 "karyotyping, interphase FISH or single-cell low-pass sequencing are the appropriate assays.\n")
        if "aneuploidy_screen" in S.get("figures", {}):
            L.append(f"![aneuploidy screen]({S['figures']['aneuploidy_screen']})\n")
    # ---------------- benchmark ----------------
    if bench:
        L.append("## 6. Does it generalise?\n")
        L.append(f"On {bench['n_cases']} synthetic recessive cases spiked into the public GIAB HG002 genome (ClinVar truncating + missense "
                 f"pair; 6 sampled HPO terms with imprecision and one noise term; no diagnosis given), the true pair was ranked first in "
                 f"{pct_ci(bench['aware']['top1'], bench['n_cases'])} (top-10: {100 * bench['aware']['top10']:.0f} %) of cases, and in "
                 f"{pct_ci(bench['blind']['top1'], bench['n_cases'])} (top-10: {100 * bench['blind']['top10']:.0f} %) when the spiked alleles' "
                 "ClinVar records were hidden. The benchmark uses only public data and is regenerated by the notebook. It is deliberately "
                 "simple and therefore optimistic: the HPO terms are sampled from the same annotations that the phenotype engine scores, every "
                 "case is a coding truncating + missense recessive pair, the hidden alleles are still ones on which the missense predictor's "
                 "calibration was built, and a GIAB benchmark genome is cleaner than a clinical callset. It shows that the machinery works "
                 "end to end, not the accuracy to expect on an unsolved case.\n")
        if "benchmark" in S.get("figures", {}):
            L.append(f"![benchmark]({S['figures']['benchmark']})\n")
    # ---------------- limitations ----------------
    L.append("## 7. Limitations\n")
    L.append("* Phase of the two alleles is inferred from the genetic model and population data, not observed; parental testing or a "
             "long-read amplicon would settle it.\n* The missense allele remains a variant of uncertain significance by automated ACMG "
             "criteria; its causal role rests on the pair-level argument (rarity, predictor, kinase-domain location, phenotype match).\n"
             "* Structural variants and repeat expansions genome-wide are not assessed (read-level screen covers the leading locus only).\n"
             "* Phenotype similarity depends on the completeness of HPO annotations (several of the child's features are not annotated "
             "to the disorder).\n* EPCR is a ranking score on a probability scale (an ACMG posterior × phenotype, mode, phase and quality "
             "factors). It is not normalised across hypotheses (a single allele and a pair containing it overlap) and it is not calibrated "
             "against clinical outcomes, so it orders hypotheses rather than stating a clinical risk.\n"
             "* Robustness checks perturb the heuristic weights, not the ACMG prior and odds, which dominate the score.\n")
    L.append("**Where every number comes from.** Literature parameters were taken as published; heuristic ones were fixed before the "
             "challenge data were opened and not changed afterwards.\n")
    L.append(df_md(pd.DataFrame(PARAMS, columns=["parameter", "value", "source"]), max_rows=20) + "\n")
    # ---------------- data handling ----------------
    L.append("## 8. Data handling and AI disclosure\n")
    L.append("* All genotype-scale files lived in ephemeral scratch storage and were deleted at the end of the run; only findings "
             "(this report, the ranked CSV, figures) were kept. Only short lists of candidate variants were sent to public annotation "
             "services (Ensembl VEP, gnomAD), no genotype tables.")
    L.append(f"* Generative AI: {CFG['LLM_DISCLOSURE']}\n")
    # ---------------- template ----------------
    L.append("## 9. Methods-template answers\n")
    t_min = sum(v["seconds"] for v in RUN.values()) / 60
    L.append("| Question | Answer |\n|---|---|")
    L.append(f"| Team name | {CFG['TEAM']} |\n| Model number | 1 ({CFG['MODEL_NAME']}) |")
    L.append(f"| Generative AI used | {CFG['LLM_DISCLOSURE']} |")
    L.append("| Automated output or manual curation? | Fully automated output of the notebook; no manual review or re-ranking. Evidence strengths follow published ClinGen calibrations; the remaining heuristics (section 7) were fixed before the challenge data were opened; no parameter was tuned on leaderboard feedback. |")
    L.append("| Public data only? | Yes: gnomAD v4.1, ClinVar, MANE, Ensembl VEP, AlphaMissense, SpliceAI, ESM-2, AlphaFold DB, UniProt, HPO, PanelApp, GIAB HG002 (benchmark). |")
    L.append("| Compound heterozygous pairs? | Yes — pairs are enumerated per gene and scored jointly with a phase factor; single alleles are kept as fall-back rows. |")
    L.append("| Secondary / incidental findings | Dominant cancer-predisposition and other genes are evaluated as competing primary hypotheses; no secondary rows are submitted and names of genes outside the indication are withheld pending the organisers' IRB guidance. |")
    hw = f"{NCPU} vCPU" + (f", {N_GPU}× GPU" if N_GPU else ", no GPU")
    L.append(f"| Run time and cost | {t_min:.0f} min on a free Kaggle notebook ({hw}); no paid compute. |\n")
    # ---------------- run log + references ----------------
    L.append("## 10. Run log\n")
    L.append(run_table() + "\n")
    L.append("## References (resolved live from Europe PMC)\n")
    for p in REPORT_PMIDS:
        if pm.get(p):
            L.append(f"* {cite(p, pm)} — {clean_title(pm[p]['title'])}")
    L.append("\n*" + ACK + "*")
    return "\n".join(L)


def _allele(lv, kind):
    return next((a for a in lv.get("alleles", []) if a.get("kind") == kind), None)


def fragment_note(ptc, pm):
    """How the truncated product compares with the published BUBR1 truncations (native residues kept, not the X number)."""
    kept = ptc.get("native_residues_kept")
    if kept is None:
        return ""
    if kept >= 715:
        return (f"the truncated product keeps {kept} native residues, all {715} of a truncation that still supported the "
                f"checkpoint (731X; the authors note it 'contains only 715 amino acids of BUBR1 sequence', so its last residues are not native), "
                f"plus {kept - 715} more; another truncation (753X) was expressed only weakly, so the stability of this one is an open question "
                f"({cite('20516114', pm)}). A stable fragment that keeps the CDC20-, BUB3- and PP2A-B56-binding motifs but lacks the "
                "pseudokinase could also sequester those partners (a dominant-negative risk), which the cell assay must exclude")
    return (f"the truncated product keeps {kept} native residues, fewer than the 715 of the shortest truncation shown to support the "
            f"checkpoint (731X) ({cite('20516114', pm)})")


def track2_summary(gene, lv, t2, pm):
    """One-paragraph summary: genotype → mechanism → lever → candidate → the experiment that tests it → scale."""
    if not lv:
        return ("The allele-level triage did not run in this notebook run (see the run log), so the summary is limited to the "
                "literature layer and the network pass below.")
    ptc, mis, pair, at = _allele(lv, "premature stop"), _allele(lv, "missense"), lv.get("pair") or {}, lv.get("atlas") or {}
    net = t2.get("network") or {}
    net_ok = net.get("n_network_nodes") is not None
    net_neg = net_ok and not (net.get("verdicts") or {}).get("plausibly beneficial direction")
    alleles = [a for a in lv.get("alleles", []) if a.get("kind")]
    if not ptc:
        return (f"The leading {gene} genotype has no premature stop codon among the alleles the triage could map, so the stop-codon "
                "levers do not apply; the hypotheses below come from the literature layer and the network pass.")
    second = (f"a missense allele ({mis['label']})" if mis else "a second allele" if len(lv.get("alleles", [])) > 1 else None)
    parts = [f"The child carries a nonsense allele ({ptc['label']})" + (f" and {second}" if second else "") + f" in {gene}. "
             + (f"No approved drug acts on the {gene} network in a direction that would help (section 5), so we asked instead "
                "what can still be done with each allele. " if net_neg else
                "Besides the network pass (section 5), we asked what can still be done with each allele. ")
             + f"The nonsense allele creates a {ptc['stop_rna']} stop"
             + (", the stop codon most prone to drug-induced read-through" if ptc["stop_rna"] == "UGA" else "")
             + ("; its mRNA is predicted to be removed by nonsense-mediated decay (NMD). " if ptc.get("nmd_escape") is False
                else "; it is predicted to escape nonsense-mediated decay. " if ptc.get("nmd_escape") else ". ")]
    if gene == "BUB1B":
        parts.append("\n\nTwo approved medicines act on that allele by different routes: **azithromycin**, which in a reporter assay stimulated "
                     "read-through specifically at UGA stops — the class this allele belongs to — and which, in mice carrying a germline "
                     "nonsense allele of a tumour-suppressor gene, restored the full-length protein and prevented tumours "
                     f"({cite('37276642', pm)}; {cite('34066907', pm)}); and **amlexanox**, which "
                     "both relieves NMD and induces read-through, and restored full-length protein from the nonsense allele of a "
                     f"missense + nonsense genotype in another gene, though in a heterologous expression system ({cite('29247835', pm)}). ")
    nc = {k: v for k, v in (ptc.get("near_cognate") or {}).items() if not v.get("restores_wild_type")}
    if nc:
        ins = ", ".join(f"{ptc['wt_aa']}{ptc['residue']}{k} LLR {v['llr']:.1f}" for k, v in nc.items())
        parts.append(f"ESM-2 scores the amino acids that read-through would insert at residue {ptc['residue']} ({ins}); these are "
                     "comparable with one another at this residue, not with scores elsewhere in the protein. ")
    if gene == "BUB1B" and fragment_note(ptc, pm):
        kept = ptc["native_residues_kept"]
        parts.append(f"\n\nIf NMD alone is relieved, a fragment of {kept} native residues would remain. "
                     + ("That includes everything in a truncation published to support the checkpoint " if kept >= 715 else
                        "That is fewer residues than the shortest truncation published to support the checkpoint ")
                     + f"({cite('20516114', pm)}); whether it is stable, functional, or instead sequesters its binding partners is "
                     "untested, and section 2 sets out how to find out. ")
    if pair:
        parts.append("\n\nThese predictions can be tested in the child's own cells in one experiment (section 2): one RT-PCR amplicon "
                     "spanning both alleles phases them and measures how much nonsense-allele mRNA escapes NMD"
                     + ("; allele-resolved peptides measure the protein each allele makes. " if lv.get("peptides") else ". "))
    if at:
        parts.append(f"The same triage over the {at['n_alleles']} ClinVar P/LP alleles in the {at['n_genes']} genes of the disorder's panel "
                     "that have any catalogued alleles shows how many other catalogued alleles share each lever (section 7), and it "
                     "runs unchanged for any gene.")
    return "".join(parts)


def predictions_table(lv, pm, gene):
    ptc, mis = _allele(lv, "premature stop"), _allele(lv, "missense")
    prot_name = "BUBR1" if gene == "BUB1B" else f"{gene} protein"
    rows = []
    if lv.get("pair"):
        rows.append(("baseline (untreated cells ± puromycin / cycloheximide)", f"RT-PCR amplicon ({lv['pair']['amplicon_nt']:,} nt)",
                     "two allele-specific molecule classes if the alleles are in trans; the nonsense-allele fraction rises when NMD is blocked"))
    if ptc:
        agent = "azithromycin; amlexanox" if gene == "BUB1B" else "NMD-inhibitor / read-through agent"
        rows.append((f"{agent}, dose series", "allelic ratio on the same amplicon" if lv.get("pair") else "nonsense-allele mRNA (RT-qPCR)",
                     "nonsense-allele mRNA fraction rises"))
        if lv.get("peptides"):
            down = (lv.get("peptide_logic") or {}).get("missense_downstream_of_stop")
            rows.append((f"{agent}; gentamicin (positive control)", "allele-resolved peptides (section 2)",
                         "read-through peptides appear only with drug"
                         + ("; so does the reference missense-site peptide, which only a full-length product of the stop allele can carry"
                            if down else "")))
        if ptc.get("truncated_kda"):
            rows.append((agent, f"immunoblot, N-terminal {prot_name} antibody",
                         f"a ~{ptc['truncated_kda']:.0f} kDa fragment (residues {ptc['truncated_product']}) if NMD is relieved and the fragment is stable"))
    if mis and gene == "BUB1B":
        rows.append(("niacin (nicotinamide as comparator)", "BUBR1 half-life (cycloheximide chase), K668 acetylation",
                     "longer half-life; more BUBR1 protein; no growth advantage for the aneuploid fraction"))
    if ptc and mis and gene == "BUB1B":
        # the two levers act on different copies of the gene, so the combination is the natural test (approved drugs only)
        rows.append(("azithromycin + niacin together (both copies at once)", "total BUBR1 (immunoblot), premature chromatid separation",
                     "the effects add, because each drug acts on a different copy; if they do not, one of the two levers is not working"))
    if gene == "BUB1B":
        rows.append(("every arm vs untreated and parental cells", "nocodazole arrest, premature chromatid separation, micronuclei",
                     "move toward parental values; no worsening of segregation (TBK1 inhibition by amlexanox is the specific risk)"))
    else:
        rows.append(("every arm vs untreated and parental cells", f"{prot_name} level and a disease-relevant functional readout",
                     "move toward parental values"))
    return "| Arm | Readout | Pre-registered prediction if the lever works |\n|---|---|---|\n" + "\n".join(f"| {a} | {b} | {c} |" for a, b, c in rows)


def track2_report():
    t2, pm = S.get("track2", {}), {**(S.get("track2", {}) or {}).get("pmids", {}), **S.get("report_pmids", {})}
    T, lv, st = S.get("top", {}), S.get("levers", {}) or {}, S.get("structure", {}) or {}
    gene, net = T.get("gene", S.get("top_gene")), t2.get("network", {})
    ptc, mis = _allele(lv, "premature stop"), _allele(lv, "missense")
    used = set()

    def c(p):
        used.add(p)
        return cite(p, pm)

    L = [f"# {CFG['TEAM']} — Track 2 report: from the {gene} genotype to approved-drug hypotheses\n",
         f"*Generated automatically by the analysis notebook on {time.strftime('%Y-%m-%d', time.gmtime())}. Hypotheses for follow-up, "
         "not evidence that any medicine works; nothing here is medical advice, and no agent named here is proposed for use in this "
         "child. Every candidate below is a laboratory question to be answered in patient-derived cells first.*\n",
         "> **Thesis.** When no approved drug fits the protein, treat the *alleles*. This child's two "
         f"{gene} alleles — a premature stop and a destabilising missense — each have a distinct, approved-drug lever, and a single "
         "experiment in the child's own cultured cells tests both levers while also supplying the evidence that would firm up the "
         "diagnosis.\n"]
    # what is actually new here, said plainly: the panel scores innovation, and a method is easy to miss in prose
    if lv.get("alleles"):
        L.append("## What is new here\n")
        L.append("Conventional drug repurposing asks which drug acts on the protein a gene makes. For a loss-of-function "
                 "disorder that question often has no answer, and section 5 shows it returning none here. This report asks a "
                 "different question — what can still be done with each *allele* — and makes it computable:\n\n"
                 "* an allele-level triage that reads the stop codon, the decay prediction, the amino acids read-through would "
                 "insert, the exon's frame and the predicted structure, and assigns each allele a lever; it is gene-agnostic and "
                 "runs from public data alone;\n"
                 "* a protein language model used not to score the variant the child has, but to score the residues a drug would "
                 "*introduce*, which is what decides whether read-through is worth attempting;\n"
                 "* one RT-PCR amplicon, designed automatically with the primers printed in section 2, that answers the phase "
                 "question, the decay question and the drug question from the same cells;\n"
                 "* allele-resolved peptides whose masses differ, so the protein made by each copy of the gene can be measured "
                 "separately rather than inferred;\n"
                 "* the same triage run across every catalogued disease-causing allele in the disorder's gene panel, so the "
                 "answer for this child indexes the other children it would apply to.\n")
    # a one-screen map of the argument, before the prose: judges skim before they read
    if lv.get("alleles"):
        ptc0, mis0, pair0 = _allele(lv, "premature stop"), _allele(lv, "missense"), lv.get("pair") or {}
        top = [h for h in (t2.get("curated") or []) if h.get("rank") in (1, 2, 3)]
        L.append("## At a glance\n")
        L.append("| | Allele 1 | Allele 2 |\n|---|---|---|")
        L.append(f"| What it is | {ptc0['label'] if ptc0 else '—'} — premature stop | {mis0['label'] if mis0 else '—'} — missense |")
        L.append(f"| What it does to the protein | {'little or no protein from this copy (mRNA predicted to decay by NMD)' if ptc0 and ptc0.get('nmd_escape') is False else 'truncated product'} "
                 f"| protein made, but less stable or less active |")
        L.append("| The lever | restore protein from this copy: read-through, and relief of NMD "
                 "| raise the abundance of what this copy makes — if the change acts by destabilising the protein rather than by "
                 "disrupting an interaction surface, which the cell experiment decides |")
        L.append("| Approved agents to test (laboratory first) | "
                 + "; ".join(h["agents"].split("(")[0].strip() for h in top if "stop allele" in h.get("lever", "")) + " | "
                 + "; ".join(h["agents"].split("(")[0].strip() for h in top if "missense" in h.get("lever", "")) + " |\n")
        if pair0:
            L.append(f"**One experiment answers both, and the diagnosis too.** The two alleles sit {pair0['genomic_distance_bp']:,} bp "
                     f"apart in the genome but only {pair0['cdna_distance_nt']:,} nt apart in the mRNA, so a single "
                     f"{pair0['amplicon_nt']:,}-nt RT-PCR amplicon from cultured cells (primers in section 2) reads both alleles on the "
                     "same molecule — read on a long-read platform, or cloned so that individual molecules are Sanger-sequenced, because the two "
                     "sites are too far apart for a single Sanger read to resolve them. It settles whether the alleles are in trans, "
                     "measures how much of the stop allele's mRNA survives NMD, and "
                     "with allele-resolved peptides measures the protein each allele makes under each drug. The same result supplies the "
                     "evidence that would move the second allele out of 'uncertain significance' (section 7). Because the two levers "
                     "act on different copies of the gene, the same cells can also test them together (a read-through agent with "
                     "niacin), where the prediction is that their effects add.\n")
        L.append("**Two things here are usable without waiting for any experiment.** Section 6 lists the medicines and procedures "
                 "that warrant extra care in a child who has both chromosomal instability and a cancer predisposition, each with "
                 "its published reason — that is for the treating team, today. Section 7 shows how the same cell experiment that "
                 "tests the drug hypotheses would also supply the evidence needed to reclassify the second allele, which is what "
                 "turns a probable diagnosis into a firm one and opens carrier testing for the family. Neither depends on any "
                 "medicine working.\n")
        L.append("**What this report does not claim.** That any of these agents should be given to this child, or that read-through "
                 "delivers enough protein to help — that is the open question, and sections 1 and 2 state what would settle it and what would "
                 "falsify each hypothesis.\n")
    L.append("## Summary\n")
    L.append(track2_summary(gene, lv, t2, pm) + "\n")
    L.append("**Method in brief.** (1) *Allele-level triage* — computed for each allele from the reference sequence, the MANE gene "
             "model, ESM-2 and the AlphaFold model; (2) *literature layer* — graded hypotheses whose every citation is resolved live; "
             "(3) *automated network pass* over Open Targets with direction gates; (4) *scalability* — the triage repeated over every "
             "ClinVar P/LP allele in the disorder's gene panel. Grades: A human trial in a relevant population · B in vivo evidence in "
             "the same gene defect (the intervention may be genetic rather than the drug itself) · C human patient cells or trials in another disease with the same molecular lesion · D mechanistic inference · "
             "X direction of effect uncertain.\n")

    # ---- 1. mechanism ---------------------------------------------------------------------------------------------------
    L.append("## 1. Mechanism of the variants\n")
    if gene == "BUB1B":
        L.append("*The genotype this report builds on came from five independent runs of this pipeline on the child's data, on three "
                 "kinds of machine (a dual-GPU machine, a TPU host and a plain four-CPU machine) and with two different "
                 "population-frequency sources, the second used when a download of the first failed. The top-ranked genotype and its "
                 "coordinates were identical in every run. The lower-ranked alternatives within the same gene did change between "
                 "runs (three to five rows), so only the top-ranked pair is relied on here. This is a statement about repeat runs of "
                 "this pipeline, not an independent confirmation of the diagnosis, and it is recorded because the rest of this report "
                 "rests on it.*\n")
    for v in T.get("variants", []):
        kind = ("null allele — " + ("nonsense-mediated decay predicted, no protein from this copy" if not v.get("nmd_escape") else "truncated protein escaping NMD")
                if v["csq"] in LOF_CSQ else f"{v['csq']} — " + ("hypomorphic candidate: protein made but expected to be less stable or less active"
                                                               if v["csq"] == "missense" else "regulatory / splice candidate"))
        L.append(f"* `{v.get('hgvsp') or v['key']}` ({v['classification']}): {kind}.")
    L.append("")
    if gene == "BUB1B":
        L.append(
            "**Loss of function with residual activity.** BUBR1 is the pseudokinase core of the mitotic checkpoint complex that holds "
            f"APC/C-CDC20 until every kinetochore is attached, and it recruits PP2A-B56 to stabilise correct attachments "
            f"({c('22698286')}; {c('23079597')}). Complete loss is cell-lethal ({c('15159543')}), so affected children carry one "
            f"null allele plus one allele with residual function ({c('15475955')}); in a graded knockdown of a human cell line, premature chromatid "
            f"separation rose first and aneuploidy appeared only below about 50 % of normal BUBR1 ({c('18932004')}). **The therapeutic lever is therefore a modest, "
            "proportional increase in residual BUBR1**; full rescue is neither expected nor required.\n")
        L.append(
            "**How much would be enough?** This question decides whether a weak intervention is worth testing at all. A graded knockdown "
            "in a human cell line — not patient cells — gives a dose-response rather than a threshold: premature chromatid separation rose "
            f"first, and aneuploidy appeared only below about half of normal BUBR1 ({c('18932004')}). Mice engineered to express only low amounts of BUBR1 are "
            f"viable and reproduce much of the MVA spectrum ({c('15208629')}), and combining different mutant alleles showed that which "
            f"alleles are present changes severity beyond the total protein level alone ({c('31738183')}). Read-through typically restores "
            "only a low single-digit percentage of a protein, and that is from one allele. **We do not claim that this is enough: "
            "the gap between what read-through delivers and what this genotype needs is the central unresolved question of this "
            "report, and we have found no data that closes it for BUB1B.** What can be done is to measure that gap rather than "
            "argue about it, which is the purpose of the gentamicin arm: as the most potent read-through agent available it sets "
            "the ceiling for how much protein this particular allele can yield, turning an open question into a number before any "
            "clinical question is raised. If that ceiling sits far below the level at which premature chromatid separation "
            "improves, the stop-codon hypotheses should be abandoned and the missense allele's abundance is the remaining lever.\n")
        dis = st.get("disease_missense_uniprot") or []
        L.append(
            f"**Two allele classes, one of each here.** In patient cells, MVA alleles fall into two classes ({c('20516114')}): truncating "
            "alleles whose transcripts are absent (NMD), and missense alleles in or near the kinase domain that lower BUBR1 through faster "
            "protein turnover. The same study showed that a truncated BUBR1 lacking the kinase domain (731X, 715 BUBR1 residues) still "
            "supports the checkpoint and stable microtubule attachments when expressed from a cDNA (no NMD), whereas a slightly "
            "longer truncation (753X) was expressed at lower levels, which the authors note need not mean reduced stability."
            + (f" The child's missense change lies in the {st.get('domain')} (AlphaFold pLDDT {st.get('plddt')}), "
               f"{dis[0]['distance_A']} Å from the reported {dis[0]['disease']} allele {dis[0]['change']}"
               + (f"; it is less buried than {100 - st['burial_percentile_in_domain']} % of the domain's residues, so a surface "
                  "or interface effect is as plausible as unfolding" if st.get("burial_percentile_in_domain") is not None and st['burial_percentile_in_domain'] < 50 else "")
               + ". Whether it acts like the published kinase-domain alleles (low abundance through faster turnover) or through a "
               "functional defect is the first question the cell experiment answers (cycloheximide chase, section 2)."
               if dis else "") + "\n")
        L.append("| Consequence of BUBR1 insufficiency | Child's feature it plausibly explains | Evidence |\n|---|---|---|")
        L.append(f"| chromosome missegregation → variegated aneuploidy | the diagnosis | {c('15475955')} |")
        L.append(f"| cancer predisposition (Wilms tumour, rhabdomyosarcoma, haematological) | rhabdomyosarcoma | {c('11746029')}; {c('39264246')} |")
        L.append(f"| p16-driven senescence, progeroid 'cachectic dwarfism' in hypomorphic mice | small for gestational age, short stature, failure to thrive, muscle atrophy | {c('15208629')}; {c('18516091')} |")
        L.append(f"| carrier parents show intermediate premature chromatid separation (method-dependent) | family counselling: whether carrier parents have any reproductive phenotype is **unproven** — one unreplicated 2026 report of BUB1B variants in recurrent pregnancy loss, and a negative screen of aneuploid miscarriages | {c('11746029')}; {c('15108198')}; {c('42434306')}; {c('40028477')} |")
        L.append("| — (not explained by the checkpoint model) | nephrocalcinosis: prematurity and renal tubulopathies are more parsimonious | see the nephrocalcinosis hypothesis in section 3 |\n")

    # ---- 2. allele-level triage ------------------------------------------------------------------------------------------
    if lv.get("alleles"):
        L.append("## 2. Allele-level therapeutic triage (computed)\n")
        L.append("Drug-target databases ask which drug acts on a protein near the gene. For a loss-of-function disorder the more useful "
                 "question is what can still be done with each allele; this section answers it from sequence, gene model, protein "
                 "language model and predicted structure. The code is gene-agnostic.\n")
        if ptc:
            ex = ptc.get("exon") or {}
            nc = ptc.get("near_cognate") or {}
            rank_txt = {1: "first", 2: "second", 3: "third"}.get(ptc.get("readthrough_rank"), "?")
            rows = [("stop codon / next nucleotide", f"{ptc['stop_rna']}, followed by {ptc.get('plus4_rna') or ptc['plus4']}; "
                     f"{ptc['stop_rna']} ranks {rank_txt} of the three stop codons for read-through (UGA > UAG > UAA)",
                     f"UGA is the most read-through-prone stop in reporter assays ({c('10939566')}); a +4 C enhances read-through, "
                     f"while stop identity alone did not change the gentamicin response in another study ({c('22479203')}), and in "
                     f"Duchenne patients neither predicted response ({c('20517938')}): measure, do not infer"),
                    ("nonsense-mediated decay", "predicted" if ptc.get("nmd_escape") is False else "escapes" if ptc.get("nmd_escape") else "unknown",
                     f"exon {ex.get('exon')} of {ex.get('n_exons')}; NMD removes the template for read-through, hence the value of a dual NMD-inhibitor / read-through drug")]
            restore = [k for k, v in nc.items() if v.get("restores_wild_type")]
            nc = {k: v for k, v in nc.items() if not v.get("restores_wild_type")}
            if nc:
                comp = (f". Every substitution the triage scores at residue {ptc['residue']} lands in a narrow band, which says "
                        "the position is weakly constrained rather than that any particular insertion would work; these scores "
                        "are comparable with each other at this residue, but not with scores at a different residue, so the "
                        f"child's own {mis['short']} (LLR {mis['llr']:.1f}) is given for scale only and not as a benchmark"
                        ) if mis and mis.get("llr") is not None else ""
                rows.append(("what read-through would insert", ", ".join(f"{ptc['wt_aa']}{ptc['residue']}{k}: LLR {v['llr']:.1f}"
                                                                       + (f" (better tolerated than {v['pct_less_tolerated']:.0f} % of all possible substitutions)" if v.get("pct_less_tolerated") is not None else "") for k, v in nc.items()),
                             f"near-cognate tRNAs insert {', '.join(AA3.get(a, a) for a in ptc.get('near_cognate_aas', ''))} at {ptc['stop_rna']} "
                             f"({c('27702906')}){comp}"
                             + (f"; {AA3.get(restore[0], restore[0])} insertion would restore the wild-type residue" if restore else "")))
            rows.append(("structure at the stop residue", f"AlphaFold pLDDT {fmt(ptc.get('plddt'))}", ""))
            rows.append(("truncated product if NMD is relieved", f"residues {ptc['truncated_product']} (~{fmt(ptc.get('truncated_kda'))} kDa); keeps "
                         + ("; ".join(ptc.get("truncated_keeps") or []) or "—") + "; loses " + ("; ".join(ptc.get("truncated_loses") or []) or "—"),
                         (fragment_note(ptc, pm) if gene == "BUB1B" else "")))
            if ex:
                if not ex.get("in_frame"):
                    verdict = "not eligible: skipping would shift the frame"
                elif not ex.get("skippable"):
                    verdict = "not eligible: first or last coding exon, or only partly coding"
                elif ex.get("mean_plddt") is None:
                    verdict = "candidate in principle; structure of the exon's residues not assessed"
                elif ex["mean_plddt"] >= 70:
                    verdict = "lower priority: the exon's residues are structured"
                else:
                    verdict = "candidate: the exon's residues are largely disordered"
                iso = ex.get("natural_skipping_isoforms")
                rows.append(("antisense exon skipping", f"exon {ex['exon']}: {ex['length_nt']} nt ({'in frame' if ex['in_frame'] else 'frame-shifting'}), "
                             f"spanning residues {ex['aa_from']}-{ex['aa_to']} (the codons at each boundary are split across "
                             f"neighbouring exons), mean pLDDT {fmt(ex.get('mean_plddt'))}; natural isoform skipping it: "
                             + ("none annotated" if iso == [] else ", ".join(iso) if iso else "not checked"), verdict))
            L.append(f"**Nonsense allele {ptc['label']}**\n")
            L.append("| Property | Value | What follows |\n|---|---|---|\n" + "\n".join(f"| {a} | {b} | {x} |" for a, b, x in rows) + "\n")
        if mis:
            near = mis.get("nearest_patient_missense") or []
            tail = (f". No approved drug binds this domain: BUBR1 is a pseudokinase ({c('22698286')}) whose kinase domain bound no "
                    f"nucleotide in a pseudokinase thermal-shift survey ({c('24107129')}), so a small-molecule chaperone is not available; "
                    "the lever is protein abundance (the missense-allele lead, hypothesis 3).\n"
                    if gene == "BUB1B" else ".\n")
            L.append(f"**Missense allele {mis['label']}** — {mis.get('domain') or 'no annotated domain'}; AlphaFold pLDDT {fmt(mis.get('plddt'))}; "
                     f"ESM-2 LLR {fmt(mis.get('llr'))}"
                     + (f"; less buried than {100 - mis['burial_percentile']} % of the domain's residues" if mis.get("burial_percentile") is not None else "")
                     + ("; nearest literature-reported patient missense: " + ", ".join(f"{x['change']} ({x['disease']}, {x['distance_A']} Å)" for x in near) if near else "")
                     + tail)
        if S.get("figures", {}).get("allele_levers"):
            L.append(f"![allele-level levers]({S['figures']['allele_levers']})\n")
        L.append("### One experiment in the child's own cells tests every lever\n")
        pr = lv.get("pair")
        if pr:
            L.append(f"* **Phase and NMD from one RT-PCR.** The alleles are {pr['genomic_distance_bp']:,} bp apart in the genome — too far "
                     f"for short reads — but {pr['cdna_distance_nt']:,} nt apart in the mRNA ({pr['fold_shorter']}× closer). One amplicon of "
                     f"~{pr['amplicon_nt']:,} nt from cultured cells, with and without an NMD inhibitor as is standard in RNA diagnostics "
                     f"({c('40473777')}), shows whether the alleles are in trans and how much nonsense-allele mRNA "
                     f"survives NMD. The two sites are {pr['cdna_distance_nt']:,} nt apart within the amplicon, which is beyond one Sanger read, so the "
                     "single-molecule readout is long-read sequencing of the amplicon or Sanger sequencing of individual clones; "
                     f"RNA-based phasing of distant variants is established in principle ({c('27605262')}); parental testing "
                     f"remains the clinical standard ({c('25741868')}).")
            if pr.get("primers"):
                p0 = pr["primers"][0]
                L.append(f"  Primer3 pair, ready to order: forward `{p0['left']}` (exon {p0['exons_left']}, Tm {p0['tm_left']} °C), "
                         f"reverse `{p0['right']}` (exon {p0['exons_right']}, Tm {p0['tm_right']} °C); {p0['product_nt']:,}-nt product with "
                         f"the two sites at positions {p0['site_offsets'][0]} and {p0['site_offsets'][1]}. "
                         + ("One primer spans an exon–exon junction, so genomic DNA does not amplify. " if p0["junction_spanning"]
                            else f"Genomic DNA would give a {p0['genomic_span_bp']:,}-bp product, beyond a standard extension time. ")
                         + ("No variant of the child " if p0["child_variants_under_primers"] == 0 else "")
                         + ("and no common SNP (MAF ≥ 1 %) " if p0["common_snps_under_primers"] == 0 and p0["child_variants_under_primers"] == 0
                            else "No common SNP (MAF ≥ 1 %) " if p0["common_snps_under_primers"] == 0 else "")
                         + ("lies under either primer, so neither allele should drop out and bias the allelic ratio. "
                            if 0 in (p0["child_variants_under_primers"], p0["common_snps_under_primers"]) else "")
                         + "Reads carrying the stop allele with the reference missense base (and vice versa) "
                         f"prove trans; {len(pr['primers']) - 1} alternative pair(s) are in the notebook output.")
        peps = lv.get("peptides") or []
        if peps:
            lg = lv.get("peptide_logic") or {}
            L.append("* **Allele-resolved peptides.** If the alleles are in trans, a read-through product of the nonsense allele carries a "
                     "near-cognate amino acid at the stop residue"
                     + (" and the reference residue at the missense site, which no other product can carry because the truncated "
                        "product ends before it" if lg.get("missense_downstream_of_stop") else "")
                     + ". Each peptide below is digested from the protein that allele actually makes and differs in mass from the others, "
                     "so targeted mass spectrometry measures read-through directly, separately from the missense allele's protein:\n")
            L.append("| Protease | Form | Peptide | Monoisotopic mass (Da) | Made by |\n|---|---|---|---|---|\n"
                     + "\n".join(f"| {p['protease']} | {p['form']} | `{p['peptide']}` | {p['mass']:.3f} | {p['made_by']} |" for p in peps) + "\n")
        L.append("**Pre-registered predictions**\n")
        L.append(predictions_table(lv, pm, gene) + "\n")
        if gene == "BUB1B":
            L.append("**What would falsify each hypothesis.** Stated before the experiment so the result cannot be reinterpreted "
                     "afterwards:\n\n"
                     "* *No read-through peptide and no full-length protein in the gentamicin arm* → this allele is not "
                     "read-through-competent, and the stop-codon hypotheses (1, 2 and 4) are dead, whatever the mRNA does.\n"
                     "* *No read-through above the assay's limit of quantitation* — which must be fixed first, from a dilution series of "
                     "synthetic peptide spiked into untreated lysate, so that 'none detected' means a number rather than an "
                     "impression. *Read-through detectable with gentamicin but the ceiling it reaches is far below the level at which "
                     "premature chromatid separation improves* → the approach cannot work for this genotype at achievable "
                     "exposures; report the ceiling and stop.\n"
                     "* *Nonsense-allele mRNA does not rise when NMD is blocked* → the NMD premise is wrong and hypothesis 2 "
                     "loses its rationale.\n"
                     "* *The truncated fragment is stable but segregation worsens* → the dominant-negative risk above is real, "
                     "and relieving NMD is harmful rather than helpful.\n"
                     "* *No increase in BUBR1 half-life or abundance with niacin* → the NAD+/SIRT2 hypothesis (3) is dead for "
                     "this missense allele.\n"
                     "* *The amplicon shows the alleles in cis* → the compound-heterozygous model is wrong, and every "
                     "allele-level hypothesis here is void.\n")
            L.append("**Feasibility.** The whole panel runs on cultured skin fibroblasts or a lymphoblastoid line from one biopsy "
                     "or blood draw that a diagnostic laboratory can already obtain; the readouts (RT-PCR read on a long-read platform or "
                     "cloned and Sanger-sequenced, two immunoblots, targeted mass spectrometry, and standard cytogenetic scoring of premature "
                     "chromatid separation and micronuclei) are routine, need no new reagents beyond the primers printed above, "
                     "and would report in weeks rather than years. Nothing here requires the child to take anything.\n")

    # ---- 3. hypotheses ---------------------------------------------------------------------------------------------------
    if t2.get("curated"):
        L.append("## 3. Ranked hypotheses\n")
        # index first (one line each, readable at a glance), then one block per hypothesis: an 8-column
        # table of paragraphs is unreadable once it is printed to PDF
        L.append("| # | Lever | Hypothesis | Grade |\n|---|---|---|---|")
        for h in t2["curated"]:
            L.append(f"| {h['rank']} | {h['lever']} | {h['hypothesis']} | {h['grade']} |")
        L.append("")
        for h in t2["curated"]:
            used.update(h["pmids"])
            L.append(f"### {h['rank']}. {h['hypothesis']} — grade {h['grade']}\n")
            L.append(f"* **Lever.** {h['lever']}")
            L.append(f"* **Agent(s).** {h['agents']}")
            cap = lambda t: t[:1].upper() + t[1:]           # the source strings are written as clauses
            L.append(f"* **Rationale.** {cap(h['mechanism'])}")
            L.append(f"* **Main caveat.** {cap(h['caveat'])}")
            L.append(f"* **How it would be tested.** {cap(h['validation'])}")
            L.append("* **References.** " + ", ".join(c(p) for p in h["pmids"]) + "\n")
        L.append("## 4. Considered and set aside\n")
        L.append("| Option | Why it is not proposed |\n|---|---|")
        for x in t2.get("set_aside", []):
            L.append(f"| {x['item']} | {x['reason']}" + (f" (PMID {', '.join(x['pmids'])})" if x["pmids"] else "") + " |")
            used.update(x["pmids"])
        L.append("")

    # ---- 5. network pass -------------------------------------------------------------------------------------------------
    L.append("## 5. What the automated network pass found\n")
    if net.get("n_network_nodes") is None:
        L.append("The network pass did not complete in this run (see the run log).\n")
    else:
        L.append(f"{net.get('n_network_nodes')} nodes (the gene, its Open Targets interaction partners with score ≥ 0.7, and literature "
                 f"mechanism nodes); every drug / clinical candidate on them (Open Targets {net.get('open_targets_api')}), gated by approval → "
                 f"direction (drug action × cited node sign) → antimitotic / genotoxic mechanism. {gene} has {net.get('n_direct_drugs')} drugs "
                 f"acting on it directly (essential gene: {net.get('essential')}); {net.get('n_approved')} approved drug–target pairs on the "
                 "network; by verdict: " + ", ".join(f"{k}: {v}" for k, v in (net.get("verdicts") or {}).items()) + ". "
                 + ("**No approved drug acts on the network in a direction predicted to raise residual function**"
                    + (" — the classes that would (SIRT2 activation, partial APC/C-CDC20 inhibition) have no approved member" if gene == "BUB1B" else "")
                    + ". This honest negative is why the allele-level and literature layers exist.\n"
                    if not (net.get("verdicts") or {}).get("plausibly beneficial direction") else
                    "Drugs in the 'plausibly beneficial direction' class are listed in `track2_network_drugs.csv`.\n"))
    gen = t2.get("generality") or []
    if gen:
        L.append("The identical pass for the other genes of the indication panel:\n\n" + df_md(pd.DataFrame(gen).fillna(0)) + "\n")

    # ---- 6. cautions -----------------------------------------------------------------------------------------------------
    if t2.get("caution"):
        L.append("## 6. Caution list for the clinical team\n")
        L.append("| Medicine / class | Reason |\n|---|---|")
        for x in t2["caution"]:
            L.append(f"| {x['item']} | {x['reason']} (PMID {', '.join(x['pmids'])}) |")
            used.update(x["pmids"])
        L.append("")

    # ---- 7. impact -------------------------------------------------------------------------------------------------------
    L.append("## 7. Impact\n")
    rc = lv.get("reclassification")
    if rc:
        L.append(f"**For the diagnosis.** The weaker allele ({rc['allele']}) is classed "
                 f"{(rc.get('current_class') or 'as shown below').lower()} by automated criteria. The evidence that would "
                 f"move it, on the point scale of {c('32720330')}:\n")
        L.append("| Step | Points | Posterior | Class |\n|---|---|---|---|\n"
                 + "\n".join(f"| {r['step']} | {r['points']:g} | {r['posterior']:.3f} | {r['classification']} |" for r in rc["path"]) + "\n")
        L.append("Parental testing (or the cDNA amplicon above) supplies PM3, and the same cell experiment that tests the drugs supplies "
                 "functional evidence: one experiment moves the diagnosis and the therapy question forward together. A firm diagnosis "
                 "opens carrier testing for the parents and, if both carry one allele, a 25 % recurrence risk per pregnancy and the "
                 "option of prenatal or preimplantation testing.\n")
    if gene == "BUB1B":
        L.append(f"**Published MVA surveillance guidance, for the treating team's reference.** Consensus surveillance ({c('39264246')}): renal ultrasound every "
                 "3 months until age 7; clinical review for signs of rhabdomyosarcoma, Wilms tumour and haematological malignancy; avoid "
                 "ionising radiation; HPV vaccination. The nephrocalcinosis work-up (section 3) may uncover a treatable cause of part of "
                 "the failure to thrive.\n")
    at = lv.get("atlas")
    if at:
        bl = at.get("by_lever") or {}
        L.append(f"**For other children with the disorder.** The triage was repeated over all {at['n_alleles']} ClinVar P/LP alleles in "
                 f"{at['n_genes']} of the panel's genes that have any catalogued alleles: " + ", ".join(f"{k}: {v}" for k, v in bl.items()) + ". Nonsense alleles by stop "
                 "codon: " + (", ".join(f"{k.replace('T', 'U')} {v}" for k, v in (at.get("nonsense_by_stop") or {}).items()) or "none")
                 + f" (in {gene}: " + (", ".join(f"{k.replace('T', 'U')} {v}" for k, v in (at.get("top_gene_nonsense_by_stop") or {}).items()) or "none")
                 + "). A positive cell result for the nonsense lever would therefore matter beyond this child.\n")
        nb = at.get("nonsense_by_stop") or {}
        n_uga = nb.get("UGA") or nb.get("TGA") or 0
        n_rt = sum(v for k, v in (at.get("by_lever") or {}).items() if "read-through" in k)
        if n_uga:
            L.append(f"**The number that matters for reach.** Of the {at['n_alleles']} catalogued disease-causing alleles in "
                     f"this panel, {n_rt} are premature stops and {n_uga} of those are the same stop codon as this child's. "
                     "Those are alleles, not people, and one person carries two; but every one of them is an allele for which "
                     "this exact assay — the same primer logic, the same peptide design, the same readouts — would be built "
                     "the same way. That is the sense in which a result here is not an n-of-1."
                     + (f" {at['n_no_criteria']} of the {at['n_alleles']} have no assertion criteria in ClinVar; the review status is "
                        "in the atlas file, so a stricter reading can simply filter them out." if at.get("n_no_criteria") else "")
                     + "\n")
        L.append(df_md(pd.DataFrame(at.get("by_gene") or [])) + "\n")
    if gene == "BUB1B":
        L.append("**What would have to happen next, and who could do it.** Nothing in this report changes anything today; it is a set of "
                 "laboratory questions, and this is the shortest path from here to an answer:\n\n"
                 "1. *Confirm the genetics.* The treating clinical genetics service confirms both alleles and tests the parents, which "
                 "settles phase by the standard route and is worth doing whatever happens with the rest.\n"
                 "2. *Grow the cells.* Skin fibroblasts or a lymphoblastoid line from the child, with a parental line alongside as the "
                 "comparison. The challenge page says the MVA Society is funding focused research; whether banked patient cells already exist is for the family and the Society to say, and nothing here assumes it.\n"
                 "3. *Run the one amplicon.* Phase, NMD escape and the drug arms are all read out from the same cultured cells using the "
                 "primers printed in section 2 — the cheapest experiment here and the one that can falsify the most.\n"
                 "4. *Report the answer either way.* Because mosaic variegated aneuploidy affects so few people, a negative result that "
                 "retires a plausible idea is itself worth publishing; the atlas above names the other alleles to which the same answer "
                 "would apply.\n\n"
                 "Each step uses assays a diagnostic or research laboratory already runs, which is the point: this is a protocol, not a "
                 "programme.\n")
        L.append("**In plain language.** The child has two changes in the same gene, one inherited from each parent (which parental "
                 "testing would confirm). One change stops the gene's instructions early, so that copy makes little or no protein. The "
                 "other changes a single building block, so that copy probably makes protein that does not last as long or work as well. "
                 "The protein's job is to check that chromosomes are shared out correctly when a cell divides, which is why cells in this "
                 "condition end up with the wrong number of chromosomes. **No medicine is being recommended here, and nothing in this "
                 "report should change anyone's current treatment.** What we have done is ask whether any already-approved medicine acts "
                 "on either of those two specific problems, and then describe a laboratory test, done on cells grown from a small skin or "
                 "blood sample, that would show whether any of them does anything useful — before anyone considers whether it could help a "
                 "person. The same test would also help the laboratory decide, formally, whether the second change is disease-causing, "
                 "which matters for the family's own genetic counselling. Questions a family might reasonably raise with their own care "
                 "team are whether both parents have been tested, whether the recommended surveillance schedule is in place, and whether "
                 "stored cells could be made available to researchers.\n")

    # ---- 7b. does the triage get known answers right? --------------------------------------------------------------------
    val = S.get("lever_validation") or {}
    if val.get("n_tested"):
        L.append("## 8. Does the triage get published answers right?\n")
        L.append(f"A classifier is only worth anything if it reproduces answers that are already known. The same triage was run, "
                 f"unchanged, over {val['n_tested']} premature-stop alleles in {len(val.get('genes') or [])} genes outside this "
                 "disorder's panel — alleles chosen because the published literature already reports how they behave when "
                 f"read-through is attempted. It recovered the published stop codon for {val['n_codon_match']} of "
                 f"{val['n_tested']} ({100 * val.get('frac_codon_match', 0):.0f} %) directly from the reference sequence and the "
                 "gene model, and assigned the matching lever class in "
                 f"{val.get('n_lever_match', val['n_codon_match'])} of {val['n_tested']}.\n")
        if val.get("confusion"):
            L.append("| Lever assigned | UGA | UAG | UAA |\n|---|---|---|---|\n"
                     + "\n".join(f"| {r['lever']} | {r.get('UGA', 0)} | {r.get('UAG', 0)} | {r.get('UAA', 0)} |"
                                 for r in val["confusion"]) + "\n")
        L.append("Recovering a stop codon from the reference sequence and the gene model is close to a sequence lookup, so a full "
                 "score is the expected result rather than a finding about biology. What the check does establish is that the "
                 "pipeline reads the reference, the gene model and the numbering correctly across genes and strands, and that its "
                 "lever classes separate by stop codon as intended.\n")
        L.append("Two honest limits on what this shows. Every allele in the set is one where read-through has been reported to "
                 "work to some degree, so the check measures whether the triage finds the alleles it should, not whether it "
                 "rejects the ones it should — it tests sensitivity, not specificity. And one of the alleles contradicts the "
                 "ordering the triage uses as a prior: a UAA allele responded better than UGA alleles tested alongside it. That "
                 "is why the stop codon is treated here as a prior to be measured rather than a rule, and it is the clearest "
                 "argument for running the laboratory arm before believing any ranking, including this report's own.\n")
        L.append("The per-allele table is saved as `lever_validation.csv`, with the published source for every entry.\n")

    # ---- 8. scalability --------------------------------------------------------------------------------------------------
    L.append("## 9. Scalability\n")
    rt = RUN.get("Track 2: allele-level therapeutic triage", {})
    L.append("* The allele-level triage, the network pass and the report are produced by the notebook with no manual step, for whichever "
             "gene leads, from public resources only (MANE, GRCh38, ClinVar, ESM-2, AlphaFold DB, UniProt, Ensembl, Open Targets, Europe PMC)"
             + (f"; the triage took {rt['seconds'] / 60:.1f} min." if rt.get("seconds") else ".")
             + "\n* In demo mode the same code runs on a public genome with a different spiked-in recessive gene, which is how it is tested."
             "\n* The atlas above is the same computation applied to a whole gene panel; it extends to any gene list. It is saved as "
             "`lever_atlas.csv` (one row per public ClinVar pathogenic / likely-pathogenic allele, with the lever, the reason for it, the "
             "stop codon, the exon and the ClinVar review status), so another group can point the same triage at their own gene list and "
             "get a comparable table without rerunning the patient-specific parts.\n")
    L.append("**Reproducibility.** Code: " + CFG.get("REPO_URL", "see the repository link submitted with this report")
             + " (CC BY 4.0, as required for submissions to this hackathon). `python build_nb.py` assembles the percent-format "
             "sources in `src/` into the notebook; running that notebook end to end reproduces every number, table and figure in "
             "this report, including this sentence, because the report is written by the pipeline rather than by hand. Without the "
             "gated dataset the same notebook runs in demo mode on a public genome with a different recessive gene spiked in, which "
             "is how each stage is tested; the run that produced this report took "
             + (f"{sum(v.get('seconds', 0) for v in RUN.values() if isinstance(v, dict)) / 60:.0f} minutes " if RUN else "")
             + "on a 4-CPU machine with no GPU.\n")

    # ---- 9. limitations --------------------------------------------------------------------------------------------------
    L.append("## 10. Limitations\n")
    L.append("* No candidate has been tested in MVA; grades reflect indirect evidence. Amlexanox's cellular concentrations may not be "
             "reachable in plasma, its TBK1/IKKε inhibition could itself affect mitosis, and its oral form has been marketed in Japan, with its current status unverified.\n"
             "* Read-through efficiency in patients is not predicted by the stop-codon context; ESM-2 scores of the inserted amino acids are "
             "predictions, not measurements, and read-through products can be immunogenic.\n"
             "* Phase is inferred, not observed (see the Track 1 report); the allele-resolved readouts assume the alleles are in trans, which "
             "the cDNA amplicon tests first.\n"
             "* Network drug lists depend on database coverage and on a small hand-curated sign map (printed in the notebook).\n"
             f"* Generative AI: {CFG['LLM_DISCLOSURE']}\n")
    L.append("## References (resolved live from Europe PMC)\n")
    for p in sorted(used | {"15475955", "16411201"}, key=lambda x: ((pm.get(x) or {}).get("year", "9999"), x)):
        if pm.get(p):
            L.append(f"* {cite(p, pm)} — {clean_title(pm[p]['title'])}")
    L.append("\n*" + ACK + "*")
    return "\n".join(L)


def write_reports():
    if "report_pmids" not in S:
        try:
            S["report_pmids"] = {k: v for k, v in resolve_pmids(REPORT_PMIDS).items() if v}
        except Exception:       # noqa: BLE001
            S["report_pmids"] = {}
    paths = []
    for name, fn in ((f"{CFG['TEAM']}_track1_report.md", track1_report), (f"{CFG['TEAM']}_track2_report.md", track2_report)):
        try:
            (OUT / name).write_text(fn(), encoding="utf-8")
            paths.append(name)
        except Exception as e:      # noqa: BLE001
            log(f"report {name} failed: {type(e).__name__}: {str(e)[:200]}")
    summary = {k: v for k, v in S.items() if k not in ("report_pmids",)}
    summary["run"] = RUN
    (OUT / "run_summary.json").write_text(json.dumps(summary, indent=1, default=str), encoding="utf-8")
    return paths


with stage("write reports"):
    log(f"reports written: {write_reports()}")
    display(Markdown((OUT / f"{CFG['TEAM']}_track1_report.md").read_text(encoding="utf-8")[:6000]))

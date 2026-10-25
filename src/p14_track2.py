# %% [markdown]
# ## 14 · Track 2 — from mechanism to approved-drug hypotheses
#
# Two layers, both reproducible from public sources:
#
# 1. **Automated network pass (gene-agnostic).** The causal gene's interaction neighbourhood (Open Targets, score ≥ 0.7,
#    cross-checked with STRING) is queried for every drug and clinical candidate. Hard gates are applied in order:
#    *approved* (maximum clinical stage = approval) → *direction* (does the drug push its target in the direction that
#    worsens or relieves a loss-of-function / hypomorphic state?) → *antimitotic or genotoxic mechanism* (a separate caution
#    class in a chromosomal-instability disorder). Direction comes from the drug's action type (ChEMBL / Open Targets)
#    combined with a short, cited sign map of the network nodes — the only hand-curated input, printed in full.
# 2. **Literature-curated layer.** Mechanism-level hypotheses that no drug–target database can reach (for example,
#    raising the abundance of a protein that has no activator). Every citation is resolved live against Europe PMC
#    and printed with its title; a citation that does not resolve is flagged rather than used.
#
# The same code is then run unchanged for the other genes of the disorder's virtual panel, to show it generalises.

# %%
OT_API = "https://api.platform.opentargets.org/api/v4/graphql"
# Sign map: +1 = the node's activity protects in a loss-of-function / hypomorphic state (inhibiting it is harmful);
# -1 = the node's activity lowers the causal protein (inhibiting it could help); 0 = indeterminate. Each entry cites
# its source. Nodes not listed are treated as indeterminate. The causal gene itself is always +1.
SIGN_MAP = {
    "BUB1B": {
        "SIRT2": (+1, "SIRT2 deacetylates BUBR1-K668 and raises BUBR1 abundance; SIRT2 overexpression extends lifespan of BubR1-hypomorphic mice", ["24825348"]),
        "CREBBP": (-1, "CBP acetylates BUBR1-K668, counteracting SIRT2 and promoting BUBR1 loss", ["24825348"]),
        "AURKB": (+1, "Aurora B phosphorylation of BUBR1-S39 is required for protective BUBR1-K250 acetylation", ["41533580"]),
        "KAT2B": (+1, "PCAF acetylates BUBR1-K250, which blocks its APC/C-dependent degradation", ["19407811"]),
        "CDC20": (-1, "APC/C-CDC20 degrades BUBR1 in mitosis", ["19407811"]),
        "PPP2R1A": (+1, "BUBR1's KARD recruits PP2A-B56 to kinetochores and its C-terminal pseudokinase domain promotes that recruitment, so "
                        "PP2A-B56 activity is part of the function the child's missense allele may impair", ["33207204", "27350047", "23079597"]),
        "HDAC1": (0, "opposite effects at K250 (must stay acetylated) and K668 (must stay deacetylated): net direction unknown", ["19407811", "24825348"]),
        "HDAC2": (0, "as HDAC1", ["19407811", "24825348"]),
        "HDAC4": (0, "as HDAC1", ["19407811", "24825348"]),
    },
}
EXTRA_TARGETS = {"BUB1B": ["SIRT2", "CREBBP", "KAT2B"]}          # mechanism nodes from the literature, outside the PPI network
CAUTION_MOA = re.compile(r"tubulin|microtubule|aurora|polo-like|TTK|kinesin|topoisomerase|proteasome|DNA alkylat|DNA cross|"
                         r"DNA polymerase|DNA synthesis|ribonucleotide reductase|thymidylate|DNA intercal", re.I)


def ot(query, variables, timeout=180):
    d = http_json("POST", OT_API, json={"query": query, "variables": variables}, timeout=timeout)
    if d.get("errors"):
        raise RuntimeError(f"Open Targets: {d['errors'][0].get('message', '')[:200]}")
    return d["data"]


def network_drug_pass(gene):
    gid = GENE_POS[gene][4]
    d = ot("""query($id: String!) { meta { apiVersion { x y z } } target(ensemblId: $id) { approvedSymbol isEssential
              pathways { pathwayId pathway }
              interactions(scoreThreshold: 0.7, page: {index: 0, size: 150}) { count rows { targetB { id approvedSymbol } score } } } }""",
           {"id": gid})
    t = d["target"] or {}
    ver = ".".join(str(v) for v in d["meta"]["apiVersion"].values())
    partners = {}
    for r in (t.get("interactions") or {}).get("rows") or []:
        tb = r.get("targetB") or {}
        if tb.get("id", "").startswith("ENSG") and tb.get("approvedSymbol"):
            partners[tb["id"]] = (tb["approvedSymbol"], max(r.get("score") or 0.0, partners.get(tb["id"], ("", 0.0))[1] or 0.0))
    nodes = {gid: (gene, 1.0, "causal gene")}
    for k, (sym, sc) in partners.items():
        nodes.setdefault(k, (sym, sc, "interaction partner"))
    for sym in EXTRA_TARGETS.get(gene, []):
        if sym in GENE_POS:
            nodes.setdefault(GENE_POS[sym][4], (sym, None, "literature mechanism node"))
    rows = []
    ids = list(nodes)
    for i in range(0, len(ids), 50):
        res = ot("""query($ids: [String!]!) { targets(ensemblIds: $ids) { id approvedSymbol
                  drugAndClinicalCandidates { rows { drug { id name maximumClinicalStage
                     mechanismsOfAction { rows { mechanismOfAction actionType targets { approvedSymbol } } } } } } } }""",
                 {"ids": ids[i:i + 50]})
        for tg in res.get("targets") or []:
            sym, score, why = nodes.get(tg["id"], (tg["approvedSymbol"], None, ""))
            for row in (tg.get("drugAndClinicalCandidates") or {}).get("rows") or []:
                dr = row["drug"]
                moas = [m for m in (dr.get("mechanismsOfAction") or {}).get("rows", [])
                        if any(x.get("approvedSymbol") == sym for x in (m.get("targets") or []))] or (dr.get("mechanismsOfAction") or {}).get("rows", [])
                for m in moas[:1]:
                    rows.append(dict(target=sym, node=why, interaction_score=score, drug=dr["name"].title(), chembl=dr["id"],
                                     stage=dr.get("maximumClinicalStage") or "", action=(m.get("actionType") or "").upper(),
                                     mechanism=m.get("mechanismOfAction") or ""))
    df = pd.DataFrame(rows, columns=["target", "node", "interaction_score", "drug", "chembl", "stage", "action", "mechanism"])
    signs = SIGN_MAP.get(gene, {})

    def classify(r):
        if r.stage != "APPROVAL":
            return "not approved"
        if CAUTION_MOA.search(r.mechanism or ""):
            return "CAUTION: antimitotic / genotoxic mechanism"
        s = +1 if r.target == gene else signs.get(r.target, (0, "", []))[0]
        a = +1 if r.action in ("ACTIVATOR", "AGONIST", "POSITIVE ALLOSTERIC MODULATOR", "STABILISER", "OPENER") else (
            -1 if r.action in ("INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR", "DEGRADER", "INVERSE AGONIST") else 0)
        if s == 0 or a == 0:
            return "indeterminate direction"
        return "plausibly beneficial direction" if s * a > 0 else "CAUTION: harmful direction"

    if len(df):
        df["verdict"] = df.apply(classify, axis=1)
    else:
        df["verdict"] = []
    summary = dict(gene=gene, open_targets_api=ver, essential=t.get("isEssential"), n_pathways=len(t.get("pathways") or []),
                   n_network_nodes=len(nodes), n_drug_target_rows=int(len(df)),
                   n_approved=int((df.stage == "APPROVAL").sum()) if len(df) else 0,
                   n_direct_drugs=int((df.target == gene).sum()) if len(df) else 0,
                   verdicts=df[df.stage == "APPROVAL"].drop_duplicates(["drug", "verdict"]).verdict.value_counts().to_dict() if len(df) else {},
                   pathways=[p["pathway"] for p in (t.get("pathways") or [])][:8])
    return df, summary


# Literature-curated hypotheses (citations resolved live below). Grades: A human trial in a relevant population;
# B in vivo in the same gene defect; C human patient cells, or human trials in another disease with the same molecular
# lesion (e.g. a nonsense allele); D mechanistic inference only; X = direction of effect uncertain or plausibly harmful.
# `lever` names the allele or pathway the hypothesis acts on; the allele-level numbers (stop codon, NMD, tolerance of the
# inserted amino acids, peptides) are computed in section 14b and printed with the hypotheses in the report.
CURATED = [
    dict(rank=1, lever="stop allele: read-through",
         hypothesis="Read-through of the nonsense allele with an approved macrolide, tested first in patient-derived cells",
         agents="azithromycin (approved worldwide; decades of long-term paediatric use in cystic fibrosis, primary ciliary dyskinesia and "
                "bronchiectasis; oral, no nephro- or ototoxicity of the aminoglycoside kind, and it accumulates to high intracellular "
                "concentrations - the exposure argument the other read-through agents lack)",
         mechanism="in a reporter assay across six premature stop codons from primary ciliary dyskinesia patients, azithromycin stimulated "
                   "read-through only at the UGA codons - the class this child's stop allele belongs to - and the effect was modest; in mice "
                   "carrying a germline nonsense allele of the tumour-suppressor Apc it read through that allele, restored full-length Apc, "
                   "lowered beta-catenin / cyclin D1 and prevented intestinal tumour formation. That is the closest precedent we identified for "
                   "read-through acting as cancer prevention in a germline nonsense cancer-predisposition genotype, which is the architecture "
                   "of this child's stop allele",
         grade="C / C",
         caveat="the grade is C, not B, because the in vivo evidence is in another gene: this report reserves B for in vivo data in "
                "the same gene defect. In that same reporter study azithromycin was less efficient than aminoglycosides, and at the two non-UGA stop codons it "
                "REDUCED read-through by more than half, so the effect is codon-specific rather than general; the most responsive reporters "
                "had a UGA followed by cytosine, whereas this child's UGA is followed by adenine, a less favourable context; viability of "
                "primary nasal epithelial cells was preserved, but azithromycin was among the compounds that caused a small, statistically "
                "significant fall in ciliary beat frequency (remaining within the normal range), whereas amlexanox and gentamicin caused none; "
                "the Apc work is a mouse tumour-prevention study in another gene, given in drinking water, where the LOWEST dose was the most "
                "effective (a non-monotonic dose-response) and the authors themselves raise a microbiome confound; read-through of BUB1B has "
                "never been tested; macrolides carry QT prolongation, antimicrobial-resistance and (with long exposure) hearing-loss concerns, "
                "so this is a laboratory hypothesis, not a proposal to treat",
         validation="patient-derived fibroblasts or LCLs: allelic ratio on the cDNA amplicon, N- and C-terminal BUBR1 immunoblots, "
                    "allele-resolved peptides, nocodazole arrest, premature chromatid separation, micronuclei",
         pmids=["37276642", "37385919", "34066907", "41015699"]),
    dict(rank=2, lever="stop allele: NMD inhibition + read-through",
         hypothesis="Rescue protein from the nonsense allele with a dual NMD-inhibitor / read-through drug",
         agents="amlexanox (oral tablet long marketed in Japan for asthma / allergic rhinitis - current marketing status not verifiable from "
                "public regulatory records, so availability is not an argument here; the US topical paste was discontinued in 2017; oral "
                "dosing was also tested in a 42-patient randomised trial in adults with type 2 diabetes, which shows oral tolerability only)",
         mechanism="the nonsense allele's mRNA is predicted to be degraded by NMD, as MVA truncating-allele transcripts are absent in "
                   "patient cells; amlexanox stabilises NMD-prone mRNAs and induces read-through, and for a compound-heterozygous "
                   "missense / nonsense patient genotype (another gene, aspartylglucosaminuria) it restored full-length protein from the "
                   "nonsense allele with a synergistic gain in activity - the same architecture as this child's pair. Two products could "
                   "result: full-length BUBR1 carrying a near-cognate amino acid at the stop, and - where only NMD is relieved - the "
                   "N-terminal fragment. A truncated BUBR1 of 715 BUBR1 residues (731X) still supported the checkpoint and stable "
                   "attachments, whereas another truncation (753X) was expressed only weakly, which the authors note need not mean "
                   "reduced stability, so whether this fragment is stable, functional, or instead "
                   "sequesters CDC20 / BUB3 / PP2A-B56 as a dominant negative is an open, testable question",
         grade="C / D",
         caveat="that rescue was measured in a heterologous expression system, not in patient primary cells; whether the concentrations that "
                "worked in cells are reachable in a person is the gating question and we leave it open rather than guess - amlexanox acted at "
                "5-25 micromolar in the studies that reported the lowest effective concentrations, and at 25-125 micromolar in "
                "others, so roughly 1.5-37 microgram/mL for its molecular weight of 298.3, and we found no peer-reviewed human "
                "oral pharmacokinetic study reporting a plasma Cmax to compare it against (the only quantitative pharmacokinetic data we could "
                "verify are from a dog study), so no fold-difference is claimed here and human exposure data would have to be obtained before "
                "any systemic route could be argued at all; it was also less efficient than aminoglycosides in a head-to-head "
                "reporter comparison; it inhibits TBK1 / IKK-epsilon and TBK1 has mitotic functions, so the net effect on segregation must be measured; "
                "and sustained NMD inhibition itself raises R-loops and DNA damage, which is an unattractive trade in a chromosomal-instability "
                "cancer-predisposition syndrome",
         validation="patient-derived fibroblasts or LCLs: allelic ratio on the cDNA amplicon, N- and C-terminal BUBR1 immunoblots, "
                    "allele-resolved peptides, nocodazole arrest, premature chromatid separation, micronuclei",
         pmids=["22938201", "29247835", "37513945", "20516114", "27702906", "26656453", "28683283", "42788676", "34066907", "41954618"]),
    dict(rank=3, lever="missense allele / pathway: BUBR1 stability",
         hypothesis="NAD+ augmentation to raise residual BUBR1 via SIRT2",
         agents="niacin (nicotinic acid; raises NAD+ by the Preiss-Handler route and does not product-inhibit sirtuins); nicotinamide as the "
                "comparator arm (long-term exposure is documented at 1.2 g/m2/day for five years in a diabetes-prevention trial whose "
                "participants were relatives aged 3-40, so it included children and growth was unaffected); nicotinamide riboside / NMN are "
                "supplements and outside scope",
         mechanism="MVA missense alleles in or near the kinase domain lower BUBR1 through faster protein turnover; SIRT2 keeps "
                   "BUBR1-K668 deacetylated and stable; NMN raised BubR1 in vivo and SIRT2 extended median lifespan of "
                   "BubR1-hypomorphic mice. Niacin corrected systemic NAD+ deficiency and improved muscle performance in adults with "
                   "mitochondrial myopathy, which is why it is the probe of choice here",
         grade="B / X",
         caveat="the in vivo support is a genetic manipulation (SIRT2 overexpression), not the drug - the same limitation for which "
                "senolytics are set aside below, and the reason this is graded B/X rather than higher; nicotinamide is itself a sirtuin product-inhibitor, so the net effect on SIRT2 activity is unproven and niacin is the "
                "cleaner probe; correcting an NAD+ deficiency does not guarantee clinical benefit - in a child with mitochondrial myopathy "
                "and NAD+ deficiency, niacin produced a metabolic response but did not halt the disease, which is the honest ceiling on this "
                "kind of hypothesis; the cancer data are not reassuring so much as absent (nicotinamide reduced keratinocyte cancers in one "
                "randomised trial but not in a second, and an NAD+ precursor suppressed liver-cancer growth in mice; none of these studies "
                "involved chromosomal instability or rhabdomyosarcoma), so an effect on tumour metabolism cannot be excluded and must be "
                "measured in the child's own cells before anything else",
         validation="patient-derived cells: BUBR1 protein and half-life (cycloheximide chase), K668 acetylation, mitotic timing, micronuclei, "
                    "and proliferation of the aneuploid fraction, with and without niacin / NAM",
         pmids=["20516114", "24825348", "39490401", "26488693", "36856616", "12297502", "42449929", "32386566", "39173541", "15043959", "36986177"]),
    dict(rank=4, lever="stop allele: laboratory benchmark only (an assay control, not a candidate)",
         hypothesis="Aminoglycoside read-through as the ex-vivo benchmark - NOT proposed for use in this child",
         agents="gentamicin (laboratory positive control in cultured patient-derived cells only)",
         mechanism="UGA is the most read-through-prone stop in reporter assays and read-through inserts Trp, Cys or Arg there; "
                   "aminoglycosides are the most potent read-through agents available, so gentamicin establishes whether this particular "
                   "allele is read-through-competent at all and calibrates how much protein the approved agents above could ever deliver",
         grade="C / X",
         caveat="systemic use is excluded here, not merely cautioned: gentamicin is nephrotoxic (acute kidney injury in 20-33 % of exposed "
                "children) and irreversibly ototoxic, the i.v. read-through trials excluded pre-existing renal impairment, and this child has "
                "nephrocalcinosis; in Duchenne patients neither the stop codon nor the +4 nucleotide predicted response and one read-through "
                "epitope was immunogenic - response must be measured, not inferred",
         validation="same cultured-cell assays as #1 and #2, as the positive control that shows whether this allele is read-through-competent",
         pmids=["10939566", "22479203", "27702906", "35234826", "38366625", "28691931", "20517938", "27848094", "39098506"]),
    dict(rank=5, lever="symptom: nephrocalcinosis",
         hypothesis="Cause-directed treatment of nephrocalcinosis after a tubular work-up",
         agents="thiazide (hypercalciuria); potassium citrate / alkali (distal RTA, hypocitraturia)",
         mechanism="independent of BUBR1: prematurity and renal tubulopathies explain most paediatric nephrocalcinosis; part of the "
                   "failure to thrive may be renal and correctable; the closest paediatric interventional evidence is indirect - in 11 "
                   "children with X-linked hypophosphataemia whose nephrocalcinosis was caused by their phosphate and vitamin D treatment, "
                   "a thiazide was followed by no radiological progression, in an uncontrolled before-and-after comparison",
         grade="A (indication) / D (MVA-specific)",
         caveat="drug choice depends on urine calcium, citrate, oxalate, pH and serum bicarbonate; this is a work-up recommendation for the "
                "treating nephrologist, not a drug proposal, and prescribing blind could harm",
         validation="standard nephrology work-up; renal function follow-up", pmids=["41615714", "41867929", "41311422", "11433085"]),
]
# Considered and set aside - kept in the report because a reviewer will ask about each of them.
SET_ASIDE = [
    dict(item="Ataluren (read-through)", reason="EU conditional authorisation not renewed (European Commission decision, 28 March "
         "2025) and never approved in the US; a UK licence is reported to have been retained for its Duchenne indication, so it is a "
         "geographically restricted comparator rather than a candidate. It does act on the mechanism of interest - it promotes insertion of "
         "near-cognate tRNAs at premature stops - but availability, not mechanism, is what rules it out here", pmids=["27702906"]),
    dict(item="Hydroxychloroquine (aneuploidy-selective chemoprevention)",
         reason="the screen compound was chloroquine, not hydroxychloroquine; in mosaic variegated aneuploidy the aneuploid cells are "
         "largely normal tissue, so selecting against them could worsen growth failure, and the aneuploidy-vulnerability field has moved on",
         pmids=["21315436", "28775166", "40054068"]),
    dict(item="Metformin (senescence / inflammation)",
         reason="the aneuploidy-selective energy-stress compound in the original screen was AICAR (not approved), not metformin; and at "
         "millimolar concentrations in vitro, under low glucose, metformin itself induced genotoxic stress (micronuclei, nucleoplasmic "
         "bridges) - a direction worth avoiding in a chromosomal-instability syndrome, even though the same study found it non-genotoxic "
         "in high glucose", pmids=["18516091", "22048312", "21315436", "35568294"]),
    dict(item="ATP-site kinase inhibitors as pharmacological chaperones of the destabilised kinase domain",
         reason="the human BUBR1 kinase domain bound no nucleotide in thermal-shift assays, so there is no pocket for such a ligand",
         pmids=["24107129", "22698286"]),
    dict(item="Antisense skipping of the exon that carries the stop", reason="assessed in the allele-level triage (frame, how structured "
         "the exon's residues are, natural isoforms); an n-of-1 ASO is not an approved drug", pmids=["20516114"]),
    dict(item="Senolytic dasatinib + quercetin", reason="p16-positive cell clearance in BubR1-hypomorphic mice was genetic (INK-ATTAC), "
         "not pharmacological; the dasatinib label warns about bone growth in children (short stature here)",
         pmids=["22048312", "25754370"]),
]
CAUTION_CURATED = [
    dict(item="Growth hormone", reason="second-neoplasm risk (RR 1.77) in a cancer-predisposition syndrome with prior rhabdomyosarcoma; growth failure is cell-autonomous", pmids=["35529328"]),
    dict(item="Vinca alkaloids / taxanes and other antimitotics", reason="BubR1-insufficient cells slip through mitotic arrest; anticipate altered response and report outcomes - NOT a reason to withhold curative therapy", pmids=["14744753", "16411201"]),
    dict(item="Topoisomerase-II poisons and other genotoxic agents", reason="appear automatically as approved drugs on the network (above); genotoxic load in a chromosomal-instability syndrome warrants monitoring for secondary malignancy", pmids=["31053147"]),
    dict(item="Ionising radiation for imaging", reason="MVA surveillance guidance advises avoiding radiation exposure: ultrasound / MRI rather than CT where possible", pmids=["39264246"]),
    dict(item="Systemic aminoglycosides", reason="AKI in 20-33 % of exposed children; read-through trials excluded pre-existing renal impairment", pmids=["27848094", "38366625"]),
    dict(item="Sustained systemic NMD inhibition", reason="blocking nonsense-mediated decay raises R-loops and DNA damage and is being developed "
         "deliberately as a way to kill genomically unstable cells; relieving NMD is the point of the amlexanox hypothesis, so this trade-off must "
         "be measured in cells (micronuclei, DNA-damage markers) before any systemic use is ever contemplated in a chromosomal-instability syndrome",
         pmids=["41954618", "42711693", "39667907"]),
    dict(item="HDAC inhibitors as mechanism-based therapy", reason="opposite-signed acetylation sites on BUBR1 make the net effect indeterminate; mitotic histone deacetylation is required for accurate segregation", pmids=["19407811", "24825348", "35922507"]),
    dict(item="Aurora B inhibitors", reason="would remove the phosphorylation required for protective BUBR1-K250 acetylation", pmids=["41533580"]),
    dict(item="Genotoxic conditioning / HSCT", reason="the one published MVA1 transplant ended in graft rejection and death; MDS with monosomy 7 is a surveillance target", pmids=["31053147"]),
]
# Verified references used in the report's prose (resolved live like the others)
TRACK2_TEXT_PMIDS = ["39264246", "11746029", "15108198", "42434306", "40028477", "40473777", "27605262", "25741868", "41257875",
                     "32720330", "16411201", "15475955", "18932004", "15159543", "22698286", "23079597", "15208629", "18516091",
                     "20516114", "21552266", "28553959", "34009673", "39169218", "31738183", "33207204", "42824587", "34066907"]


def resolve_pmids(pmids):
    out = {}
    for p in sorted(set(pmids)):
        try:
            r = http_json("GET", "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                          params={"query": f"EXT_ID:{p} AND SRC:MED", "format": "json", "resultType": "lite", "pageSize": 1}, timeout=60)
            res = (r.get("resultList") or {}).get("result") or []
            out[p] = (dict(title=res[0].get("title", ""), journal=res[0].get("journalTitle", ""), year=res[0].get("pubYear", ""),
                           authors=(res[0].get("authorString") or "").split(",")[0]) if res else None)
        except Exception:       # noqa: BLE001
            out[p] = None
        time.sleep(0.15)
    return out


with stage("Track 2: drug-repurposing evidence"):
    if not CFG["RUN_TRACK2"]:
        raise RuntimeError("disabled in CFG")
    g = HYPS[0]["gene"]
    T2_DF, T2_SUM = network_drug_pass(g)
    S["track2"] = dict(network=T2_SUM)
    log(f"Track 2 network pass for {g}: {T2_SUM}")
    appr = T2_DF[T2_DF.stage == "APPROVAL"].drop_duplicates(["drug", "target"]).sort_values(["verdict", "target", "drug"])
    appr.to_csv(OUT / "track2_network_drugs.csv", index=False)
    md(f"**Approved drugs acting on the {g} network ({T2_SUM['n_network_nodes']} nodes; Open Targets {T2_SUM['open_targets_api']})**\n\n"
       + df_md(appr[["target", "node", "drug", "action", "mechanism", "verdict"]], max_rows=40))
    md("**Sign map used for direction (the only hand-curated input)**\n\n" + df_md(pd.DataFrame(
        [dict(node=k, sign=v[0], rationale=v[1], pmids=", ".join(v[2])) for k, v in SIGN_MAP.get(g, {}).items()])))
    # generality: identical code for the other genes of the stated indication's panel
    gen = []
    for g2 in [x for x in S["pheno"].get("dx_panel", []) if x != g and x in GENE_POS][:6]:
        try:
            _, s2 = network_drug_pass(g2)
            gen.append(dict(gene=g2, network_nodes=s2["n_network_nodes"], drugs_on_gene=s2["n_direct_drugs"],
                            approved_rows=s2["n_approved"], **{k: v for k, v in s2["verdicts"].items()}))
        except Exception as e:      # noqa: BLE001
            gen.append(dict(gene=g2, error=str(e)[:80]))
    S["track2"]["generality"] = gen
    if gen:
        md("**Same automated pass for the other genes of the indication panel**\n\n" + df_md(pd.DataFrame(gen).fillna(0)))
    # literature layer, citations resolved live
    all_p = [p for h in CURATED for p in h["pmids"]] + [p for c in CAUTION_CURATED + SET_ASIDE for p in c["pmids"]] + \
            [p for v in SIGN_MAP.get(g, {}).values() for p in v[2]] + TRACK2_TEXT_PMIDS
    PMIDS = resolve_pmids(all_p)
    bad = [p for p, v in PMIDS.items() if not v]
    S["track2"].update(curated=CURATED if g == "BUB1B" else [], caution=CAUTION_CURATED if g == "BUB1B" else [],
                       set_aside=SET_ASIDE if g == "BUB1B" else [],
                       pmids={p: v for p, v in PMIDS.items() if v}, unresolved_pmids=bad)
    log(f"citations resolved live: {len(PMIDS) - len(bad)} of {len(PMIDS)}" + (f"; unresolved (not used): {bad}" if bad else ""))
    if g == "BUB1B":
        md("**Ranked hypotheses (literature layer)**\n\n" + df_md(pd.DataFrame(CURATED)[["rank", "lever", "hypothesis", "agents", "grade", "caveat"]]))
        md("**Considered and set aside**\n\n" + df_md(pd.DataFrame(SET_ASIDE)[["item", "reason"]]))
        md("**Caution list for the clinical team**\n\n" + df_md(pd.DataFrame(CAUTION_CURATED)[["item", "reason"]]))
    else:
        md(f"The literature-curated layer is specific to BUBR1 biology and is not applied to {g}.")

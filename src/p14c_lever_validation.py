# %% [markdown]
# ## 14c · Positive control: does the triage classifier recover published answers?
#
# Stage 14b *asserts* that each allele maps to a particular lever. An assertion is not a method. This stage is the
# positive control. A hand-curated table of **published premature-stop alleles in genes outside this child's panel** —
# alleles for which the read-through literature already reports what happens — is pushed through the *same* functions
# that triage the child's own alleles (`cds_model`, `coding_change`, `classify_clinvar`), and the predictions are
# compared with the published answer.
#
# * **What is held out.** Only the gene symbol and the HGVS coding change go in. The stop codon is recomputed from the
#   GRCh38 reference and the MANE Select gene model; the published stop codon is the answer key and never an input.
# * **Why the stop codon is the right thing to score.** It is the single quantity the whole read-through branch turns
#   on: `STOP_ORDER` ranks UGA above UAG and UAA, `NEAR_COGNATE` decides which amino acids read-through would insert,
#   and the +4 nucleotide modifies both. Get the codon wrong and every downstream claim about that allele is wrong.
# * **Provenance.** Every PMID below was checked against Europe PMC (`EXT_ID:<pmid> AND SRC:MED`) while the table was
#   assembled, and the allele is named in the cited paper's title or abstract. Candidates whose allele or citation
#   could not be pinned down are listed in `unverified_dropped` rather than quietly discarded.
# * **Public data only.** The child's VCF, genotypes and phenotype are not read here. The inputs are the public
#   reference genome, the public MANE gene model and the literature table in the next cell.
#
# The agreement figure is deliberately a *necessary* condition, not a sufficient one: it shows the classifier reads
# the genome correctly, not that read-through will work in this child. The honest limits are reported with it.

# %%
# Published premature-stop alleles with read-through literature, in genes that are not in this child's panel.
# `expected_stop` is the answer key: it follows from the published coding change and the genetic code, and is never
# read from the reference. `hgvsc` is written on the MANE Select transcript the pipeline uses, which for MECP2 is not
# the transcript the original report numbered against (see `note`).
PUBLISHED_PTC = [
    # --- cystic fibrosis (CFTR, NM_000492.4 / ENST00000003084; the legacy CF numbering is the MANE numbering) --------
    dict(gene="CFTR", hgvsc="c.178G>T", hgvsp="p.Glu60Ter", expected_stop="TAG", pmids=["36143233"],
         readthrough="achieved in cells (G418 or ELX-02 with CFTR modulators); an N-terminal allele that escapes NMD"),
    dict(gene="CFTR", hgvsc="c.263T>G", hgvsp="p.Leu88Ter", expected_stop="TGA", pmids=["36143233"],
         readthrough="achieved in cells (G418 or ELX-02 with CFTR modulators); an N-terminal allele that escapes NMD"),
    dict(gene="CFTR", hgvsc="c.366T>A", hgvsp="p.Tyr122Ter", expected_stop="TAA", pmids=["17394637", "28376702", "35487895"],
         readthrough="achieved; gentamicin read-through was at least five-fold that of the three UGA alleles tested "
                     "alongside it, and this allele does not need UPF2/UPF3 for its decay"),
    dict(gene="CFTR", hgvsc="c.1624G>T", hgvsp="p.Gly542Ter", expected_stop="TGA",
         pmids=["18272502", "24251786", "38725656", "17394637"],
         readthrough="achieved in a CF mouse model, in patient bronchial cells and by mass spectrometry of the product"),
    dict(gene="CFTR", hgvsc="c.1657C>T", hgvsp="p.Arg553Ter", expected_stop="TGA", pmids=["24251786"],
         readthrough="achieved, but the weakest of the four common CFTR stop alleles tested with NB124"),
    dict(gene="CFTR", hgvsc="c.3484C>T", hgvsp="p.Arg1162Ter", expected_stop="TGA",
         pmids=["24251786", "38725656", "17394637"],
         readthrough="achieved; among the most responsive CFTR alleles to synthetic aminoglycosides"),
    dict(gene="CFTR", hgvsc="c.3587C>G", hgvsp="p.Ser1196Ter", expected_stop="TGA", pmids=["38725656"],
         readthrough="achieved; Cys, Arg and Trp incorporation at the UGA confirmed by tandem mass spectrometry"),
    dict(gene="CFTR", hgvsc="c.3846G>A", hgvsp="p.Trp1282Ter", expected_stop="TGA",
         pmids=["24251786", "17290305", "17394637"],
         readthrough="achieved, but only in patients whose cells leave enough nonsense transcript for NMD to spare"),
    # --- Hurler syndrome (IDUA, NM_000203.5 / ENST00000514224) -------------------------------------------------------
    dict(gene="IDUA", hgvsc="c.208C>T", hgvsp="p.Gln70Ter", expected_stop="TAG", pmids=["11159948"],
         readthrough="achieved in patient fibroblasts (gentamicin), restoring 2.8 % of normal iduronidase activity"),
    dict(gene="IDUA", hgvsc="c.1205G>A", hgvsp="p.Trp402Ter", expected_stop="TAG", pmids=["11159948"],
         readthrough="achieved in patient fibroblasts (gentamicin), restoring 2.8 % of normal iduronidase activity"),
    # --- Rett syndrome (MECP2) ---------------------------------------------------------------------------------------
    dict(gene="MECP2", hgvsc="c.538C>T", hgvsp="p.Arg180Ter", expected_stop="TGA", pmids=["21120444"],
         readthrough="achieved in HeLa cells and in ear fibroblasts of a knock-in mouse (gentamicin; NB54 and NB84 "
                     "more effective)",
         note="published as p.Arg168Ter on NM_004992.4 (MANE Plus Clinical); the MANE Select transcript "
              "NM_001110792.2 numbers the same residue 180, so the coding change is restated on it"),
    # --- recessive dystrophic epidermolysis bullosa (COL7A1, NM_000094.4 / ENST00000681320) --------------------------
    dict(gene="COL7A1", hgvsc="c.751C>T", hgvsp="p.Gln251Ter", expected_stop="TAG", pmids=["25155989"],
         readthrough="achieved in patient keratinocytes (geneticin, gentamicin, paromomycin), restoring full-length "
                     "type VII collagen"),
    dict(gene="COL7A1", hgvsc="c.1732C>T", hgvsp="p.Arg578Ter", expected_stop="TGA", pmids=["25155989", "38366625"],
         readthrough="achieved in patient cells and in patients' skin after intravenous gentamicin"),
    # --- Herlitz junctional epidermolysis bullosa (LAMB3, NM_000228.3 / ENST00000356082) -----------------------------
    dict(gene="LAMB3", hgvsc="c.870T>A", hgvsp="p.Cys290Ter", expected_stop="TGA", pmids=["29946029"],
         readthrough="achieved in laminin-beta-3-null keratinocytes (gentamicin), restoring laminin 332 assembly"),
    dict(gene="LAMB3", hgvsc="c.1903C>T", hgvsp="p.Arg635Ter", expected_stop="TGA", pmids=["29946029"],
         readthrough="achieved in laminin-beta-3-null keratinocytes (gentamicin), restoring laminin 332 assembly"),
    # --- tumour-suppressor nonsense allele (TP53, NM_000546.6 / ENST00000269305) -------------------------------------
    dict(gene="TP53", hgvsc="c.637C>T", hgvsp="p.Arg213Ter", expected_stop="TGA", pmids=["31120902"],
         readthrough="achieved in cancer cells (G418, strongly potentiated by mefloquine), restoring functional p53"),
]

# Candidates assembled and then excluded, kept in the record so that the agreement figure cannot be read as a
# selection of the cases that happened to work.
VALIDATION_DROPPED = [
    dict(gene="DMD", allele="eight premature stop codons of Duchenne / Becker patients, two of them refractory to "
                            "aminoglycosides", pmids=["14991821"],
         reason="the paper reports the read-through level per allele but does not name the alleles, so no coding "
                "change could be attached to a citation"),
    dict(gene="DMD", allele="the mdx UAA stop codon", pmids=["10939566"],
         reason="a mouse allele; the human orthologous position is not the reported one"),
    dict(gene="COL7A1", allele="p.Arg163Ter", pmids=["25155989"],
         reason="c.487C>T does not give an arginine codon on MANE Select NM_000094.4 (codon 163 is CTG), so the "
                "published label and the coding change could not be reconciled"),
    dict(gene="CFTR", allele="the six common CFTR premature stops surveyed for inserted amino acids", pmids=["41015699"],
         reason="the six alleles are not named in the abstract, so they were not added as separate entries"),
]

# The priors the classifier encodes, and the published work they come from. Checked, not tuned: these thresholds are
# fixed before the run and no leaderboard or panel feedback touches them.
VALIDATION_BASIS = [
    dict(prior="UGA is read through more readily than UAG or UAA (`STOP_ORDER`)", pmids=["10939566", "41015699"]),
    dict(prior="a cytosine at +4 raises read-through; the order is C > U > A >= G (`plus4` in the triage)",
         pmids=["10939566", "22479203"]),
    dict(prior="read-through inserts the amino acids of near-cognate tRNAs, set by the stop codon (`NEAR_COGNATE`)",
         pmids=["25733896", "27702906"]),
    dict(prior="nonsense-mediated decay governs how much transcript is left for read-through to act on",
         pmids=["17290305", "35487895"]),
]

_PubRow = collections.namedtuple("_PubRow", "gene chrom pos ref alt mc")        # the fields `classify_clinvar` reads


def _hgvsc_snv(h):
    """(CDS position, reference base, alternate base) of a coding substitution written `c.<n><ref>><alt>`."""
    m = re.fullmatch(r"c\.(\d+)([ACGT])>([ACGT])", h or "")
    return (int(m.group(1)), m.group(2), m.group(3)) if m else None


def _hgvsp_ter(h):
    """(residue, one-letter wild-type amino acid) of a protein change written `p.<Xxx><n>Ter`."""
    inv = {b: a for a, b in AA3.items()}
    m = re.fullmatch(r"p\.([A-Za-z]{3})(\d+)Ter", h or "")
    return (int(m.group(2)), inv.get(m.group(1))) if m else None


def validate_ptc(e, fa, exclude=()):
    """Run one published allele through the production triage and compare with the published answer. Returns a row
    with `computed_stop=None` and a `skipped` reason whenever the allele cannot be placed on the gene model; never
    raises, so one bad row cannot cost the rest of the control."""
    g = e["gene"]
    row = dict(gene=g, hgvsc=e["hgvsc"], hgvsp=e.get("hgvsp", ""), expected_stop=e["expected_stop"],
               expected_stop_rna=e["expected_stop"].replace("T", "U"), computed_stop=None, computed_stop_rna="",
               lever=None, expected_lever=("read-through: UGA stop" if e["expected_stop"] == "TGA"
                                           else "read-through: UAG / UAA stop"),
               transcript="", residue=None, wt_aa="", plus4_rna="", in_frame_exon=None, match=False,
               lever_match=False, readthrough=e.get("readthrough", ""), pmids=list(e.get("pmids", [])),
               note=e.get("note", ""), skipped=None)
    if g in exclude:        # a gene of this child's own panel would stop being an independent control
        row["skipped"] = "gene is in this child's diagnostic panel, so it is not an independent control"
        return row
    if g not in GENE_TX or g not in GENE_POS:
        row["skipped"] = "gene absent from the MANE gene model"
        return row
    pr = _hgvsc_snv(e["hgvsc"])
    if pr is None:
        row["skipped"] = "coding change is not a simple substitution"
        return row
    c, ref_t, alt_t = pr
    m = cds_model(GENE_TX[g][0], fa)
    row["transcript"] = TX[m["tid"]].get("version") or m["tid"]
    if not 1 <= c <= len(m["seq"]):
        row["skipped"] = f"CDS position {c} is beyond the {len(m['seq'])} nt MANE coding sequence"
        return row
    plus = TX[m["tid"]]["strand"] == "+"
    pos = m["gpos"][c - 1]              # HGVS c. alleles are in transcript orientation; the triage wants genomic ones
    ref_g, alt_g = (ref_t, alt_t) if plus else (ref_t.translate(COMP), alt_t.translate(COMP))
    cc = coding_change(m, pos, ref_g, alt_g)
    if cc is None:
        row["skipped"] = "the reference base at this CDS position does not match the published coding change"
        return row
    pp = _hgvsp_ter(e.get("hgvsp", ""))
    if pp and (cc["residue"], cc["wt_aa"]) != pp:       # a second, independent check on the coordinate mapping
        row["skipped"] = (f"computed residue {cc['wt_aa']}{cc['residue']} does not match the published "
                          f"{pp[1]}{pp[0]}")
        return row
    r = _PubRow(gene=g, chrom=GENE_POS[g][0], pos=pos, ref=ref_g, alt=alt_g, mc="SO:0001587|nonsense")
    lever, stop, in_frame = classify_clinvar(r, m)      # the production classifier, unchanged
    row.update(computed_stop=stop, computed_stop_rna=(stop or "").replace("T", "U"), lever=lever,
               residue=int(cc["residue"]), wt_aa=cc["wt_aa"], plus4_rna=cc["plus4"].replace("T", "U"),
               in_frame_exon=bool(in_frame), match=bool(stop == e["expected_stop"]),
               lever_match=bool(lever == row["expected_lever"]))
    return row


with stage("Track 2: triage validation on published alleles"):
    if not CFG["RUN_TRACK2"]:
        raise RuntimeError("disabled in CFG")
    from pyfaidx import Fasta
    fav = Fasta(RES["ref"], as_raw=True, sequence_always_upper=True)
    own = {x for x in ((S.get("pheno") or {}).get("dx_panel") or []) if x} | {(S.get("levers") or {}).get("gene")}
    VROWS = [validate_ptc(e, fav, exclude=own) for e in PUBLISHED_PTC]      # plain dicts: already JSON-serialisable
    V = pd.DataFrame(VROWS)
    done = V[V.skipped.isna()]
    n_tested, n_match = int(len(done)), int(done.match.sum())
    n_lever = int(done.lever_match.sum())

    # ---- agreement, per gene and as a lever-vs-published-codon table ------------------------------------------------
    by_gene = [dict(gene=g, n_tested=int((d.skipped.isna()).sum()), n_codon_match=int(d.match.sum()),
                    n_skipped=int(d.skipped.notna().sum()),
                    published_codons=", ".join(sorted(set(d.expected_stop_rna))))
               for g, d in V.groupby("gene", sort=True)]
    conf = pd.crosstab(done.lever, done.expected_stop_rna) if n_tested else pd.DataFrame()
    if n_tested:
        conf = conf.reindex(index=[c for c in LEVER_ORDER if c in conf.index],
                            columns=[c for c in ("UGA", "UAG", "UAA") if c in conf.columns], fill_value=0)
    S["lever_validation"] = dict(
        n_entries=int(len(PUBLISHED_PTC)), n_tested=n_tested, n_codon_match=n_match,
        frac_codon_match=round(n_match / n_tested, 3) if n_tested else 0.0,
        n_lever_match=n_lever, n_skipped=int(len(V) - n_tested),
        genes=sorted(set(V.gene)), n_pmids=len({p for e in PUBLISHED_PTC for p in e["pmids"]}),
        n_pmids_cited=len({p for d in PUBLISHED_PTC + VALIDATION_DROPPED + VALIDATION_BASIS for p in d["pmids"]}),
        published_codon_counts=done.expected_stop_rna.value_counts().to_dict() if n_tested else {},
        by_gene=by_gene,
        confusion=conf.reset_index().to_dict("records") if n_tested else [],
        entries=VROWS,
        unverified_dropped=VALIDATION_DROPPED, basis=VALIDATION_BASIS,
        pmid_source="Europe PMC REST (EXT_ID:<pmid> AND SRC:MED), checked when the table was assembled")
    VOUT = V.assign(pmids=V.pmids.map(lambda p: "; ".join(p)), skipped=V.skipped.fillna(""))
    VOUT.to_csv(OUT / "lever_validation.csv", index=False)    # reusable product: the control set and its predictions

    md(f"**Positive control: {n_tested} published premature-stop alleles in {done.gene.nunique() if n_tested else 0} "
       f"genes outside this child's panel**\n\nOnly the gene symbol and the HGVS coding change were supplied; the stop "
       f"codon was recomputed from GRCh38 and the MANE Select model by the same `cds_model` / `coding_change` / "
       f"`classify_clinvar` path that triages the child's alleles.\n\n"
       + df_md(VOUT[["gene", "hgvsc", "hgvsp", "transcript", "residue", "expected_stop_rna", "computed_stop_rna",
                     "plus4_rna", "lever", "match", "pmids", "skipped"]], max_rows=len(VOUT)))
    if n_tested:
        md(f"**Computed lever against the published stop codon** (agreement "
           f"{n_match}/{n_tested} on codon identity, {n_lever}/{n_tested} on lever class)\n\n"
           + df_md(conf.reset_index()))
    md("**What this does and does not establish**\n\n"
       "The control shows that the triage places a published coding change on the right codon of the right "
       "transcript and assigns the lever its own rules imply - the necessary condition for every read-through claim "
       "downstream. It does not show that read-through will work in this child. Two honest limits: every allele that "
       "could be tied to a named citation is read-through-positive in the literature, with graded rather than binary "
       "efficiency, so the set cannot measure specificity against true negatives; and CFTR p.Tyr122Ter is a UAA "
       "allele whose reported gentamicin read-through exceeded that of three UGA alleles tested beside it, which is "
       "why `STOP_ORDER` is treated as a prior to be measured and not a rule - the +4 cytosine at that allele is the "
       "likely reason, and the triage reports the +4 nucleotide for exactly this purpose.\n\n"
       + df_md(pd.DataFrame([dict(prior=b["prior"], pmids="; ".join(b["pmids"])) for b in VALIDATION_BASIS]))
       + "\n\nSaved as `lever_validation.csv` (public alleles and public "
       "literature only). Candidates excluded while the table was assembled: "
       + "; ".join(f"{d['gene']} {d['allele']} ({d['reason']})" for d in VALIDATION_DROPPED) + ".")
    log(f"triage validation: {n_match}/{n_tested} published stop codons recovered "
        f"({S['lever_validation']['frac_codon_match']:.0%}), lever class {n_lever}/{n_tested}, "
        f"{S['lever_validation']['n_skipped']} skipped, {len(VALIDATION_DROPPED)} candidates dropped unverified, "
        f"codons tested: {json.dumps(S['lever_validation']['published_codon_counts'])}")

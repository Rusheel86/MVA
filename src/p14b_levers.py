# %% [markdown]
# ## 14b · Allele-level therapeutic triage (gene-agnostic)
#
# Drug–target databases answer "which drug acts on a protein near this gene?". For a loss-of-function disorder the more
# useful question is allele-specific: **what can still be done with each of the child's two alleles?** This stage answers
# it from the reference sequence, the MANE gene model, ESM-2 and the AlphaFold model, for whichever gene leads:
#
# * **Premature stop codons** — stop-codon identity and the following (+4) nucleotide, the main sequence determinants of
#   drug-induced read-through; whether nonsense-mediated decay removes the transcript; which amino acids read-through would
#   insert (near-cognate tRNAs) and how ESM-2 rates them at that residue; whether the stop lies in an exon whose skipping
#   keeps the reading frame (antisense exon skipping), how structured that exon's residues are, and whether a natural
#   isoform already skips it.
# * **Missense alleles** — domain, AlphaFold confidence, burial, distance to literature-reported patient alleles.
# * **The pair** — the distance between the alleles in the mRNA (one RT-PCR amplicon from patient RNA phases them
#   directly) and allele-specific tryptic peptides that tell apart the protein made from each allele, so a read-through
#   product can be measured in the child's own cells.
# * **The path to a firm diagnosis** — which evidence would move the weaker allele from VUS to (likely) pathogenic.
# * **Scalability** — the same triage over every ClinVar P/LP allele in the genes of the disorder's panel.

# %%
STOP_CODONS = {"TAA", "TAG", "TGA"}
STOP_ORDER = {"TGA": 1, "TAG": 2, "TAA": 3}           # read-through permissiveness, most to least permissive
NEAR_COGNATE = {"TGA": "WCR", "TAG": "QYK", "TAA": "QYK"}     # amino acids inserted at each stop during induced read-through
COMP = str.maketrans("ACGTN", "TGCAN")
CODON = {a + b + c: aa for (a, b, c), aa in
         zip(itertools.product("TCAG", repeat=3), "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
MONO = dict(G=57.02146, A=71.03711, S=87.03203, P=97.05276, V=99.06841, T=101.04768, C=103.00919, L=113.08406, I=113.08406,
            N=114.04293, D=115.02694, Q=128.05858, K=128.09496, E=129.04259, M=131.04049, H=137.05891, F=147.06841,
            R=156.10111, Y=163.06333, W=186.07931)
LEVER_ORDER = ["read-through: UGA stop", "read-through: UAG / UAA stop", "exon skipping: lesion in an in-frame exon",
               "truncated protein escapes NMD", "missense: protein-level", "splice: splice-switching candidate",
               "start-loss: translation must restart downstream", "no allele-specific lever"]


def cds_model(tid, fa):
    """Coding sequence (5'→3'), the genomic position of every coding base, and its CDS segment index."""
    t = TX[tid]
    seq, gpos, seg = [], [], []
    for i, (s, e) in enumerate(t["cds"]):
        sub = fa[t["chrom"]][s - 1:e]
        if t["strand"] == "+":
            seq.append(sub)
            gpos.extend(range(s, e + 1))
        else:
            seq.append(sub[::-1].translate(COMP))
            gpos.extend(range(e, s - 1, -1))
        seg.extend([i] * (e - s + 1))
    return dict(tid=tid, seq="".join(seq), gpos=gpos, seg=seg, idx={p: i for i, p in enumerate(gpos)})


def coding_change(m, pos, ref, alt):
    """Codon-level consequence of a coding SNV (None when it is not a coding SNV or the reference does not match)."""
    t = TX[m["tid"]]
    if len(ref) != 1 or len(alt) != 1 or pos not in m["idx"]:
        return None
    i = m["idx"][pos]
    plus = t["strand"] == "+"
    if m["seq"][i] != (ref if plus else ref.translate(COMP)):
        return None
    k = i // 3
    wt = m["seq"][3 * k:3 * k + 3]
    mut = wt[:i % 3] + (alt if plus else alt.translate(COMP)) + wt[i % 3 + 1:]
    nxt = m["seq"][3 * k + 3] if 3 * k + 3 < len(m["seq"]) else ""
    return dict(c=i + 1, residue=k + 1, wt_codon=wt, mut_codon=mut, wt_aa=CODON.get(wt, "?"), mut_aa=CODON.get(mut, "?"),
                plus4=nxt, segment=m["seg"][i])


def exon_info(m, i):
    """What skipping CDS segment i would do: frame, residues removed, and whether it is a whole internal exon."""
    t = TX[m["tid"]]
    s, e = t["cds"][i]
    n = len(t["cds"])
    length = e - s + 1
    first = sum(b - a + 1 for a, b in t["cds"][:i])
    exon_no = next((k + 1 for k, (xs, xe) in enumerate(t["exons"]) if xs <= s and e <= xe), None)
    whole = (s, e) in set(map(tuple, t["exons"]))
    return dict(exon=exon_no, n_exons=len(t["exons"]), length_nt=length, in_frame=length % 3 == 0,
                skippable=bool(whole and 0 < i < n - 1), aa_from=first // 3 + 1, aa_to=(first + length - 1) // 3 + 1,
                genomic=(s, e), flank_up=t["cds"][i - 1] if i > 0 else None, flank_dn=t["cds"][i + 1] if i < n - 1 else None)


def natural_skipping(gene, ex):
    """Annotated protein-coding Ensembl transcripts that skip this exon but keep both flanking exons."""
    if not (ex["flank_up"] and ex["flank_dn"]):
        return []
    d = http_json("GET", f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}", params={"expand": 1},
                  headers={"Content-Type": "application/json"}, timeout=60)
    hits = []
    for tr in d.get("Transcript", []):
        if tr.get("biotype") != "protein_coding":
            continue
        exs = [(x["start"], x["end"]) for x in tr.get("Exon", [])]

        def ov(a, b):
            return any(xs <= b and xe >= a for xs, xe in exs)
        if not ov(*ex["genomic"]) and ov(*ex["flank_up"]) and ov(*ex["flank_dn"]):
            hits.append(tr["id"])
    return hits


def esm_row(gene, residue):
    """ESM-2 log-likelihood ratios (AA20 order) for every substitution at one residue, and the protein-wide values
    (None when only single positions were scored on CPU; a missing position is then scored on CPU on demand)."""
    ctx = globals().get("ESM_CTX")
    if not ctx or ctx["gene"] != gene:
        return None, None
    k = residue - 1 - ctx["off"]
    if ctx.get("llr") is not None:
        if not 0 <= k < len(ctx["llr"]):
            return None, None
        flat = ctx["llr"].ravel()
        return ctx["llr"][k], flat[flat != 0]
    rows = ctx.setdefault("rows", {})
    if residue not in rows and 0 <= k < len(ctx.get("wseq", "")):
        try:
            rows[residue] = esm_single(ctx["wseq"], k)
        except Exception as e:      # noqa: BLE001
            log(f"ESM-2 single-position scoring of residue {residue} skipped ({str(e)[:100]})")
    return rows.get(residue), None


def af_ctx(gene):
    """The AlphaFold/UniProt context, only when it belongs to this gene and uses the MANE residue numbering."""
    ctx = globals().get("AF_CTX")
    return ctx if ctx and ctx["gene"] == gene and ctx.get("numbering_matches_mane") else None


def af_at(gene, residues):
    ctx = af_ctx(gene)
    if not ctx:
        return None
    vals = [ctx["ca"][r][1] for r in residues if r in ctx["ca"]]
    return round(float(np.mean(vals)), 1) if vals else None


def features_overlapping(gene, a, b):
    ctx = af_ctx(gene)
    if not ctx:
        return []
    return [f"{f['name']} ({f['start']}-{f['end']})" for f in ctx["feats"] if f["start"] and f["end"] and f["start"] <= b and f["end"] >= a]


PROTEASES = {        # cleavage rules (C-terminal side unless noted); the first that yields observable peptides is used
    "trypsin": lambda s, i: s[i] in "KR" and s[i + 1] != "P",
    "Glu-C": lambda s, i: s[i] == "E",
    "chymotrypsin": lambda s, i: s[i] in "FWY" and s[i + 1] != "P",
    "Asp-N": lambda s, i: s[i + 1] == "D",          # cleaves N-terminal to Asp
    "Lys-C": lambda s, i: s[i] == "K",
}


def digest_peptide(seq, residue, protease="trypsin"):
    """(peptide, first residue, last residue) containing a 1-based residue after complete digestion with one protease."""
    rule = PROTEASES[protease]
    cuts = [0] + [i + 1 for i in range(len(seq) - 1) if rule(seq, i)] + [len(seq)]
    for a, b in zip(cuts, cuts[1:]):
        if a < residue <= b:
            return seq[a:b], a + 1, b
    return "", 0, 0


def observable(p):
    return 7 <= len(p) <= 30


def site_peptides(forms, residue, avoid=None, min_delta_da=0.02):
    """forms: {label: (protein sequence, made_by)}, each the protein actually made by that allele / product. Returns the first
    protease whose peptides covering the residue are MS-observable, do not also span the other allele's site (`avoid`), and
    differ pairwise in monoisotopic mass by at least `min_delta_da`."""
    for prot_ in PROTEASES:
        peps = {k: digest_peptide(sq, residue, prot_) for k, (sq, _) in forms.items()}
        masses = sorted(peptide_mass(x[0]) for x in peps.values())
        if (all(observable(x[0]) for x in peps.values())
                and not any(avoid and lo <= avoid <= hi for _, lo, hi in peps.values())
                and all(b - a >= min_delta_da for a, b in zip(masses, masses[1:]))):
            return [dict(protease=prot_, form=k, peptide=peps[k][0], mass=peptide_mass(peps[k][0]), made_by=forms[k][1]) for k in forms]
    return []


def peptide_mass(p):
    return round(sum(MONO.get(a, 0.0) for a in p) + 18.01056, 3)


def tx_model(tid, fa):
    """Spliced transcript (exons 5'→3'), the genomic position and exon number of every base, and the first coding base."""
    t = TX[tid]
    seq, gpos, exn = [], [], []
    for k, (s, e) in enumerate(t["exons"]):
        sub = fa[t["chrom"]][s - 1:e]
        if t["strand"] == "+":
            seq.append(sub)
            gpos.extend(range(s, e + 1))
        else:
            seq.append(sub[::-1].translate(COMP))
            gpos.extend(range(e, s - 1, -1))
        exn.extend([k + 1] * (e - s + 1))
    first = t["cds"][0][0] if t["strand"] == "+" else t["cds"][0][1]
    return dict(seq="".join(seq), gpos=gpos, exon=exn, cds0=gpos.index(first), chrom=t["chrom"])


def positions_with_variants(chrom, positions, source):
    """How many of the given genomic positions carry a variant: in the child's VCF ('child', count only — no genotypes
    leave this function) or a public variant with minor-allele frequency >= 1 % (Ensembl REST, reference coordinates)."""
    if not positions:
        return 0
    a, b = min(positions), max(positions)
    if source == "child":
        raw = chrom if S.get("vcf_chr_style", True) else chrom[3:]
        txt = sh(f"bcftools query -r {raw}:{a}-{b} -f '%POS\\n' '{S['vcf_raw']}'", name="primer_site_child", capture=True)
        return len({int(x) for x in txt.split()} & set(positions))
    d = http_json("GET", f"https://rest.ensembl.org/overlap/region/human/{chrom.replace('chr', '')}:{a}-{b}",
                  params={"feature": "variation"}, headers={"Content-Type": "application/json"}, timeout=60)
    return len({p for x in d if (x.get("minor_allele_freq") or 0) >= 0.01
                for p in range(int(x["start"]), int(x["end"]) + 1)} & set(positions))


def design_rtpcr_primers(m, c1, c2, fa, min_gap=60, flank=450, n=3):
    """Primer3 primers for one RT-PCR amplicon spanning both coding positions (c1 < c2), each primer >= `min_gap` nt from
    the nearer allele so that both sites are read by Sanger or amplicon sequencing. Pairs are ranked by: no variant of the
    child under a primer (allele dropout would distort the allelic ratio), no common public SNP, a primer spanning an
    exon-exon junction (no product from genomic DNA), then the Primer3 penalty."""
    ensure_modules({"primer3": "primer3-py"})
    import primer3
    tm = tx_model(m["tid"], fa)
    s = tm["seq"]
    if s[tm["cds0"]:tm["cds0"] + len(m["seq"])] != m["seq"]:
        raise RuntimeError("exon and CDS models disagree")
    i1, i2 = tm["cds0"] + c1 - 1, tm["cds0"] + c2 - 1
    lo, hi = max(0, i1 - flank), min(len(s), i2 + flank + 1)
    design = getattr(primer3.bindings, "design_primers", None) or primer3.bindings.designPrimers
    span = i2 - i1 + 2 * min_gap + 1
    r = design(dict(SEQUENCE_ID="amplicon", SEQUENCE_TEMPLATE=s[lo:hi], SEQUENCE_TARGET=[i1 - lo - min_gap, span]),
               dict(PRIMER_OPT_SIZE=22, PRIMER_MIN_SIZE=18, PRIMER_MAX_SIZE=27, PRIMER_OPT_TM=60.0, PRIMER_MIN_TM=57.0,
                    PRIMER_MAX_TM=63.0, PRIMER_MIN_GC=40.0, PRIMER_MAX_GC=60.0, PRIMER_MAX_POLY_X=4, PRIMER_PAIR_MAX_DIFF_TM=3.0,
                    PRIMER_NUM_RETURN=25, PRIMER_PRODUCT_SIZE_RANGE=[[span, span + 2 * (flank - min_gap)]]))
    pairs = []
    for k in range(r.get("PRIMER_PAIR_NUM_RETURNED", 0)):
        ls, ll = r[f"PRIMER_LEFT_{k}"]
        rs, rl = r[f"PRIMER_RIGHT_{k}"]
        left = range(lo + ls, lo + ls + ll)
        right = range(lo + rs - rl + 1, lo + rs + 1)
        ex_l, ex_r = [tm["exon"][j] for j in left], [tm["exon"][j] for j in right]

        def spans(ex):
            return len(set(ex)) > 1 and min(collections.Counter(ex).values()) >= 4      # >= 4 nt on each side of the junction
        gp = [tm["gpos"][j] for j in list(left) + list(right)]
        pairs.append(dict(left=r[f"PRIMER_LEFT_{k}_SEQUENCE"], right=r[f"PRIMER_RIGHT_{k}_SEQUENCE"],
                          tm_left=round(r[f"PRIMER_LEFT_{k}_TM"], 1), tm_right=round(r[f"PRIMER_RIGHT_{k}_TM"], 1),
                          product_nt=r[f"PRIMER_PAIR_{k}_PRODUCT_SIZE"], exons_left="-".join(map(str, sorted(set(ex_l)))),
                          exons_right="-".join(map(str, sorted(set(ex_r)))), junction_spanning=spans(ex_l) or spans(ex_r),
                          genomic_span_bp=abs(tm["gpos"][right[-1]] - tm["gpos"][left[0]]) + 1,
                          site_offsets=(i1 - left[0] + 1, i2 - left[0] + 1), _gpos=gp,
                          penalty=round(r[f"PRIMER_PAIR_{k}_PENALTY"], 2)))
    out = []
    for p in sorted(pairs, key=lambda p: (not p["junction_spanning"], p["penalty"])):
        gp = p.pop("_gpos")
        p["child_variants_under_primers"] = positions_with_variants(tm["chrom"], gp, "child") if S.get("vcf_raw") else None
        try:
            p["common_snps_under_primers"] = positions_with_variants(tm["chrom"], gp, "public")
        except Exception as e:      # noqa: BLE001 - REST unavailable: report as unchecked
            p["common_snps_under_primers"] = None
            log(f"public SNP check of primers skipped ({str(e)[:80]})")
        if not p["child_variants_under_primers"] and not p["common_snps_under_primers"]:
            out.append(p)
        if len(out) == n:
            break
    return out


def reclassification_path(v, partner, pm3_now):
    """Evidence that would move the weaker allele's automated ACMG class, applied cumulatively on the point scale of
    Tavtigian 2020. PM3 follows the pipeline's ClinGen SVI rule: in trans with a Pathogenic partner = moderate (2),
    with a Likely pathogenic partner = supporting (1)."""
    pm3_trans = 2.0 if partner["p_path"] >= 0.99 else 1.0 if partner["p_path"] >= 0.90 else 0.0
    steps = [("current automated evidence (phase unknown" + (f", PM3 {pm3_now:+g}" if pm3_now else "") + ")", pm3_now),
             (f"+ PM3 confirmed in trans (partner allele {'Pathogenic' if pm3_trans == 2 else 'Likely pathogenic' if pm3_trans == 1 else 'below LP'}; "
              "parental testing or the cDNA amplicon)", pm3_trans - pm3_now),
             ("+ PP4: phenotype highly specific for the gene's disorder (ClinGen caveats for heterogeneous disorders apply)", 1),
             ("+ PS3_moderate: allele-specific readout in the child's cells (peptides from this allele's protein; its half-life) "
              "shows reduced protein or function", 2)]
    rows, extra = [], 0.0
    for label, pts in steps:
        extra += pts
        tot = v["points"] + extra
        rows.append(dict(step=label, points=tot, posterior=round(points_to_posterior(tot), 3), classification=points_class(tot)))
    return rows


def classify_clinvar(r, m):
    """Allele-specific lever class of one ClinVar P/LP record."""
    mc = r.mc or ""
    seg = m["seg"][m["idx"][r.pos]] if r.pos in m["idx"] else (m["seg"][m["idx"][r.pos + 1]] if r.pos + 1 in m["idx"] else None)
    ex = exon_info(m, seg) if seg is not None else None
    in_frame_exon = bool(ex and ex["skippable"] and ex["in_frame"])
    if "nonsense" in mc:
        cc = coding_change(m, r.pos, r.ref, r.alt)
        stop = cc["mut_codon"] if cc and cc["mut_codon"] in STOP_CODONS else None
        lever = ("read-through: UGA stop" if stop == "TGA" else "read-through: UAG / UAA stop" if stop
                 else "exon skipping: lesion in an in-frame exon" if in_frame_exon else "no allele-specific lever")
        return lever, stop, in_frame_exon
    if "frameshift" in mc:
        esc, _ = nmd_status(r.gene, r.pos)
        return ("exon skipping: lesion in an in-frame exon" if in_frame_exon else "truncated protein escapes NMD" if esc
                else "no allele-specific lever"), None, in_frame_exon
    if re.search(r"splice_(donor|acceptor)", mc):
        return "splice: splice-switching candidate", None, in_frame_exon
    if "initiator_codon" in mc:      # a lost start codon is also annotated as missense, but then no protein is made from it
        return "start-loss: translation must restart downstream", None, False
    if "missense" in mc:
        return "missense: protein-level", None, False
    return "no allele-specific lever", None, False


LEVER_RATIONALE = {
    "read-through: UGA stop": "UGA is the stop codon most responsive to read-through agents in reporter assays, so this allele is a "
                             "candidate for the read-through hypotheses - to be measured, not assumed",
    "read-through: UAG / UAA stop": "a premature stop, but not the most read-through-prone codon: lower prior than UGA",
    "exon skipping: lesion in an in-frame exon": "the lesion sits in an exon whose length is a multiple of three, so skipping it would "
                                                "keep the reading frame (an antisense approach, not an approved drug)",
    "truncated protein escapes NMD": "predicted to escape nonsense-mediated decay, so a truncated protein is likely made: its stability "
                                     "and any dominant-negative effect are the questions",
    "splice: splice-switching candidate": "a canonical splice-site lesion, the class splice-switching oligonucleotides address",
    "missense: protein-level": "no nucleic-acid-level lever; the handle is protein abundance or function",
    "start-loss: translation must restart downstream": "the start codon is lost, so protein can only be made if translation restarts "
                                                       "at a later start site; none of the levers in this triage applies",
    "no allele-specific lever": "none of the allele-level levers in this triage apply",
}


def _short_allele(a, limit=24):
    """ClinVar holds some very long indel alleles; keep the file readable without losing the length."""
    a = str(a)
    return a if len(a) <= limit else f"{a[:limit]}...({len(a)} bp)"


def lever_atlas(genes, fa):
    rows = []
    for g in genes:
        if g not in GENE_TX or g not in GENE_POS:
            continue
        m = cds_model(GENE_TX[g][0], fa)
        cv = CLINVAR_PLP[CLINVAR_PLP.gene == g].drop_duplicates(["chrom", "pos", "ref", "alt"])
        for r in cv.itertuples():
            lever, stop, inframe = classify_clinvar(r, m)
            seg = m["seg"][m["idx"][r.pos]] if r.pos in m["idx"] else None
            ex = exon_info(m, seg) if seg is not None else None
            cc = coding_change(m, r.pos, r.ref, r.alt) if len(str(r.ref)) == 1 and len(str(r.alt)) == 1 else None
            rows.append(dict(gene=g, chrom=r.chrom, pos=int(r.pos), ref=_short_allele(r.ref), alt=_short_allele(r.alt),
                             key=f"{r.chrom}:{r.pos} {_short_allele(r.ref)}>{_short_allele(r.alt)}",
                             protein_change=(f"{cc['wt_aa']}{cc['residue']}" + (cc.get("mut_aa") or "")) if cc and cc.get("residue") else "",
                             exon=(ex or {}).get("exon") or "", n_exons=(ex or {}).get("n_exons") or "",
                             molecular_consequence=(r.mc or "").replace("|", ";"),
                             clinical_significance=r.clnsig, review_status=r.clnrevstat,
                             lever=lever, lever_rationale=LEVER_RATIONALE.get(lever, ""),
                             stop_codon_rna=(stop or "").replace("T", "U"), in_frame_exon=inframe,
                             note=("nonsense-tagged indel: no single stop codon computed" if "nonsense" in (r.mc or "") and not stop else "")))
    return pd.DataFrame(rows, columns=["gene", "chrom", "pos", "ref", "alt", "key", "protein_change", "exon", "n_exons",
                                       "molecular_consequence", "clinical_significance", "review_status",
                                       "lever", "lever_rationale", "stop_codon_rna", "in_frame_exon", "note"])


def lever_figure(gene, alleles, prot_len, ptc, mis):
    """Panel A: the protein, its domains, AlphaFold confidence and where the alleles fall.
    Panel B: ESM-2 tolerance of every substitution at the stop-codon residue (what read-through would insert)."""
    ctx = af_ctx(gene)
    row = esm_row(gene, ptc["residue"])[0] if ptc else None
    nrow = 2 if row is not None else 1
    fig, axes = plt.subplots(nrow, 1, figsize=(7.4, 5.8 if nrow == 2 else 3.0), gridspec_kw=dict(height_ratios=[1.15, 1][:nrow]))
    axes = np.atleast_1d(axes)
    ax = axes[0]
    if ctx:
        res = sorted(ctx["ca"])
        ax.plot(res, [ctx["ca"][r][1] for r in res], color=VIZ["ink2"], linewidth=1.0, label="AlphaFold confidence (pLDDT)")
        for f in ctx["feats"]:
            if f["type"] == "Domain" and f["start"] and f["end"]:
                ax.axvspan(f["start"], f["end"], ymin=0, ymax=0.06, color=VIZ["neutral"], linewidth=0)
                ax.text((f["start"] + f["end"]) / 2, 7, f["name"], ha="center", va="bottom", fontsize=7.5, color=VIZ["ink2"])
        dis = ctx.get("disease_missense") or []
        if dis:
            ax.scatter([d["residue"] for d in dis], [ctx["ca"][d["residue"]][1] for d in dis], s=30, color=VIZ["s3"], zorder=3,
                       edgecolors=VIZ["surface"], linewidths=1.5, label=f"literature-reported patient missense (UniProt, n={len(dis)})")
    if ptc and ptc.get("exon"):
        ex = ptc["exon"]
        ax.axvspan(ex["aa_from"], ex["aa_to"], color=VIZ["s2"], alpha=0.12, linewidth=0)
        ax.text((ex["aa_from"] + ex["aa_to"]) / 2, 103, f"exon {ex['exon']}", ha="center", va="bottom", fontsize=7.5, color=VIZ["ink2"])
    if ptc:
        ax.axvline(ptc["residue"], color=VIZ["s2"], linewidth=2, label=f"{ptc['label']} ({ptc['stop_rna']} stop)")
    if mis:
        ax.axvline(mis["residue"], color=VIZ["s1"], linewidth=2, label=mis["label"])
    ax.set_xlim(0, prot_len + 1)
    ax.set_ylim(0, 112)
    ax.set_xlabel(f"{gene} protein residue")
    ax.set_ylabel("pLDDT")
    ax.set_title(f"{gene}: where the child's alleles fall")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    if row is not None:
        ax = axes[1]
        wt = ptc["wt_aa"]
        items = sorted(((a, float(row[AA20.index(a)])) for a in AA20 if a != wt), key=lambda x: x[1])
        near = set(NEAR_COGNATE.get(ptc["stop_dna"], ""))
        cols = [VIZ["s2"] if a in near else VIZ["neutral"] for a, _ in items]
        ax.bar(range(len(items)), [x for _, x in items], color=cols, width=0.72)
        ax.set_xticks(range(len(items)))
        ax.set_xticklabels([f"{wt}{ptc['residue']}{a}" for a, _ in items], rotation=90, fontsize=7.5)
        for j, (a, x) in enumerate(items):
            if a in near:
                ax.text(j, x - 0.25, f"{x:.1f}", ha="center", va="top", fontsize=7.5, color=VIZ["ink"])
        if mis and mis.get("llr") is not None:
            ax.axhline(mis["llr"], color=VIZ["s1"], linewidth=1.5)
            ax.text(len(items) - 0.5, mis["llr"] + 0.25,
                    f"the child's {mis['short']} (LLR {mis['llr']:.1f}) — a different residue, shown for scale only",
                    ha="right", va="bottom",
                    fontsize=7.5, color=VIZ["ink2"])
        ax.axhline(0, color=VIZ["axis"], linewidth=0.8)
        ax.set_ylabel("ESM-2 LLR")
        ax.set_title(f"What read-through would insert at residue {ptc['residue']}: near-cognate amino acids (orange) "
                     f"vs all others", fontsize=9.5)
    fig.tight_layout()
    return save_fig(fig, "allele_levers")


with stage("Track 2: allele-level therapeutic triage"):
    if not CFG["RUN_TRACK2"]:
        raise RuntimeError("disabled in CFG")
    from pyfaidx import Fasta
    fa = Fasta(RES["ref"], as_raw=True, sequence_always_upper=True)
    TOP = HYPS[0]
    g = TOP["gene"]
    m = cds_model(GENE_TX[g][0], fa)
    t = TX[m["tid"]]
    prot = PROT.get(t["protein"], "").rstrip("*")
    L_out, ptc, mis = [], None, None
    for v in TOP["variants"]:
        cc = coding_change(m, int(v["pos"]), v["ref"], v["alt"])
        a = dict(key=f"{v['chrom']}:{int(v['pos'])} {v['ref']}>{v['alt']}", hgvsc=v.get("hgvsc", ""), hgvsp=v.get("hgvsp", ""),
                 csq=v["csq"], c=cc["c"] if cc else None)
        if cc and cc["mut_codon"] in STOP_CODONS:
            ex = exon_info(m, cc["segment"])
            seg_res = range(ex["aa_from"], ex["aa_to"] + 1)
            ex.update(mean_plddt=af_at(g, seg_res), features=features_overlapping(g, ex["aa_from"], ex["aa_to"]))
            try:
                ex["natural_skipping_isoforms"] = natural_skipping(g, ex) if ex["skippable"] and ex["in_frame"] else None
            except Exception as e:      # noqa: BLE001
                ex["natural_skipping_isoforms"] = None
                log(f"Ensembl isoform lookup skipped ({str(e)[:100]})")
            row, flat = esm_row(g, cc["residue"])
            nc = {}
            if row is not None:
                for aa in NEAR_COGNATE[cc["mut_codon"]]:
                    x = float(row[AA20.index(aa)])
                    nc[aa] = dict(llr=round(x, 2), pct_less_tolerated=round(float((flat < x).mean()) * 100, 1) if flat is not None else None,
                                  restores_wild_type=aa == cc["wt_aa"])
            ptc = a
            a.update(kind="premature stop", residue=cc["residue"], wt_aa=cc["wt_aa"], stop_dna=cc["mut_codon"],
                     stop_rna=cc["mut_codon"].replace("T", "U"), plus4=cc["plus4"],
                     readthrough_rank=STOP_ORDER[cc["mut_codon"]], near_cognate_aas=NEAR_COGNATE[cc["mut_codon"]],
                     plus4_rna=cc["plus4"].replace("T", "U"), native_residues_kept=cc["residue"] - 1,
                     nmd_escape=v.get("nmd_escape"), plddt=af_at(g, [cc["residue"]]), near_cognate=nc, exon=ex,
                     label=v.get("hgvsp") or a["key"], truncated_product=f"1-{cc['residue'] - 1}",
                     truncated_kda=round(peptide_mass(prot[:cc["residue"] - 1]) / 1000, 1) if prot else None,
                     truncated_keeps=features_overlapping(g, 1, cc["residue"] - 1),
                     truncated_loses=features_overlapping(g, cc["residue"], len(prot)))
        elif cc and v["csq"] == "missense":
            row, flat = esm_row(g, cc["residue"])
            st = S.get("structure", {}) if (af_ctx(g)
                                            and (S.get("structure") or {}).get("residue") == cc["residue"]) else {}
            mis = a
            a.update(kind="missense", residue=cc["residue"], wt_aa=cc["wt_aa"], mut_aa=cc["mut_aa"],
                     llr=round(float(row[AA20.index(cc["mut_aa"])]), 2) if row is not None else None,
                     plddt=af_at(g, [cc["residue"]]), domain=st.get("domain"), burial_percentile=st.get("burial_percentile_in_domain"),
                     nearest_patient_missense=(st.get("disease_missense_uniprot") or [])[:3],
                     label=v.get("hgvsp") or a["key"], short=f"{cc['wt_aa']}{cc['residue']}{cc['mut_aa']}")
        L_out.append(a)
    S["levers"] = dict(gene=g, transcript=t.get("version") or m["tid"], protein_length=len(prot), alleles=L_out)
    md("**Allele-level triage**\n\n" + df_md(pd.DataFrame([{k: x.get(k) for k in ("key", "hgvsp", "kind", "residue", "stop_rna", "plus4",
                                                                                   "nmd_escape", "plddt", "llr")} for x in L_out])))

    # ---- the pair: cDNA phasing amplicon, allele-resolved peptides, path to reclassification -----------------------------
    if TOP["model"] == "compound_het" and all(x.get("c") for x in L_out):
        c1, c2 = sorted(x["c"] for x in L_out)
        gd = abs(int(TOP["variants"][0]["pos"]) - int(TOP["variants"][1]["pos"]))
        S["levers"]["pair"] = dict(cdna_distance_nt=c2 - c1, genomic_distance_bp=gd, amplicon_nt=c2 - c1 + 300,
                                   fold_shorter=round(gd / max(c2 - c1, 1), 1))
        log(f"the alleles are {gd:,} bp apart in the genome but {c2 - c1:,} nt apart in the mRNA: one RT-PCR amplicon "
            f"of ~{c2 - c1 + 300:,} nt from the child's RNA phases them")
        try:
            prim = design_rtpcr_primers(m, c1, c2, fa)
            if prim:
                S["levers"]["pair"]["primers"] = prim
                S["levers"]["pair"]["amplicon_nt"] = prim[0]["product_nt"]
                md("**RT-PCR primers for the phasing / NMD-escape amplicon** (Primer3; no variant of the child and no common SNP "
                   "under either primer)\n\n" + df_md(pd.DataFrame(prim)))
        except Exception as e:      # noqa: BLE001
            log(f"RT-PCR primer design skipped ({type(e).__name__}: {str(e)[:120]})")
    if ptc and mis and prot and TOP["model"] == "compound_het":
        def put(seq, res, aa):
            return seq[:res - 1] + aa + seq[res:]
        mis_prot = put(prot, mis["residue"], mis["mut_aa"])          # missense allele: wild type at the stop site
        rt = [aa for aa in NEAR_COGNATE[ptc["stop_dna"]] if aa != ptc["wt_aa"]]
        restores = ptc["wt_aa"] in NEAR_COGNATE[ptc["stop_dna"]]
        downstream = mis["residue"] > ptc["residue"]
        rt_any = put(prot, ptc["residue"], rt[0]) if rt else prot  # stop allele, full length: wild type at the missense site
        peps = site_peptides({f"{mis['wt_aa']}{mis['residue']} (reference)":
                              (rt_any, "the stop allele — only by read-through, since the truncated product ends before this residue"
                               if downstream else "the stop allele (truncated product if NMD is relieved, or read-through)"),
                              mis["short"]: (mis_prot, "the missense allele")}, mis["residue"], avoid=ptc["residue"])
        forms = {f"{ptc['wt_aa']}{ptc['residue']} (reference)":
                 (mis_prot, "the missense allele" + (" (or read-through restoring the wild-type residue)" if restores else ""))}
        for aa in rt:
            forms[f"{ptc['wt_aa']}{ptc['residue']}{aa}"] = (put(prot, ptc["residue"], aa), f"read-through of the stop allele inserting {aa}")
        peps += site_peptides(forms, ptc["residue"], avoid=mis["residue"]) if rt else []
        S["levers"]["peptides"] = peps
        S["levers"]["peptide_logic"] = dict(missense_downstream_of_stop=downstream, near_cognate_restores_wild_type=restores)
        if peps:
            md("**Allele-resolved peptides for targeted mass spectrometry**\n\n" + df_md(pd.DataFrame(peps), floatfmt="{:.3f}"))
    if TOP["model"] == "compound_het":
        weak = min(TOP["variants"], key=lambda v: v["p_path"])
        strong = max(TOP["variants"], key=lambda v: v["p_path"])
        wk = f"{weak['chrom']}:{int(weak['pos'])} {weak['ref']}>{weak['alt']}"
        cls = next((x.get("classification") for x in (S.get("top") or {}).get("variants", []) if x.get("key") == wk), None)
        S["levers"]["reclassification"] = dict(allele=weak.get("hgvsp") or wk, current_class=cls,
                                               path=reclassification_path(weak, strong, float(TOP.get("pm3") or 0.0)))
        md("**Path to reclassifying the weaker allele**\n\n" + df_md(pd.DataFrame(S["levers"]["reclassification"]["path"])))

    # ---- scalability: the same triage over every ClinVar P/LP allele in the disorder's panel ---------------------------
    panel = [g] + [x for x in (S.get("pheno") or {}).get("dx_panel", []) if x != g]
    ATLAS = lever_atlas(panel, fa)
    if len(ATLAS):
        tab = ATLAS.groupby(["gene", "lever"]).size().unstack(fill_value=0)
        tab = tab[[c for c in LEVER_ORDER if c in tab.columns]]
        tab["total"] = tab.sum(1)
        tab = tab.sort_values("total", ascending=False)
        non = ATLAS[ATLAS.stop_codon_rna != ""]
        S["levers"]["atlas"] = dict(n_alleles=int(len(ATLAS)), n_genes=int(ATLAS.gene.nunique()),
                                    by_gene=tab.reset_index().to_dict("records"),
                                    by_lever=ATLAS.lever.value_counts().reindex([c for c in LEVER_ORDER if c in set(ATLAS.lever)]).to_dict(),
                                    nonsense_by_stop=non.stop_codon_rna.value_counts().to_dict(),
                                    top_gene_nonsense_by_stop=non[non.gene == g].stop_codon_rna.value_counts().to_dict(),
                                    n_no_criteria=int(ATLAS.review_status.fillna("").str.contains("no_assertion|no_classification").sum()))
        ATLAS.to_csv(OUT / "lever_atlas.csv", index=False)       # reusable product: one row per public ClinVar allele
        md(f"**Allele-specific lever atlas: {len(ATLAS)} ClinVar P/LP alleles in {ATLAS.gene.nunique()} genes of the panel**\n\n"
           + df_md(tab.reset_index()) + "\n\nSaved as `lever_atlas.csv` (public ClinVar alleles only).")
    try:
        lever_figure(g, L_out, len(prot), ptc, mis)
    except Exception as e:      # noqa: BLE001
        log(f"lever figure skipped ({type(e).__name__}: {str(e)[:120]})")
    log(f"allele-level triage: {json.dumps({k: v for k, v in S['levers'].items() if k in ('pair', 'atlas')}, default=str)[:600]}")

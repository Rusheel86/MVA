# %% [markdown]
# ## 10 · Deep annotation of the short list and final ranking
#
# Only the variants that appear in the leading hypotheses (a few dozen alleles) are sent to public annotation services:
#
# * **Ensembl VEP REST** — HGVS on the MANE transcript, LOFTEE, pre-computed SpliceAI, CADD, REVEL, literature links;
# * **gnomAD v4.1 joint** (remote tabix) — exact allele counts, group-maximum filtering AF, homozygotes, and a check
#   for common overlapping indels (the classic way a mis-represented common allele masquerades as a novel one);
# * **Genomics England PanelApp** — is the leading gene a green gene on NHS virtual panels, and with which mode of inheritance;
# * **chance-bystander analysis** (gnomAD gene-level counts) — how often would an unrelated person carry, by coincidence,
#   a rare variant at least as damaging as the weaker allele of the leading pair?
#
# The evidence model is then re-run with the richer annotations and the submission file is rewritten.

# %%
import requests

HTTP = requests.Session()
HTTP.headers.update({"User-Agent": "duorank-notebook/1.0 (research; rare-disease hackathon)"})


def http_json(method, url, tries=4, timeout=120, **kw):
    last = None
    for i in range(tries):
        try:
            r = HTTP.request(method, url, timeout=timeout, **kw)
            if r.status_code == 200:
                return r.json()
            last = f"HTTP {r.status_code}: {r.text[:150]}"
            if r.status_code in (400, 404):
                break
        except Exception as e:      # noqa: BLE001
            last = f"{type(e).__name__}: {str(e)[:150]}"
        time.sleep(3 * (i + 1))
    raise RuntimeError(f"{url.split('?')[0]} failed ({last})")


def shortlist_variants(hyps, recs, n_hyp=40, n_genes=3):
    keep = {}
    for h in hyps[:n_hyp]:
        for v in h["variants"]:
            keep[(vkey(v), v["gene"])] = v
    top_genes = list(dict.fromkeys(h["gene"] for h in hyps[:20]))[:n_genes]
    for v in recs:
        if v["gene"] in top_genes:
            keep[(vkey(v), v["gene"])] = v
    return list(keep.values()), top_genes


def vep_annotate(variants):
    """Ensembl VEP REST (GRCh38) for a short list; returns {vkey: raw VEP record}."""
    uniq = {vkey(v): v for v in variants}
    strings = {f"{k[0][3:] if k[0] != 'chrM' else 'MT'} {k[1]} . {k[2]} {k[3]} . . .": k for k in uniq}
    out, items = {}, list(strings)
    for i in range(0, len(items), 150):
        res = http_json("POST", "https://rest.ensembl.org/vep/homo_sapiens/region",
                        headers={"Content-Type": "application/json", "Accept": "application/json"},
                        params={"hgvs": 1, "mane": 1, "canonical": 1, "variant_class": 1, "SpliceAI": 1, "CADD": 1, "AlphaMissense": 1,
                                "REVEL": 1, "LoF": 1, "numbers": 1}, json={"variants": items[i:i + 150]}, timeout=600)
        for o in res:
            k = strings.get(o.get("input"))
            if k:
                out[k] = o
    return out


def apply_vep(recs, vep):
    n = 0
    for v in recs:
        o = vep.get(vkey(v))
        if not o:
            continue
        tcs = o.get("transcript_consequences") or []
        mane = [t for t in tcs if t.get("mane_select") and t.get("gene_symbol") == v["gene"]] or [t for t in tcs if t.get("mane_select")]
        t = (mane or [t for t in tcs if t.get("canonical")] or [None])[0]
        v["rsids"] = [c.get("id") for c in o.get("colocated_variants", []) if str(c.get("id", "")).startswith("rs")]
        v["pubmed"] = sorted({p for c in o.get("colocated_variants", []) for p in (c.get("pubmed") or [])})[:12]
        if not t:
            continue
        sp = t.get("spliceai") or {}
        if sp:
            v["spliceai"] = max(float(sp.get(x) or 0) for x in ("DS_AG", "DS_AL", "DS_DG", "DS_DL"))
            v["spliceai_detail"] = {x: sp.get(x) for x in ("DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL")}
        for src, dst in (("revel", "revel"), ("cadd_phred", "cadd")):
            if t.get(src) is not None:
                v[dst] = float(t[src])
        am = (t.get("alphamissense") or {}).get("am_pathogenicity")
        if am is not None and _num(v.get("am")) is None:
            v["am"] = float(am)
        v["hgvsc"] = t.get("hgvsc", "")
        if t.get("hgvsp"):
            v["hgvsp"] = t["hgvsp"].split(":")[-1]
        v["mane_refseq"] = t.get("mane_select", "")
        v["exon"], v["intron"] = t.get("exon", ""), t.get("intron", "")
        v["loftee"] = t.get("lof", "")
        v["sift"], v["polyphen"] = t.get("sift_prediction", ""), t.get("polyphen_prediction", "")
        v["vep_consequence"] = ",".join(t.get("consequence_terms") or [])
        n += 1
    return n


GNOMAD_FMT = ("%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%INFO/AC_joint\\t%INFO/AN_joint\\t%INFO/AF_joint\\t%INFO/nhomalt_joint"
              "\\t%INFO/AF_grpmax_joint\\t%INFO/grpmax_joint\\t%INFO/fafmax_faf95_max_joint\\n")


def gnomad_exact(variants):
    """gnomAD v4.1 joint: exact records for the short list (remote tabix) + common overlapping indels within 30 bp."""
    out, by_chr = {}, collections.defaultdict(dict)
    for v in variants:
        by_chr[v["chrom"]][vkey(v)] = v
    for chrom, vs in by_chr.items():
        if chrom not in PRIMARY:
            continue
        regions = ",".join(f"{chrom}:{max(1, k[1] - 30)}-{k[1] + len(k[2]) + 30}" for k in vs)
        txt = sh(f"bcftools query -r '{regions}' -f '{GNOMAD_FMT}' '{GNOMAD_V4_JOINT.format(chrom=chrom)}'", name=f"gnomad_{chrom}",
                 capture=True, timeout=900)
        rows = {}
        for line in txt.splitlines():
            c = line.split("\t")
            if len(c) >= 12:
                rows[(c[0], int(c[1]), c[2], c[3])] = c
        for k in vs:
            c = rows.get(k)
            rec = dict(present=c is not None, source="gnomAD v4.1 joint (remote tabix)")
            if c is not None:
                rec.update(filter=c[4], ac=_num(c[5]), an=_num(c[6]), af=_num(c[7]), nhomalt=_num(c[8]), af_grpmax=_num(c[9]),
                           grpmax=c[10] if c[10] != "." else "", faf95=_num(c[11]))
            if len(k[2]) != len(k[3]):       # indel: a common indel overlapping the same span suggests a representation artefact
                lo, hi = k[1], k[1] + len(k[2])
                ov = [(_num(r[7]) or 0.0) for kk, r in rows.items() if kk != k and len(kk[2]) != len(kk[3])
                      and kk[1] <= hi and kk[1] + len(kk[2]) >= lo]
                rec["common_overlap"] = max(ov, default=0.0)
            out[k] = rec
    return out


def gnomad_graphql(variants):
    """Fallback when remote tabix is unavailable: the public gnomAD GraphQL API (exact matches only)."""
    q = ("query($v: String!) { variant(variantId: $v, dataset: gnomad_r4) { exome { ac an homozygote_count } "
         "genome { ac an homozygote_count } joint { ac an homozygote_count faf95 { popmax popmax_population } } } }")
    out = {}
    for k in {vkey(v) for v in variants}:
        try:
            d = http_json("POST", "https://gnomad.broadinstitute.org/api", json={"query": q, "variables": {"v": f"{k[0][3:]}-{k[1]}-{k[2]}-{k[3]}"}})
            g = (d.get("data") or {}).get("variant")
            rec = dict(present=g is not None, source="gnomAD v4.1 (GraphQL API)")
            if g:
                j = g.get("joint") or {}
                ac = j.get("ac") if j else sum((g.get(x) or {}).get("ac", 0) for x in ("exome", "genome"))
                an = j.get("an") if j else sum((g.get(x) or {}).get("an", 0) for x in ("exome", "genome"))
                rec.update(ac=ac, an=an, af=(ac / an if an else None), nhomalt=j.get("homozygote_count"),
                           faf95=(j.get("faf95") or {}).get("popmax"), grpmax=(j.get("faf95") or {}).get("popmax_population") or "")
            out[k] = rec
            time.sleep(0.4)
        except Exception:       # noqa: BLE001
            continue
    return out


def apply_gnomad(recs, gn):
    for v in recs:
        g = gn.get(vkey(v))
        if not g:
            continue
        v["gnomad"] = g
        if g["present"]:
            vals = [x for x in (g.get("faf95"), g.get("af")) if x is not None and x > 0]
            v["af"] = max(vals) if vals else None
            v["nhomalt"] = g.get("nhomalt")
        else:
            v["af"], v["nhomalt"] = None, None
        if g.get("common_overlap"):
            v["common_overlap"] = g["common_overlap"]


def am_region(chrom, start, end):
    """All AlphaMissense rows in a region -> {(pos, ref, alt): (score, protein_variant, transcript)}."""
    if "alphamissense" not in RES:
        return {}
    txt = sh(f"tabix '{RES['alphamissense']}' {chrom}:{start}-{end}", name="am_region", capture=True)
    out = {}
    for line in txt.splitlines():
        c = line.split("\t")
        if len(c) >= 9:
            out[(int(c[1]), c[2], c[3])] = (float(c[8]), c[7], c[6])
    return out


def bystander_analysis(gene, weaker):
    """Chance-bystander model for the weaker allele of a pair, from gnomAD v4.1 gene-level allele counts."""
    q = ("query($s: String!) { gene(gene_symbol: $s, reference_genome: GRCh38) { variants(dataset: gnomad_r4) { "
         "variant_id consequence exome { ac an } genome { ac an } } } }")
    d = http_json("POST", "https://gnomad.broadinstitute.org/api", json={"query": q, "variables": {"s": gene}}, timeout=180)
    variants = ((d.get("data") or {}).get("gene") or {}).get("variants") or []
    if not variants:
        raise RuntimeError("no gnomAD gene-level data")
    chrom, start, end = GENE_POS[gene][0], GENE_POS[gene][1], GENE_POS[gene][2]
    am = am_region(chrom, start, end)
    cls_map = {"missense": {"missense_variant"}, "stop_gained": {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"},
               "frameshift": {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"},
               "inframe_indel": {"inframe_deletion", "inframe_insertion"}}
    wanted = cls_map.get(weaker["csq"], None)
    q_all = q_strict = 0.0
    n_all = n_strict = 0
    am_w = _num(weaker.get("am"))
    band_lo = max((thr for thr, _ in AM_BANDS if am_w is not None and am_w >= thr), default=None)
    for g in variants:
        if wanted is not None and g["consequence"] not in wanted:
            continue
        ac = sum((g.get(x) or {}).get("ac") or 0 for x in ("exome", "genome"))
        an = sum((g.get(x) or {}).get("an") or 0 for x in ("exome", "genome"))
        af = ac / an if an else 0.0
        if not 0 < af < 1e-4:
            continue
        q_all += af
        n_all += 1
        if weaker["csq"] == "missense" and band_lo is not None:
            _, pos, ref, alt = g["variant_id"].split("-")
            s = am.get((int(pos), ref, alt))
            if s and s[0] >= band_lo:
                q_strict += af
                n_strict += 1
    res = dict(gene=gene, variant_class=weaker["csq"], n_rare_variants=n_all, cumulative_af=q_all,
               carrier_probability=1 - (1 - q_all) ** 2)
    if weaker["csq"] == "missense" and band_lo is not None:
        res.update(am_band_threshold=band_lo, n_rare_at_least_as_damaging=n_strict, cumulative_af_strict=q_strict,
                   carrier_probability_strict=1 - (1 - q_strict) ** 2)
    # P(an unrelated person carries such an allele by chance): the conservative estimate uses every rare variant of the
    # same class; the strict one (only alleles at least as damaging) is reported as a sensitivity bound.
    pb = res["carrier_probability"]
    res["posterior_second_allele"] = {f"prior_{p}": round(p / (p + (1 - p) * pb), 4) for p in (0.2, 0.5, 0.8)}
    return res


def acmg_table(v, extra_points=0.0, extra_label=None):
    rows = [dict(criterion=c, points=p, basis=txt) for c, p, txt in v["evidence"]]
    if extra_points:
        rows.append(dict(criterion=extra_label or "PM3", points=extra_points, basis="partner allele (likely) pathogenic; phase unknown"))
    total = v["points"] + extra_points
    post = posterior(v, extra_points)
    by_post = ("Pathogenic" if post >= 0.99 else "Likely pathogenic" if post >= 0.90 else "Uncertain significance" if post > 0.10
               else "Likely benign" if post > 0.001 else "Benign")
    # the stronger of the posterior class (which carries ClinVar floors / caps) and the point-threshold class
    label = min(by_post, points_class(total), key=ACMG_CLASSES.index) if post > 0.10 else by_post
    return pd.DataFrame(rows), total, post, label


ACMG_CLASSES = ["Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "Benign"]


def points_class(points):
    """ACMG/AMP class from the naturally scaled point system (Tavtigian 2020): P >= 10, LP 6-9, VUS 0-5, LB -1..-6, B <= -7."""
    return ("Pathogenic" if points >= 10 else "Likely pathogenic" if points >= 6 else "Uncertain significance" if points >= 0
            else "Likely benign" if points >= -6 else "Benign")


def rescore(recs):
    for v in recs:
        inh = set((GENES_DX.get(v["gene"]) or {}).get("inheritance_all") or [])
        lof_ok = v["gene"] in LOF_GENES or bool(inh & {"AR", "XLR", "XL"})
        v["points"], v["floor"], v["cap"], v["evidence"] = variant_evidence(v, lof_ok=lof_ok)
        v["p_path"] = posterior(v)


with stage("deep annotation of the short list"):
    SHORT, TOP_GENES = shortlist_variants(HYPS, RECS)
    log(f"short list: {len({vkey(v) for v in SHORT})} alleles; leading genes: {TOP_GENES}")
    try:
        VEP = vep_annotate(SHORT)
        log(f"Ensembl VEP REST: {apply_vep(RECS, VEP)} variant rows annotated")
        S["vep_ok"] = True
    except Exception as e:      # noqa: BLE001
        log(f"Ensembl VEP REST unavailable ({str(e)[:160]}); continuing with offline annotations")
        S["vep_ok"] = False
    try:
        GN = gnomad_exact(SHORT)
    except Exception as e:      # noqa: BLE001
        log(f"gnomAD remote tabix unavailable ({str(e)[:160]}); using the GraphQL API")
        GN = gnomad_graphql(SHORT)
    apply_gnomad(RECS, GN)
    S["gnomad_exact_n"] = len(GN)
    rescore(RECS)
    HYPS = build_hypotheses(RECS, GENES_DX, PHENO_TOP, W)
    S["ranking_stage"] = "after deep annotation (VEP, gnomAD v4.1 exact)"
    S["top_gene"] = HYPS[0]["gene"]
    md("**Top hypotheses after deep annotation**\n\n" + df_md(public_table(hyp_table(HYPS, 10))))

# %%
def analyse_leader():
    """Evidence tables, panel status and chance-bystander analysis for the current leading hypothesis."""
    global TOP
    TOP = HYPS[0]
    S["top_key"] = sorted(TOP["key"])
    S["top_gene"] = TOP["gene"]
    S["top"] = dict(gene=TOP["gene"], model=TOP["model"], phase=TOP["phase"], epcr=round(TOP["epcr"], 4), variants=[])
    info = GENES_DX.get(TOP["gene"]) or {}
    S["top"].update(disease=info.get("disease_name"), dx_disease=info.get("dx_disease"), inheritance=info.get("inheritance"),
                    pheno_score=round(info.get("pheno", 0.0), 3), pheno_rank_blind=GENE_RANK_BLIND.get(TOP["gene"]))
    ordered = sorted(TOP["variants"], key=lambda v: -v["p_path"])
    for i, v in enumerate(ordered):
        extra = TOP.get("pm3", 0.0) if (TOP["model"] == "compound_het" and i == 1) else (1.0 if TOP["model"] in ("homozygous", "hemizygous") else 0.0)
        tab, total, post, label = acmg_table(v, extra, "PM3" if TOP["model"] == "compound_het" else "PM3 (homozygous)")
        g = v.get("gnomad") or {}
        S["top"]["variants"].append(dict(
            key=f"{v['chrom']}:{int(v['pos'])} {v['ref']}>{v['alt']}", chrom=v["chrom"], pos=int(v["pos"]), ref=v["ref"], alt=v["alt"],
            csq=v["csq"], hgvsc=v.get("hgvsc", ""), hgvsp=v.get("hgvsp", ""), mane=v.get("mane_refseq", ""), exon=v.get("exon", ""),
            zyg=v["zyg"], ad=f"{v.get('ad_ref')}/{v.get('ad_alt')}", gq=v.get("gq"), filt=v.get("filt"), pid=v.get("pid"),
            clnsig=v.get("clnsig", ""), clnrevstat=v.get("clnrevstat", ""), rsids=v.get("rsids", []), pubmed=v.get("pubmed", []),
            am=_num(v.get("am")), revel=_num(v.get("revel")), cadd=_num(v.get("cadd")), spliceai=_num(v.get("spliceai")),
            loftee=v.get("loftee", ""), nmd_escape=v.get("nmd_escape"), frac_lost=_num(v.get("frac_lost")), aa_pos=v.get("aa_pos"),
            gnomad=g, points=total, posterior=round(post, 4), classification=label, evidence=tab.to_dict("records")))
        md(f"**{v['gene']} {v.get('hgvsc') or ''} {v.get('hgvsp') or ''}** — `{v['chrom']}:{int(v['pos'])} {v['ref']}>{v['alt']}` · "
           f"{total:g} points · posterior {post:.3f} · **{label}** (automated)\n\n" + df_md(tab))
    try:
        res = http_json("GET", f"https://panelapp.genomicsengland.co.uk/api/v1/genes/{TOP['gene']}/?format=json").get("results", [])
        S["panelapp"] = [dict(panel=x["panel"]["name"], id=x["panel"]["id"], version=x["panel"].get("version", ""),
                              mode=x.get("mode_of_inheritance", "")) for x in res if str(x.get("confidence_level")) == "3"]
        log(f"PanelApp: {TOP['gene']} is a green gene on {len(S['panelapp'])} Genomics England panels")
    except Exception as e:      # noqa: BLE001
        log(f"PanelApp lookup skipped ({str(e)[:120]})")
    if TOP["model"] == "compound_het":
        try:
            S["bystander"] = bystander_analysis(TOP["gene"], ordered[1])
            b = S["bystander"]
            log(f"chance-bystander: P(an unrelated person carries a rare {b['variant_class']} variant in {b['gene']}) = "
                f"{b['carrier_probability']:.4f}" + (f"; at least as damaging by AlphaMissense: {b['carrier_probability_strict']:.4f}"
                                                     if "carrier_probability_strict" in b else "")
                + f"; posterior that the weaker allele is the second causal allele: {b['posterior_second_allele']}")
        except Exception as e:      # noqa: BLE001
            log(f"chance-bystander analysis skipped ({str(e)[:160]})")


with stage("leading hypothesis: evidence tables, panel status, chance-bystander analysis"):
    analyse_leader()

# %%
with stage("write + validate submission (after deep annotation)", required=True):
    SUB_ROWS, parsed = write_submission(HYPS)
    chk = verify_in_source([v for h in SUB_ROWS for v in h["variants"]])
    S["source_check"] = chk
    if not DEMO:
        assert all(c["found"] and not c["synthetic"] for c in chk), f"submitted allele not found in the source VCF: {chk}"
    else:
        S["demo_score"] = official_score(parsed, frozenset(SPIKE_CASES[0]["variants"]))
        log(f"DEMO self-score against the spiked truth (organisers' scorer): {S['demo_score']}")
    md("**Submission rows**\n\n" + df_md(pd.read_csv(SUBMISSION_CSV, dtype=str, keep_default_na=False)))

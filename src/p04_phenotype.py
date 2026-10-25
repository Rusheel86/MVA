# %% [markdown]
# ## 4 · Phenotype engine
#
# Every Mendelian disease in the HPO annotation corpus is scored against the proband's terms with an
# information-content (Resnik) best-match similarity; diseases roll up to genes. Two rankings are kept:
#
# * **HPO-only** — blind to the diagnosis; this is what would be available for an undiagnosed child.
# * **HPO + stated clinical indication** — genes whose disease name matches `CFG["CLINICAL_DIAGNOSIS"]` form a
#   virtual panel (the way a clinical laboratory applies a panel for the referral indication).
#
# Terms that describe relatives are kept but down-weighted.

# %%
class HPO:
    def __init__(self, obo_path):
        self.name, self.parents, self.alt, self.version = {}, collections.defaultdict(set), {}, "n/a"
        cur = None
        with open(obo_path, encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("data-version:"):
                    self.version = line.split(":", 1)[1].strip()
                    continue
                if line == "[Term]":
                    cur = {}
                    continue
                if line.startswith("["):
                    cur = None
                    continue
                if cur is None:
                    continue
                if line.startswith("id: "):
                    cur["id"] = line[4:]
                elif line.startswith("name: "):
                    self.name[cur["id"]] = line[6:]
                elif line.startswith("is_a: "):
                    self.parents[cur["id"]].add(line[6:].split(" ! ")[0].strip())
                elif line.startswith("alt_id: "):
                    self.alt[line[8:].strip()] = cur["id"]
                elif line.startswith("replaced_by: "):
                    self.alt[cur["id"]] = line[13:].strip()
        self._anc = {}

    def canon(self, t):
        return self.alt.get(t, t)

    def ancestors(self, t):
        t = self.canon(t)
        if t in self._anc:
            return self._anc[t]
        out = {t}
        for p in self.parents.get(t, ()):
            out |= self.ancestors(p)
        self._anc[t] = out
        return out


INHERITANCE_TERMS = {
    "HP:0000007": "AR", "HP:0000006": "AD", "HP:0001417": "XL", "HP:0001419": "XLR", "HP:0001423": "XLD",
    "HP:0001427": "MT", "HP:0001450": "YL", "HP:0012275": "AD", "HP:0012274": "AD", "HP:0001466": "AD",
    "HP:0025352": "AD",
}


class PhenoEngine:
    def __init__(self, obo, hpoa, genes_to_disease):
        self.hpo = HPO(obo)
        self.d_name, self.d_terms, self.d_inh = {}, collections.defaultdict(set), collections.defaultdict(set)
        with open(hpoa, encoding="utf-8") as f:
            ix = None
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                if ix is None:
                    ix = {k: i for i, k in enumerate(c)}
                    continue
                did, hid, aspect = c[ix["database_id"]], self.hpo.canon(c[ix["hpo_id"]]), c[ix["aspect"]]
                self.d_name[did] = c[ix["disease_name"]]
                if c[ix["qualifier"]].upper() == "NOT":
                    continue
                if aspect == "P":
                    self.d_terms[did].add(hid)
                elif aspect == "I" and hid in INHERITANCE_TERMS:
                    self.d_inh[did].add(INHERITANCE_TERMS[hid])
        self.g2d, self.d2g = collections.defaultdict(set), collections.defaultdict(set)
        with open(genes_to_disease, encoding="utf-8") as f:
            ix = {k: i for i, k in enumerate(f.readline().rstrip("\n").split("\t"))}
            for line in f:
                c = line.rstrip("\n").split("\t")
                if c[ix["association_type"]] == "POLYGENIC":       # susceptibility loci are not Mendelian explanations
                    continue
                self.g2d[c[ix["gene_symbol"]]].add(c[ix["disease_id"]])
                self.d2g[c[ix["disease_id"]]].add(c[ix["gene_symbol"]])
        n = len(self.d_terms)
        cnt = collections.Counter()
        self.d_closure = {}
        for did, ts in self.d_terms.items():
            clo = set()
            for t in ts:
                clo |= self.hpo.ancestors(t)
            self.d_closure[did] = clo
            cnt.update(clo)
        self.ic = {t: -math.log(k / n) for t, k in cnt.items()}
        self.n_diseases = n

    def term_ic(self, q):
        q = self.hpo.canon(q)
        if q in self.ic:
            return self.ic[q]
        return max((self.ic.get(a, 0.0) for a in self.hpo.ancestors(q)), default=0.0)

    def _best(self, q, closure):
        best = 0.0
        for a in self.hpo.ancestors(q):
            if a in closure:
                v = self.ic.get(a, 0.0)
                if v > best:
                    best = v
        return best

    def score_diseases(self, terms, weights=None):
        """score = 0.8 * (share of the patient's information content recovered by the disease)
                 + 0.2 * (share of the disease's information content present in the patient)."""
        terms = [self.hpo.canon(t) for t in terms]
        w = {self.hpo.canon(k): v for k, v in (weights or {}).items()}
        denom = sum(w.get(q, 1.0) * self.term_ic(q) for q in terms) or 1.0
        q_closure = set()
        for q in terms:
            q_closure |= self.hpo.ancestors(q)
        out = {}
        for did, clo in self.d_closure.items():
            fwd = sum(w.get(q, 1.0) * self._best(q, clo) for q in terms) / denom
            if fwd <= 0:
                continue
            rev_den = sum(self.ic.get(d, 0.0) for d in self.d_terms[did])
            rev = sum(self._best(d, q_closure) for d in self.d_terms[did]) / rev_den if rev_den else 0.0
            out[did] = (0.8 * fwd + 0.2 * rev, fwd, rev)
        return out

    def score_genes(self, terms, weights=None, diagnosis=None):
        ds = self.score_diseases(terms, weights)
        ranked = sorted(ds, key=lambda d: -ds[d][0])
        pct = {did: 1.0 - i / max(len(ranked) - 1, 1) for i, did in enumerate(ranked)}
        pat = (diagnosis or "").strip().lower()
        genes = {}
        for g, dids in self.g2d.items():
            scored = sorted((d for d in dids if d in ds), key=lambda d: -ds[d][0])
            dx = [d for d in dids if pat and pat in self.d_name.get(d, "").lower()]
            best = scored[0] if scored else (dx[0] if dx else None)
            inh_all = set()
            for d in dids:
                inh_all |= self.d_inh.get(d, set())
            # mode of inheritance: of the best-matching disease (or of the diagnosis-matched disease); gene-level union as fallback
            inh = set()
            for d in (dx or scored[:1]):
                inh |= self.d_inh.get(d, set())
            genes[g] = dict(pheno=ds[best][0] if best in ds else 0.0, pct=pct.get(best, 0.0), disease=best,
                            disease_name=self.d_name.get(best), inheritance=sorted(inh or inh_all),
                            inheritance_all=sorted(inh_all), dx_match=bool(dx),
                            dx_disease=self.d_name.get(dx[0]) if dx else None)
        return genes, ds


with stage("phenotype engine (load HPO)", required=True):
    PE = PhenoEngine(RES["hp_obo"], RES["hpoa"], RES["genes_to_disease"])
    S["hpo_version"] = PE.hpo.version
    log(f"HPO release {PE.hpo.version}: {PE.n_diseases:,} diseases with phenotype annotations, {len(PE.g2d):,} genes")

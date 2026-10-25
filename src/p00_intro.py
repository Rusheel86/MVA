# %% [markdown]
# # DuoRank — phenotype-driven compound-heterozygote discovery for a rare-disease genome
#
# **Rare Disease, Real Kid: MVA Hackathon 2026** (Sage Bionetworks · MVA Society · Hugging Face · BEACON) — Track 1 (variant prediction) with a Track 2 (drug repurposing) evidence module.
#
# This notebook is a fully automated, public-data-only pipeline. Nothing about the answer is hard-coded: the ranked
# submission is derived from the proband's VCF and HPO terms, and the same code runs unchanged on any other single-sample
# GRCh38 VCF plus a list of HPO terms.
#
# | § | Stage | Evidence it adds |
# |---|---|---|
# | 1–3 | Environment, tools (micromamba / bioconda), public reference resources | — |
# | 4 | Phenotype engine: HPO information-content similarity against every annotated Mendelian disease | gene-level phenotype fit |
# | 5 | Case data (gated Hugging Face dataset via the Kaggle Secret `HF_TOKEN`; public demo case without it) | — |
# | 6 | Genome QC, sex, autozygosity, chromosome-level aneuploidy screen with detection limit | data quality, cytogenetic context |
# | 7 | Evidence model: ACMG/AMP Bayesian points, genotype hypotheses, EPCR | — |
# | 8 | Normalise + annotate: gnomAD v4.1, ClinVar, MANE consequences, NMD rule, AlphaMissense | rarity, consequence, prior assertions |
# | 9 | Ranking, second-hit rescue, robustness analysis, submission CSV checked with the organisers' scorer | ranked causal hypotheses |
# | 10 | Deep annotation of the short list (Ensembl VEP, exact gnomAD, PanelApp, chance-bystander model) | refined ranking |
# | 11 | GPU models: SpliceAI on the leading genes, ESM-2 saturation mutagenesis, AlphaFold context | orthogonal evidence |
# | 12 | Phase of the leading pair | cis / trans |
# | 13 | Spike-in benchmark on a public genome | generalisation |
# | 14 | Track 2: mechanism → approved-drug hypotheses, citations resolved live | repurposing evidence |
# | 14b | Track 2: allele-level therapeutic triage (read-through, NMD, exon skipping, cDNA phasing, allele-resolved peptides, ClinVar lever atlas) | repurposing evidence |
# | 15 | Reports (Track 1 methods write-up, Track 2 report) generated from the results | deliverables |
# | 16 | Read-level evidence from the raw FASTQ (unique-minimizer recruitment, re-alignment, re-calling, depth / SV screen) | independent confirmation |
# | 17 | Final reports, deletion of all scratch data, run log | — |
#
# **How to run.** Settings → Accelerator **GPU T4 x2**, Internet **on**. Add-ons → Secrets → add `HF_TOKEN` (a Hugging Face
# token of an account that has been granted access to `SageBio/mva-hackathon-2026-data`) and tick it for this notebook.
# Then **Save Version → Save & Run All** (runs in the background; ~1 h for everything except the FASTQ module, which adds
# ~2–5 h and stops itself after `READ_LEVEL_MAX_HOURS`). Without the secret the notebook runs end-to-end in **demo mode**
# on public data. Outputs (download from the Output tab): `<team>_<model>.csv` (Track 1 submission),
# `<team>_track1_report.md`, `<team>_track2_report.md`, `figures/`, `run_summary.json`.
#
# **Data handling.** The genome is controlled-access. Everything that carries genotype-scale data is written to ephemeral
# scratch storage (`/kaggle/tmp`) and deleted at the end of the run; only findings (ranked candidates, figures, reports)
# are written to `/kaggle/working`. Genes outside the stated indication are not named in displayed tables or reports.
# Keep this notebook **private**.

# %%
CFG = dict(
    # ---- identity (used in output file names; see the submission instructions) -------------------------------------
    TEAM="rusheel86",
    MODEL_NAME="duorank-v1",
    # ---- case inputs ------------------------------------------------------------------------------------------------
    HF_DATASET="SageBio/mva-hackathon-2026-data",
    VCF_FILE="WGS_EX2312012_HGWCNDSX7.vcf.gz",
    PHENOTYPE_DOCX="Challenge_Clinical_Phenotype_1.docx",
    FASTQ_PATTERN="WGS_EX2312012_HGWCNDSX7_S16_L00{lane}_R{read}_001.fastq.gz",
    FASTQ_LANES=[1, 2, 3, 4],
    PROBAND_ID="PROBAND01",
    # Clinical indication as stated by the organisers. It is used as a transparent prior on top of the HPO terms
    # (both rankings are reported). Set to "" for an undiagnosed case.
    CLINICAL_DIAGNOSIS="mosaic variegated aneuploidy",
    HPO_TERMS=[],                 # leave empty to parse them from the phenotype document
    # ---- model parameters -------------------------------------------------------------------------------------------
    MAX_AF=0.01,                  # population allele-frequency ceiling for candidate variants
    FLANK_SECOND_HIT=5000,        # bp around a gene searched for non-coding second hits
    TOP_GENES_DEEP=30,            # genes taken forward to deep annotation / second-hit rescue
    # ---- optional modules -------------------------------------------------------------------------------------------
    RUN_GPU_MODELS=True,
    RUN_READ_LEVEL=True,          # FASTQ module (the long one)
    RUN_BENCHMARK=True,
    RUN_TRACK2=True,
    READ_LEVEL_MAX_HOURS=7.0,     # the FASTQ module stops recruiting new lanes after this many hours
    INCLUDE_SECONDARY_ROWS=False, # secondary / incidental rows in the CSV (off: organisers' IRB guidance is pending)
    MASK_INCIDENTAL=True,         # withhold gene names outside the indication panel / leading hypothesis from tables and reports
    REPO_URL="https://github.com/Rusheel86/MVA",   # public after the hackathon closes, as the rules require
    # ---- disclosure (required by the methods template) -------------------------------------------------------------
    LLM_DISCLOSURE=("Anthropic Claude Code (Claude Opus 5.5, Claude Opus 5 and Claude Sonnet 5.5), Claude Pro plan, was used to design and "
                    "write this pipeline, and to draft the text that it renders into this report. Every literature claim "
                    "was checked against the cited record before it was included. The model-training opt-out was enabled "
                    "and no output ratings or feedback were submitted."),
)

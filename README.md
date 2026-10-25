# DuoRank: phenotype-driven compound-heterozygote discovery and allele-level drug triage (MVA Hackathon 2026)

Code for a fully automated, public-data-only pipeline that ranks genotype hypotheses (homozygous, compound-heterozygous,
heterozygous, hemizygous) from a single-sample GRCh38 VCF plus HPO terms, validates the leading hypothesis with orthogonal
evidence (VEP, gnomAD v4.1, SpliceAI, ESM-2, AlphaFold, read-level re-alignment), and derives approved-drug hypotheses from
the disrupted mechanism (Track 2).

The Track 2 idea is to **treat the alleles** when no approved drug fits the protein. Each causal allele is triaged from sequence,
gene model, protein language model and predicted structure (stop-codon identity, decay prediction, the amino acids read-through
would insert, exon frame, structure), and a single RT-PCR amplicon spanning both alleles is designed automatically. The triage is
checked against published alleles in genes outside the disorder (`p14c`), and repeated over every ClinVar pathogenic or likely
pathogenic allele in the disorder's gene panel (`lever_atlas.csv`).

No patient data is stored in this repository.

## Layout

| Path | Content |
|---|---|
| `src/p00_intro.py` ... `src/p17_finish.py` | notebook sections as percent-format Python (`# %%` cells), one file per section |
| `src/p14_track2.py`, `src/p14b_levers.py`, `src/p14c_lever_validation.py` | Track 2: literature layer and network pass, allele-level triage, validation on published alleles |
| `build_nb.py` | assembles `src/*.py` into the Kaggle notebook (syntax-checks every cell) |
| `kaggle/` | the built notebook and its Kaggle settings |
| `results/` | the reports, figures, `lever_atlas.csv`, `lever_validation.csv` and the submitted prediction file from the run on the challenge data |
| `submission/` | notes for the Track 2 judges |

## Reproduce

1. `python build_nb.py` (needs `nbformat`).
2. `kaggle kernels push -p kaggle`, or upload the `.ipynb` to Kaggle.
3. In the Kaggle editor: Settings, Accelerator **None** (CPU; about 75 minutes, uses no accelerator quota) or GPU T4 x2 / TPU v5e-8,
   Internet **on**; Add-ons, Secrets, `HF_TOKEN` (Hugging Face token of an account with access to `SageBio/mva-hackathon-2026-data`),
   ticked for the notebook.
4. **Save Version, Save & Run All.** Without the secret the notebook runs end to end in demo mode on public data
   (GIAB HG002 with spiked-in ClinVar alleles); a run started by `kaggle kernels push` never sees the secret.

The reports are written by the pipeline itself, so rerunning the notebook regenerates every number, table and figure in them.

Outputs (Kaggle Output tab): `<team>_<model>.csv` (Track 1 submission, validated with the organisers' scoring code),
`<team>_track1_report.md`, `<team>_track2_report.md`, `lever_atlas.csv`, `lever_validation.csv`, `track2_network_drugs.csv`,
`figures/`, `run_summary.json`.

## Reusing the triage on another gene or disorder

`lever_atlas.csv` has one row per public ClinVar pathogenic or likely pathogenic allele: gene, coordinates, protein change, exon,
consequence, ClinVar significance and review status, the lever assigned, the reason, and the stop codon. To build the same table
for another gene list, run the notebook sections through `p14b` and call `lever_atlas([...gene symbols...], fa)` in the same
namespace; `classify_clinvar` assigns the lever for a single allele. Nothing in it is specific to MVA.

## Data handling

Genotype-scale files are written only to ephemeral scratch (`/kaggle/tmp`) and deleted at the end of the run; only findings are
kept. Only short lists of candidate variants are sent to public annotation services. Gene names outside the stated clinical
indication are withheld from tables and reports (`CFG["MASK_INCIDENTAL"]`).

## Public resources

gnomAD v4.1 (Zenodo 22297079 echtvar archive; GCS sites VCFs), ClinVar, MANE v1.5, GRCh38 no-alt analysis set,
AlphaMissense (Zenodo 10813168), HPO, Ensembl VEP REST, Genomics England PanelApp, Open Targets Platform, STRING,
UniProt, AlphaFold DB, Europe PMC, SpliceAI (installed at run time; models CC BY-NC 4.0), ESM-2 (Hugging Face),
GIAB HG002 (benchmark).

## Licence

Code, notebook, reports and figures: CC BY 4.0, as required for submissions to the hackathon.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

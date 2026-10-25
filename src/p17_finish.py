# %% [markdown]
# ## 17 · Final reports, data deletion and run log
#
# The reports are regenerated with every result of the run. All scratch storage that held patient data (VCF, reads,
# annotated variant tables, the re-aligned reads) is deleted; `/kaggle/working` keeps findings only.

# %%
with stage("final reports"):
    log(f"reports written: {write_reports()}")

with stage("delete scratch data"):
    if DEMO:
        log("demo mode: scratch holds public data only; deleting anyway")
    os.chdir(OUT)
    shutil.rmtree(WORK, ignore_errors=True)
    S["scratch_deleted"] = not WORK.exists()
    log(f"scratch directory deleted: {S['scratch_deleted']}")

print("\n" + "=" * 100)
print(f"Mode: {S.get('mode')} · total run time {(time.time() - T0) / 60:.0f} min")
print("=" * 100)
for k, v in RUN.items():
    print(f"{v['status']:>8}  {v['seconds'] / 60:6.1f} min  {k}" + (f"  — {v['note'][:120]}" if v["note"] else ""))
print("\nFiles in the output directory (download these for the submission):")
for p in sorted(OUT.rglob("*")):
    if p.is_file():
        print(f"  {p.relative_to(OUT)}  ({p.stat().st_size / 1024:.0f} KB)")
if SUBMISSION_CSV.exists():
    print(f"\nTrack 1 submission file: {SUBMISSION_CSV.name}")
    print(SUBMISSION_CSV.read_text())
if DEMO:
    print("DEMO run finished. To analyse the challenge genome: Add-ons → Secrets → add HF_TOKEN (a token of a Hugging Face account "
          "that has been granted access to the gated dataset), tick it for this notebook, then Save Version → Save & Run All.")

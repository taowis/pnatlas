# Decisions

Each non-obvious choice, with its reason. Newest last.

## 2026-11-02 — Stage 0

- **Scope: peripheral nervous system first, MASLD/MASH on hold.** A literature
  and data check on 2026-11-02 found each MASH layer already covered by
  2025-2026 papers, the largest staged single-cell cohorts under controlled
  access and only about 33 open spatial MASH biopsies. For ganglia and nerve,
  harmonised DRG atlases exist (Bhuiyan 2024, PMID 38905344; Bhuiyan 2025
  preprint, PMID 41279342), but no joint ganglia + nerve reference, no GTEx
  tibial nerve deconvolution with a nerve reference and no cell-state nerve
  eQTL colocalisation were found. Details in `docs/novelty_check.md`.
- **Separate from PNPrank.** PNPrank is submitted on its own; pnatlas is a
  follow-up that uses PNPrank as an application (stage 6). Merging them would
  hold back a finished manuscript.
- **PNPrank as a dependency, not copied code.** Added at stage 6, pinned to
  the PNPrank release tag used for submission, so the back-test code and
  splits are identical between versions.
- **Reuse processed atlases** from the PNPrank data root (`pnprank_root` in
  `config/paths.yaml`) instead of downloading them again.
- **Nerve layer starts from Heming 2025** (Nat Commun, PMID 40849297; GEO
  GSE285983/GSE285984), building on an earlier local re-analysis.
- **Workflow engine: Snakemake**, same layout as PNPrank (one rules file per
  stage). Packaging with `pyproject.toml`, `src/` layout, dependencies in
  `environment.yml`, `pip install -e . --no-deps`.
- **Environment:** spatial and deconvolution tools (cell2location, squidpy,
  scvi-tools) are added at the stage that needs them, possibly as a separate
  environment because of their torch dependency.
- **Gadi:** heavy integration runs on NCI Gadi. A Snakemake PBS profile is
  added with the first Gadi job, using project qn73 in committed files.
- **Random seed:** `20261102` (scaffold date).
- **Disk budget:** 300 GB under `data_root` (same as PNPrank); ask before
  exceeding.
- **Private notes:** `private_notes/` is kept locally only and is
  gitignored here.
- **Runbook and logs:** every command and manual step is written in
  `docs/runbook/` and the complete output of each step is kept in
  `logs/<stage>/`, so the analysis can be redone by hand. Both are part of the
  public record and go on the public release whitelist. They carry no run dates,
  times or durations; timestamps in tool logs are removed before committing.

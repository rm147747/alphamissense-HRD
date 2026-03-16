# Changelog

All notable changes to this project are documented in this file.

## [1.1.0] — 2026-03-16

### Repository Cleanup (Audit-Ready)
- **Renamed notebooks** for consistent, descriptive naming:
  - `Notebook1_AlphaMissense_PRAD_HRR (3).ipynb` → `Notebook1_Data_Acquisition.ipynb`
  - `Notebook2_Statistical_Analysis.ipynb` → `Notebook2_Concordance_VUS.ipynb`
  - `Notebook3_Validation_mCRPC.ipynb` → `Notebook3_mCRPC_Validation.ipynb`
  - `Notebook4_PanCancer_FINAL.ipynb` → `Notebook4_PanCancer_Survival.ipynb`
- **Removed 3 DEBUG cells** from Notebook 3 (cells 14–16: ad-hoc patient-ID matching diagnostics; not part of the analysis pipeline)
- **Removed stale figures** (pre-correction versions):
  - `Fig_pancancer_KM_pooled.png/pdf` → replaced by corrected version
  - `Fig_pancancer_forest_plot.png/pdf` → replaced by corrected version
- **Renamed `_FIXED` figures** to clean names (the FIX is now the canonical version)
- **Removed `.vscode/settings.json`** (IDE config, not relevant to analysis)
- **Updated README.md** to reflect actual repository contents
- **Added METHODS.md** with detailed methodological documentation
- **Added CHANGELOG.md** (this file)

### Manuscript Correction
- **HR discrepancy fixed**: Manuscript text previously reported HR = 0.87 (p = 0.091), which corresponded to a ridge-penalized pooled Cox sensitivity analysis (Notebook 4, Cell 10). The primary analysis (true stratified Cox, Notebook 4, Cell 6) produces HR = 0.848 (rounded to 0.85), 95% CI 0.716–1.005, p = 0.057. All five instances in the manuscript (Abstract, Section 3.4 ×2, Discussion, Limitations) were corrected.
- **Root cause**: During manuscript drafting, numbers from the sensitivity analysis (ridge λ=0.01) were accidentally substituted for the primary analysis. The `manifest.json` and `pancancer_cox_results.csv` always contained the correct values.
- **Verification**: `manifest.json` reports `stratified_cox_hr: 0.848` — this matches the Notebook 4 output and the corrected manuscript.

## [1.0.0] — 2026-02-16

### Initial Release
- Complete analysis pipeline: 4 notebooks, 25 HRR genes, 31 TCGA tumor types
- ClinVar concordance: κ = 0.733 (bootstrap 95% CI 0.712–0.754)
- VUS reclassification: 90.1% (66,912 / 74,246)
- Pan-cancer stratified Cox: HR = 0.85, 95% CI 0.72–1.01, p = 0.057
- mCRPC feasibility: SU2C/PCF cohort (n = 18, underpowered)
- All data sourced from public databases (TCGA via cBioPortal, ClinVar, AlphaMissense)

# Changelog

## [1.2.0] — 2026-03-16

### New Notebooks
- **Notebook5_Synthetic_Validation.ipynb**: Pipeline unit test (87.3% accuracy, biallelic sens 68%, spec 99%), power analysis (need N=5,000-7,500 for 80% power), adversarial stress test (12/12 scenarios survive p<0.05)
- **Notebook6_Robustness_Analyses.ipynb**: Age-adjusted Cox (HR unchanged at 0.844), Schoenfeld PH test (p=0.870, assumption holds), bootstrap HR (median 0.802, CI 0.669-0.942), permutation test (p<0.001), ridge-penalized Cox sweep

### New Results
- `results/reviewer_robustness_analyses.csv` — 9 robustness tests
- `results/synthetic_module1.csv` — pipeline unit test metrics
- `results/synthetic_module2.csv` — power analysis by sample size
- `results/synthetic_module3.csv` — 12 adversarial stress test scenarios

### Documentation
- Updated README.md to reflect 6 notebooks and all results files

## [1.1.0] — 2026-03-16

### Repository Cleanup (Audit-Ready)
- Renamed notebooks for consistent naming (removed spaces, parentheses, "FINAL")
- Removed 3 DEBUG cells from Notebook 3
- Removed stale pre-correction figure duplicates (_FIXED suffix)
- Removed `.vscode/settings.json`
- Added METHODS.md and CHANGELOG.md

### Manuscript Correction
- HR corrected from 0.87 (sensitivity analysis) to 0.85 (primary analysis) in 5 manuscript instances. Repository code always contained the correct value (0.848). See METHODS.md for details.

## [1.0.0] — 2026-02-16

### Initial Release
- 4 notebooks (NB1-NB4), 25 HRR genes, 31 TCGA tumor types
- ClinVar kappa = 0.733, VUS reclassification 90.1%
- Stratified Cox HR = 0.85, meta-analytic HR = 0.86
- All data from public databases (TCGA, ClinVar, AlphaMissense)

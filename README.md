# AlphaMissense-Guided Reclassification of Variants of Uncertain Significance in Homologous Recombination Repair Genes: A Pan-Cancer Validation Study Across 31 TCGA Tumor Types

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)

## Overview

This repository contains the complete analysis pipeline for the manuscript:

> **AlphaMissense-Guided Reclassification of Variants of Uncertain Significance in Homologous Recombination Repair Genes: A Pan-Cancer Validation Study Across 31 TCGA Tumor Types**
>
> Raphael Brandao, Indianara V. Brandão, Felipe Batalini and Vamsi K. Velcheti

## Key Findings

| Metric | Value | Source |
|--------|-------|--------|
| ClinVar agreement (Cohen's κ) | 0.733 (95% CI 0.712–0.754) | Manuscript Table 2 |
| Sensitivity / Specificity | 75.9% / 95.1% | Manuscript Table 2 |
| VUS reclassification rate | 90.1% (66,912 / 74,246) | Manuscript Section 3.3 |
| Pan-cancer stratified Cox HR | 0.85 (95% CI 0.72–1.01, p = 0.057) | Manuscript Section 3.4 |
| Meta-analytic HR (FE, I² = 0%) | 0.86 (95% CI 0.73–1.02) | Manuscript Section 3.4 |
| Pooled KM log-rank p (unstratified) | 0.0096 | Manuscript Figure 4 |
| Tumor types analyzed | 31 TCGA PanCancer Atlas | Manuscript Section 3.1 |
| Patients with HRR missense | 1,939 (4,301 variants) | Manuscript Section 3.1 |
| Survival cohort | 1,757 patients, 638 events, 20 tumor types | Manuscript Section 3.4 |

> **Note:** The manuscript and pipeline CSVs are now synchronized. The stratified Cox HR = 0.85 (p = 0.057) is the primary analysis; the meta-analytic HR = 0.86 (CI 0.73–1.02) is reported as a robustness check. Minor rounding differences (e.g., 0.848 in CSV → 0.85 in text) are intentional.

## Repository Structure

```
├── Notebook1_AlphaMissense_PRAD_HRR (3).ipynb  # Data acquisition + AM annotation (see note below)
├── Notebook2_Statistical_Analysis.ipynb          # ClinVar concordance + VUS reclassification
├── Notebook3_Validation_mCRPC.ipynb              # SU2C/PCF mCRPC feasibility analysis
├── Notebook4_PanCancer_FINAL.ipynb               # Pan-cancer TCGA analysis (31 tumor types)
├── requirements.txt                              # Python dependencies (pinned versions)
├── manifest.json                                 # Reproducibility manifest (seeds, dates)
├── results/                                      # Analysis outputs
│   ├── pancancer_hrr_variants.csv                # All 4,301 annotated variants (31 studies)
│   ├── pancancer_cox_results.csv                 # Per-tumor Cox results (20 tumor types)
│   ├── annotated_hrr_variants.csv                # PRAD-only variants from Notebook 1 (n=52)
│   ├── patient_hrr_summary.csv                   # Patient-level summary from Notebook 1
│   ├── classification_spec.json                  # AM classification thresholds
│   ├── Table_S_sensitivity.csv                   # Sensitivity analyses (event thresholds, RE)
│   └── Table_S_tumor_imbalance.csv               # AM-pathogenic vs benign by tumor type
├── figures/                                      # Publication-quality figures (PDF + PNG)
│   ├── Fig_pancancer_forest_plot_FIXED.{pdf,png} # Figure 3: Forest plot (manuscript version)
│   ├── Fig_pancancer_KM_pooled_FIXED.{pdf,png}   # Figure 4: Pooled KM (manuscript version)
│   ├── Fig1_AM_distribution.{pdf,png}             # Figure 1: Variant distribution by tumor type
│   ├── Fig2_concordance.{pdf,png}                 # Figure 2: ClinVar concordance matrix
│   └── (additional exploratory figures)
└── LICENSE
```

> **Naming conventions:** Files with `_FIXED` suffix are the final manuscript versions. Files without this suffix are earlier iterations retained for audit trail. Figures numbered `Fig1`–`Fig6` are from earlier Notebook 2 (PRAD-focused); the manuscript figures are the `Fig_pancancer_*_FIXED` files.

### About "PRAD" in Notebook 1

Notebook 1 is named `Notebook1_AlphaMissense_PRAD_HRR` because the project began as a prostate cancer–focused study. This notebook downloads **TCGA-PRAD** data only (52 HRR missense variants in 40 patients) and serves as the pilot analysis. **Notebook 4** subsequently expands the pipeline to all 31 TCGA PanCancer Atlas tumor types (4,301 variants, 1,939 patients). The PRAD-specific output (`results/annotated_hrr_variants.csv`) is an intermediate artifact from the pilot phase and is **not** used in the manuscript's pan-cancer analyses.

## Reproducibility

### Quick Start

```bash
git clone https://github.com/rm147747/alphamissense-prostate-hrr.git
cd alphamissense-prostate-hrr

pip install -r requirements.txt

# Run notebooks in order: NB1 → NB2 → NB3 → NB4
jupyter notebook
```

### Re-run from Scratch (Clean State)

```bash
rm -rf results/ figures/
mkdir -p results figures

# Then execute notebooks in order: NB1 → NB2 → NB3 → NB4
# Each notebook will recreate its outputs from upstream data.
```

> **Note:** NB1 downloads data from cBioPortal API (requires internet). Subsequent notebooks depend on outputs from prior ones — running out of order will fail with a clear FileNotFoundError.

### Environment

- **Python:** 3.12.3
- **Random seed:** 42 (set globally in every notebook)
- **Data access date:** February 16, 2026
- **Key packages:** pandas 2.3.3, lifelines 0.30.1, scikit-learn 1.8.0, scipy 1.17.0

All package versions are pinned in `requirements.txt`. A reproducibility manifest is in `manifest.json`.

### Notebook Execution Order

Notebooks must be run in order (1 → 2 → 3 → 4), as each depends on outputs from the previous:

1. **Notebook 1** (`Notebook1_AlphaMissense_PRAD_HRR (3).ipynb`): Downloads TCGA-PRAD mutations, filters HRR missense variants, annotates with AlphaMissense. Produces PRAD pilot results.
2. **Notebook 2** (`Notebook2_Statistical_Analysis.ipynb`): Compares AlphaMissense with ClinVar, computes concordance metrics (κ = 0.733), reclassifies VUS.
3. **Notebook 3** (`Notebook3_Validation_mCRPC.ipynb`): Applies pipeline to SU2C/PCF mCRPC cohort (exploratory feasibility, n = 18).
4. **Notebook 4** (`Notebook4_PanCancer_FINAL.ipynb`): Pan-cancer analysis across 31 TCGA tumor types. Produces all pan-cancer results, forest plot, and pooled KM curves.

### GitHub Codespaces

This repository is configured for GitHub Codespaces. Click the green "Code" button → "Codespaces" → "Create codespace on main" for a ready-to-use environment.

## Data Sources

| Source | Access | Date |
|--------|--------|------|
| TCGA PanCancer Atlas (31 studies) | cBioPortal API | Feb 2026 |
| AlphaMissense | Public database (Zenodo) | Feb 2026 |
| ClinVar | NCBI FTP | Feb 2026 |
| SU2C/PCF mCRPC | cBioPortal | Feb 2026 |

## HRR Gene Panel (25 genes)

- **Cohort A** (PROfound, established): BRCA1, BRCA2, ATM
- **Cohort B** (PROfound, expanded): PALB2, BRIP1, BARD1, CDK12, CHEK1, CHEK2, FANCL, RAD51B, RAD51C, RAD51D, RAD54L
- **Extended DDR/FA pathway**: FANCA, FANCC, FANCD2, FANCE, FANCF, FANCG, NBN, MRE11, RAD50, ATR, ATRX

## AlphaMissense Classification Rule

AlphaMissense provides **pre-computed categorical labels** (`benign`, `ambiguous`, `pathogenic`) directly from the model output. This pipeline uses these labels as-is — **no custom thresholds are applied**.

| Classification | Score range | Source |
|---------------|------------|--------|
| Benign | am_pathogenicity < 0.34 | Cheng et al. Science 2023 |
| Ambiguous | 0.34 ≤ score ≤ 0.564 | Cheng et al. Science 2023 |
| Pathogenic | score > 0.564 | Cheng et al. Science 2023 |

**Patient-level derivation:** `has_am_pathogenic = True` if ≥1 variant has `am_class == "pathogenic"`.

## Statistical Methods

- **Primary analysis:** Stratified Cox model (strata = tumor type; separate baseline hazards; ridge penalization λ = 0.01)
- **Meta-analysis:** Fixed-effect inverse-variance + REML random-effects + Hartung–Knapp CI + prediction interval
- **Concordance:** Cohen's κ with 2,000-iteration bootstrap 95% CI
- **Sensitivity:** Event threshold (≥3/5/10), ridge penalty sweep (λ = 0.001–0.5), leave-one-tumor-out meta-analysis

## Scope & Ethics

This repository is **paper-only** — it reproduces the analyses described in the manuscript. It is not a supported CLI tool or library.

This study used only publicly available, de-identified datasets (TCGA via cBioPortal, ClinVar, AlphaMissense). IRB approval was not required.

## Limitations & Intended Use

> **⚠️ This analysis is hypothesis-generating and NOT intended for clinical decision-making.**

- AlphaMissense predictions are **computational annotations**, not clinical reclassifications per ACMG/AMP standards.
- The VUS "reclassification" reported here is a **computational triage** — it does not replace expert curation, functional assays, or clinical-grade variant interpretation.
- Survival analysis is exploratory. The proportional hazards assumption was not formally assessed (Schoenfeld residuals). The stratified Cox model controls for tumor type but does not adjust for age, stage, or treatment.
- **Biallelic inactivation / LOH status is not assessed.** The current analysis does not distinguish monoallelic from biallelic HRR variants. This is a major limitation — HRR functional deficiency typically requires biallelic loss, and the observed survival effect is likely diluted by monoallelic passengers. Integration with LOH data and mutational signature 3 (BRCAness) is planned.
- For clinical use, AlphaMissense scores should be considered as **PP3/BP4-level supporting evidence** within the ACMG/AMP framework, not as standalone determinants.
- Prospective clinical validation is required before integration into treatment decisions.

## Known Issues (Pre-Submission)

| Issue | Status | Resolution |
|-------|--------|------------|
| ~~Stratified Cox HR differed between CSV and manuscript~~ | **Resolved** | Manuscript corrected to match pipeline primary analysis (HR = 0.85, p = 0.057) |
| PH assumption not tested via Schoenfeld residuals | Identified | Will add in final version |
| Notebook 1 named "PRAD" but paper is pan-cancer | By design | NB1 = pilot; NB4 = pan-cancer expansion (see note above) |
| Some intermediate files referenced in prior README do not exist in repo | Identified | Non-critical; these are generated at runtime by notebooks |
| Figures contain both `_FIXED` and non-FIXED versions | By design | `_FIXED` = manuscript versions; others retained for audit |

## Citation

If you use this code or data, please cite:

```
Moreira, R.B.; Brandão, I.V.; Velcheti, V.K. AlphaMissense-Guided Reclassification of
Variants of Uncertain Significance in Homologous Recombination Repair Genes: A Pan-Cancer
Validation Study Across 31 TCGA Tumor Types. [Journal] 2026.
```

## License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

## Contact

- **Raphael B. Moreira, MD** — rb@firstsaude.com — São Camilo Hospital Network / UNIFESP, São Paulo
- **Vamsi K. Velcheti, MD** — Mayo Clinic, Jacksonville, FL

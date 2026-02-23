# AlphaMissense-Guided Reclassification of HRR Variants of Uncertain Significance: A Pan-Cancer Benchmarking Study

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)

## Overview

This repository contains the complete analysis pipeline for the manuscript:

> **AlphaMissense-Guided Reclassification of Variants of Uncertain Significance in Homologous Recombination Repair Genes: A Pan-Cancer Benchmarking Study**
>
> Raphael B. Moreira and Vamsi K. Velcheti
>


## Key Findings

| Metric | Value |
|--------|-------|
| ClinVar agreement (Cohen's κ) | 0.733 (95% bootstrap CI 0.712–0.754) |
| VUS reclassification rate | 90.1% (66,912 / 74,246) |
| Pan-cancer stratified Cox HR | 0.85 (95% CI 0.72–1.01, p = 0.057) |
| Meta-analytic HR (FE, I² = 0%) | 0.86 (95% CI 0.73–1.02) |
| Tumor types analyzed | 31 TCGA PanCancer Atlas |
| Patients with HRR missense | 1,939 |

## Repository Structure

```
├── Notebook1_AlphaMissense_PRAD_HRR.ipynb   # Data acquisition + annotation
├── Notebook2_Concordance_VUS.ipynb           # ClinVar concordance + VUS reclassification
├── Notebook3_mCRPC_Validation.ipynb          # SU2C/PCF mCRPC feasibility analysis
├── Notebook4_PanCancer_v2_FIX.ipynb          # Pan-cancer TCGA analysis (31 tumor types)
├── requirements.txt                          # Python dependencies (pinned versions)
├── manifest.json                             # Reproducibility manifest (seeds, dates, checksums)
├── data/                                     # Generated data (created by notebooks)
│   ├── processed/                            # AlphaMissense lookup tables
│   └── pancancer/                            # Per-study intermediate files
├── results/                                  # Analysis outputs
│   ├── pancancer_hrr_variants.csv            # All 4,301 annotated variants
│   ├── pancancer_cox_results.csv             # Per-tumor Cox results
│   ├── Table_S_sensitivity.csv               # Sensitivity analyses
│   └── Table_S_tumor_imbalance.csv           # Tumor-type distribution by AM status
├── figures/                                  # Publication-quality figures
└── LICENSE
```

## Reproducibility

### Quick Start

```bash
# Clone
git clone https://github.com/rm147747/alphamissense-prostate-hrr.git
cd alphamissense-prostate-hrr

# Install dependencies
pip install -r requirements.txt

# Run notebooks in order
jupyter notebook Notebook1_AlphaMissense_PRAD_HRR.ipynb
```

### Re-run from Scratch (Clean State)

To ensure full reproducibility without leftover artifacts:

```bash
rm -rf results/ figures/ data/processed/ data/pancancer/
mkdir -p results figures data/processed data/pancancer

# Then execute notebooks in order: NB1 → NB2 → NB3 → NB4
# Each notebook will recreate its outputs from upstream data.
```

> **Note:** NB1 downloads data from cBioPortal (requires internet). Subsequent notebooks depend on outputs from prior ones — running out of order will fail with a clear FileNotFoundError.

### Environment

- **Python:** 3.12.3
- **Random seed:** 42 (set globally in every notebook)
- **Data access date:** February 16, 2026
- **Key packages:** pandas 2.3.3, lifelines 0.30.1, scikit-learn 1.8.0, scipy 1.17.0

All package versions are pinned in `requirements.txt`. A complete reproducibility manifest including data access dates and random seeds is in `manifest.json`.

### Notebook Execution Order

Notebooks must be run in order (1 → 2 → 3 → 4), as each depends on outputs from the previous:

1. **Notebook 1:** Downloads TCGA-PRAD mutations, filters HRR missense variants, annotates with AlphaMissense
2. **Notebook 2:** Compares AlphaMissense with ClinVar, computes concordance metrics, reclassifies VUS
3. **Notebook 3:** Applies pipeline to SU2C/PCF mCRPC cohort (exploratory feasibility)
4. **Notebook 4:** Pan-cancer analysis across 31 TCGA tumor types with survival analysis

### GitHub Codespaces

This repository is configured for GitHub Codespaces. Click the green "Code" button → "Codespaces" → "Create codespace on main" for a ready-to-use environment.

## Data Sources

| Source | Access | Date |
|--------|--------|------|
| TCGA PanCancer Atlas | cBioPortal API | Feb 2026 |
| AlphaMissense | Public database | Feb 2026 |
| ClinVar | NCBI | Feb 2026 |
| SU2C/PCF mCRPC | cBioPortal | Feb 2026 |

## HRR Gene Panel (25 genes)

- **Cohort A** (established): BRCA1, BRCA2, ATM
- **Cohort B** (PROfound): PALB2, BRIP1, BARD1, CDK12, CHEK1, CHEK2, FANCL, RAD51B, RAD51C, RAD51D, RAD54L
- **Extended DDR**: FANCA, FANCC, FANCD2, FANCE, FANCF, FANCG, NBN, MRE11, RAD50, ATR, ATRX

## AlphaMissense Classification Rule

AlphaMissense provides **pre-computed categorical labels** (`benign`, `ambiguous`, `pathogenic`) directly from the model output. This pipeline uses these labels as-is — **no custom thresholds are applied**.

| Classification | Score range | Source |
|---------------|------------|--------|
| Benign | am_pathogenicity < 0.34 | AlphaMissense TSV column 4 |
| Ambiguous | 0.34 ≤ score ≤ 0.564 | AlphaMissense TSV column 4 |
| Pathogenic | score > 0.564 | AlphaMissense TSV column 4 |

The thresholds were calibrated by [Cheng et al. (Science 2023)](https://doi.org/10.1126/science.adg7492) to achieve 90% precision on ClinVar pathogenic variants. A full classification specification is available in `results/classification_spec.json`.

**Patient-level derivation:** `has_am_pathogenic = True` if ≥1 variant has `am_class == "pathogenic"` (defined in Notebook 1, Cell 22).

## Statistical Methods

- **Primary analysis:** True stratified Cox model (strata = tumor type; separate baseline hazards)
- **Meta-analysis:** Fixed-effect inverse-variance + REML random-effects + Hartung–Knapp CI + prediction interval
- **Concordance:** Cohen's κ with 2,000-iteration bootstrap 95% CI
- **Sensitivity:** Event threshold (≥3/5/10), ridge penalty sweep (λ = 0.001–0.5)

## Scope & Ethics

This repository is **paper-only** — it reproduces the analyses described in the manuscript. It is not a supported CLI tool or library.

This study used only publicly available, de-identified datasets (TCGA via cBioPortal, ClinVar, AlphaMissense). IRB approval was not required.

## Limitations & Intended Use

> **⚠️ This analysis is hypothesis-generating and NOT intended for clinical decision-making.**

- AlphaMissense predictions are **computational annotations**, not clinical reclassifications per ACMG/AMP standards.
- The VUS "reclassification" reported here is a **computational triage** — it does not replace expert curation, functional assays, or clinical-grade variant interpretation.
- Survival analysis is **univariate** (no adjustment for age, stage, treatment, or other confounders) and should be interpreted as associative, not causal.
- For clinical use, AlphaMissense scores should be considered as **PP3/BP4-level supporting evidence** within the ACMG/AMP framework, not as standalone determinants.
- Prospective clinical validation is required before integration into treatment decisions.

## Intermediate Artifacts

Notebooks must be run in order. Key intermediate files:

| File | Produced by | Consumed by |
|------|------------|-------------|
| `results/annotated_hrr_variants.csv` | Notebook 1 | Notebook 2 |
| `results/patient_hrr_summary.csv` | Notebook 1 | Notebook 2 |
| `data/processed/alphamissense_hrr_genes.csv` | Notebook 1 | Notebooks 3, 4 |
| `results/concordance_data.csv` | Notebook 2 | Notebook 4 |
| `results/concordance_results.csv` | Notebook 2 | — |
| `results/vus_reclassification.csv` | Notebook 2 | — |
| `results/pancancer_hrr_variants.csv` | Notebook 4 | — |
| `results/pancancer_cox_results.csv` | Notebook 4 | — |

## Citation

If you use this code or data, please cite:

```
Moreira, R.B.; Velcheti, V.K. AlphaMissense-Guided Reclassification of Variants of
Uncertain Significance in Homologous Recombination Repair Genes: A Pan-Cancer
Benchmarking Study. Diagnostics 2026, X, XX. https://doi.org/10.3390/diagnosticsXXXXXXXX
```

## License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/). If you use this code or data, please cite the associated publication.

## Contact

- **Raphael B. Moreira** — rb@firstsaude.com — São Camilo Hospital
- **Vamsi K. Velcheti** — Mayo Clinic, Jacksonville, FL

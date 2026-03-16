# AlphaHRD: Multi-Layer Pipeline for HRR VUS Functional Classification

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Overview

> **AlphaHRD: Integrating AlphaMissense, Loss of Heterozygosity, and Genomic Instability Scores for Functional Classification of HRR Variants of Uncertain Significance Across 31 Cancer Types**
>
> Raphael B. Moreira, Indianara V. Brandao, Felipe Batalini and Vamsi K. Velcheti

## Key Findings

| Metric | Value |
|--------|-------|
| ClinVar kappa | 0.733 (CI 0.712-0.754) |
| VUS reclassified | 90.1% (66,912/74,246) |
| Stratified Cox HR | 0.85 (CI 0.72-1.01, p=0.057) |
| Biallelic vs Mono HRDsum | 28.0 vs 12.4 (p=6.8e-20) |
| HRD+ biallelic vs mono | 24.8% vs 0.0% (p=9.9e-24) |
| TRUE_HRD identified | 38 patients |
| Tumors analyzed | 31 TCGA PanCancer Atlas |

## Notebooks

| NB | Purpose | Output |
|----|---------|--------|
| NB1 Data Acquisition | TCGA download, HRR filter, AM+ClinVar | annotated_hrr_variants.csv |
| NB2 Concordance VUS | ClinVar kappa, VUS, PRAD survival | Figures 1-5 |
| NB3 mCRPC Validation | SU2C/PCF feasibility (n=18) | Figure 6 |
| NB4 PanCancer Survival | 31-tumor Cox, meta-analysis | pancancer_hrr_variants.csv |
| NB5 Synthetic Validation | Unit test 87%, power, 12/12 stress | synthetic_module1-3.csv |
| NB6 Robustness | Schoenfeld p=0.87, permutation p<0.001 | reviewer_robustness.csv |

## Documentation

- [METHODS.md](METHODS.md) - Detailed methodology for audit
- [CHANGELOG.md](CHANGELOG.md) - Version history and corrections
- manifest.json - Reproducibility manifest (seed=42, dates)

## Reproducibility

Python 3.12.3, seed 42, data Feb 16 2026. Run NB1-NB6 in order.
Codespaces ready. All deps in requirements.txt.

## Limitations

Hypothesis-generating, NOT for clinical decisions. Univariate survival
(age-adjusted HR unchanged). PH verified (Schoenfeld p=0.87). Bootstrap
CI 0.669-0.942. Permutation p<0.001. Needs PARPi-treated validation.

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Contact

- Raphael B. Moreira, MD - rb@firstsaude.com
- Vamsi K. Velcheti, MD - Mayo Clinic, Jacksonville, FL

# AlphaHRD: Multi-Layer Pipeline for HRR VUS Functional Classification

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Overview

> **AlphaHRD: Integrating AlphaMissense, Loss of Heterozygosity, and Genomic Instability Scores for Functional Classification of HRR Variants of Uncertain Significance Across 31 Cancer Types**
>
> Raphael B. Moreira, Indianara V. Brandao, Felipe Batalini and Vamsi K. Velcheti

A computational pipeline that integrates AlphaMissense pathogenicity scores, allele-specific LOH from ASCAT, genomic scar scores (HRDsum), structural validation (Boltz-2), BRCA1 promoter methylation, and immune microenvironment features to functionally classify HRR variants of uncertain significance across 31 TCGA pan-cancer types.

## Key Findings

| Metric | Value |
|--------|-------|
| ClinVar concordance (kappa) | 0.733 (CI 0.712-0.754) |
| VUS reclassified | 90.1% (66,912/74,246) |
| Unadjusted Cox HR | 0.803 (CI 0.683-0.945, p=0.008) |
| Germline carriers HR | 0.617 (CI 0.473-0.806, p<0.001) |
| RMST delta (tau=81mo) | +5.2 months (CI 1.5-8.5, p=0.002) |
| Biallelic vs Mono HRDsum | 28.0 vs 12.4 (p=6.8e-20) |
| HRD+ biallelic vs mono | 24.8% vs 0.0% (p=9.9e-24) |
| TRUE_HRD identified | 38 patients |
| Purity-filtered HR (>=0.5) | 0.699 (CI 0.573-0.853, p<0.001) |
| E-value | 1.80 (CI bound 1.31) |
| Tumors analyzed | 31 TCGA PanCancer Atlas |

## Notebooks (execution order)

| # | Script | Purpose | Key Output |
|---|--------|---------|------------|
| NB1 | Data Acquisition | TCGA download, HRR filter, AM+ClinVar | annotated_hrr_variants.csv |
| NB2 | Concordance VUS | ClinVar kappa, VUS reclassification | Figures 1-5 |
| NB3 | mCRPC Validation | SU2C/PCF feasibility (n=18) | Figure 6 |
| NB4 | PanCancer Survival | 31-tumor Cox, LOH, HRDsum, Boltz-2 | pancancer_hrr_variants.csv |
| NB5 | Synthetic Validation | Unit test 87%, power, 12/12 stress | synthetic_module1-3.csv |
| NB6 | Robustness | Schoenfeld, bootstrap, permutation, ridge | reviewer_robustness.csv |
| **NB7** | **Extended Robustness** | **Germline/somatic, stage-adj, RMST, purity, E-value** | **robustness_5analyses_FINAL.csv** |
| **NB8** | **BRCA1 Methylation** | **Promoter methylation as independent HRD mechanism** | **analysis6_results.csv** |
| **NB9** | **TMB + Immune + FoldX** | **Thorsson 2018 integration, FoldX preparation** | **immune_features_comparison.csv** |

## 8 Robustness Analyses

| # | Analysis | Result | Status |
|---|----------|--------|--------|
| 1 | Germline vs Somatic | Germline HR=0.62 (p<0.001); Somatic HR=0.91 (ns) | Done |
| 2 | Stage-adjusted Cox | HR attenuates to 0.87 (p=0.16) | Done |
| 3 | RMST | delta=+3.1 to +7.3 months, all p<0.01 | Done |
| 4 | Purity filter | HR strengthens: 0.80 to 0.70 | Done |
| 5 | E-value | E=1.80 (CI bound=1.31) | Done |
| 6 | BRCA1 methylation | Independent of AM status (OR~1.0) | Done |
| 7 | FoldX ddG | Script ready, requires local FoldX | Pending |
| 8 | TMB + Immune | TMB 3x higher; immune infiltrate identical | Done |

## Documentation

- [METHODS.md](METHODS.md) - Detailed methodology for every analysis
- [CHANGELOG.md](CHANGELOG.md) - Version history and corrections
- [AUDIT.md](AUDIT.md) - Complete audit trail with data provenance and AI disclosure
- manifest.json - Reproducibility manifest with key results

## Reproducibility

Python 3.12.3, seed 42. All deps in requirements.txt. Data access: Feb 16 + Mar 16-17, 2026.
Run NB1-NB9 in order. NB7-NB9 depend on NB4 outputs. Codespaces ready.

## Limitations

Hypothesis-generating, NOT for clinical decisions. Stage adjustment attenuates HR.
E-value=1.80 is moderate. TCGA is treatment-naive. Needs PARPi-treated validation.

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Contact

- Raphael B. Moreira, MD - rb@firstsaude.com
- Vamsi K. Velcheti, MD - Mayo Clinic, Jacksonville, FL

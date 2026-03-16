# Detailed Methods

This document describes every analytical decision for audit and peer review.

## 1. Study Design

Retrospective pan-cancer computational benchmarking study. Not pre-registered. All analyses are hypothesis-generating.

## 2. Data Sources

### 2.1 TCGA PanCancer Atlas
- **API:** `https://www.cbioportal.org/api` (Feb 16, 2026)
- **Studies:** 31 TCGA PanCancer Atlas (`_tcga_pan_can_atlas_2018`)
- **Code:** NB1 Cell 6, NB4 Cell 3

### 2.2 AlphaMissense
- **Source:** Per-protein TSV from `https://alphamissense.hegelab.org`
- **Publication:** Cheng et al., Science 2023, doi:10.1126/science.adg7492
- **Classification:** Pre-computed labels (benign < 0.34, ambiguous 0.34-0.564, pathogenic > 0.564). No custom thresholds applied.
- **Code:** NB1 Cell 12

### 2.3 ClinVar
- **Source:** NCBI E-utilities API (Feb 16, 2026)
- **Inclusion:** Pathogenic/Likely pathogenic and Benign/Likely benign only
- **Code:** NB1 Cell 14

### 2.4 SU2C/PCF mCRPC
- **Study:** `prad_su2c_2019` (Abida et al., PNAS 2019)
- **Code:** NB3 Cell 3

## 3. Gene Panel (25 genes)

Cohort A (established): BRCA1, BRCA2, ATM. Cohort B (PROfound): PALB2, BRIP1, BARD1, CDK12, CHEK1, CHEK2, FANCL, RAD51B, RAD51C, RAD51D, RAD54L. Extended DDR: FANCA, FANCC, FANCD2, FANCE, FANCF, FANCG, NBN, MRE11, RAD50, ATR, ATRX.

## 4. Variant Filtering

Missense mutations only (`Variant_Classification == "Missense_Mutation"`). Non-missense excluded because AlphaMissense is designed specifically for missense classification.

## 5. Patient Classification

`has_am_pathogenic = True` if patient carries >= 1 variant with `am_class == "pathogenic"`.

## 6. Concordance Analysis

Cohen's kappa with 2,000-iteration bootstrap 95% CI. AlphaMissense pathogenic aligned with ClinVar Pathogenic/Likely pathogenic. Ambiguous variants excluded from binary concordance.

## 7. Survival Analysis

### 7.1 Primary
Stratified Cox PH model: `OS ~ has_am_pathogenic, strata = tumor_type`. Univariate. No adjustment for age, stage, treatment.

### 7.2 Meta-Analysis
Per-tumor log(HR) combined via fixed-effect, REML random-effects, Hartung-Knapp CI, and prediction interval. Tumor-specific estimates are NOT independent studies — meta-analysis is treated as robustness aggregation.

### 7.3 Sensitivity
Event threshold (>= 3/5/10), ridge penalty (lambda = 0.001-0.5), leave-one-gene-out.

## 8. Robustness Analyses (Notebook 6)

Pre-addressing reviewer concerns with 5 tests:

### 8.1 Age-Adjusted Cox
Model: `OS ~ has_am_pathogenic + age_scaled, strata = tumor`. Age pulled from cBioPortal clinical data (98% available). Result: HR changes by < 0.1%, indicating minimal confounding by age.

### 8.2 Schoenfeld Residuals
Tests proportional hazards assumption. Uses `lifelines.statistics.proportional_hazard_test`. Result: p = 0.870 (assumption holds).

### 8.3 Bootstrap HR
500 resamples with replacement. Simple Cox model refitted each time. Provides non-parametric 95% CI. Result: median HR = 0.802, CI 0.669-0.942 (excludes 1.0).

### 8.4 Permutation Test
500 random shuffles of AM-pathogenic labels. Tests null hypothesis that label assignment is random. Result: p < 0.001 (0/500 permutations produced HR as extreme as observed).

### 8.5 Ridge-Penalized Cox
L2 penalty sweep (lambda = 0.001, 0.01, 0.1, 0.5) on stratified Cox. Tests robustness to small-sample bias. Result: HR stable at 0.84-0.85 for low penalty; expected attenuation toward null at high penalty.

## 9. Synthetic Validation (Notebook 5)

### 9.1 Pipeline Unit Test
5,000 synthetic patients with known ground truth. Simulated AM errors (sensitivity 76%, specificity 95%) and LOH errors (sensitivity 90%, FPR 5%). Classification accuracy: 87.3%. Biallelic detection: sensitivity 68%, specificity 99%. HRDsum enrichment recovered at p = 1.85e-21.

### 9.2 Power Analysis
Simulated survival studies at N = 500-10,000. Assumptions: TRUE_HRD prevalence 2%, HR = 0.70, median OS control 60 months. Result: current TCGA (N ~ 2,000) has ~40% power; minimum N for 80% power = 5,000-7,500 patients.

### 9.3 Adversarial Stress Test
12 scenarios with systematically injected errors (AM FPR 5-15%, LOH FPR 10-30%, combined worst cases). Core finding (biallelic HRD enrichment) survives all 12/12 scenarios at p < 0.05, including extreme scenario with 50% AM error + 30% LOH error (p = 0.005).

## 10. Role of AI Tools

### What AI assisted with
- Code generation, manuscript formatting, literature search
- Pipeline automation and figure generation

### What was human-led
- Study design, gene panel selection, statistical method choices
- Interpretation of results, clinical context
- All data acquisition decisions
- Review and validation of all outputs

### Verification
Every result traceable to a specific notebook cell. `manifest.json` contains key results for cross-verification. All notebooks can be re-run from scratch.

## 11. Reproducibility

- Python 3.12.3, random seed 42
- Packages pinned in `requirements.txt`
- Data access date: Feb 16, 2026
- Notebooks run in order: NB1 -> NB2 -> NB3 -> NB4 -> NB5 -> NB6
- NB5 and NB6 depend on NB4 outputs

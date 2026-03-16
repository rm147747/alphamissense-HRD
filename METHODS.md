# Detailed Methods

This document provides a granular description of every analytical decision in this study. It is intended to serve as an audit trail demonstrating that all results derive from the code in this repository, executed on publicly available data, with pre-specified methods.

## 1. Study Design

**Type:** Retrospective pan-cancer computational benchmarking study.

**Objective:** Evaluate AlphaMissense as a triage tool for variants of uncertain significance (VUS) in 25 homologous recombination repair (HRR) genes, using ClinVar as the reference standard and TCGA overall survival as an exploratory clinical endpoint.

**Registration:** This study was not pre-registered. All analyses are considered hypothesis-generating.

## 2. Data Sources

### 2.1 TCGA PanCancer Atlas (via cBioPortal)

- **API endpoint:** `https://www.cbioportal.org/api`
- **Access date:** February 16, 2026
- **Studies:** All 31 TCGA PanCancer Atlas studies (suffix `_tcga_pan_can_atlas_2018`)
- **Data types:** Somatic mutations (MAF format), clinical data (OS_MONTHS, OS_STATUS)
- **Code:** `Notebook1_Data_Acquisition.ipynb`, Cell 6 (single-study PRAD); `Notebook4_PanCancer_Survival.ipynb`, Cell 3 (all 31 studies)
- **Reproducibility note:** cBioPortal API returns live data. Results may differ if TCGA annotations are updated. The `manifest.json` records the access date. Intermediate CSVs in `results/` serve as frozen snapshots.

### 2.2 AlphaMissense

- **Source:** Pre-computed per-protein TSV files from `https://alphamissense.hegelab.org`
- **Publication:** Cheng et al., Science 2023, doi:10.1126/science.adg7492
- **Genes loaded:** 25 HRR genes (see Section 3 below)
- **Variants loaded:** 554,363 missense variant predictions across all 25 genes
- **Code:** `Notebook1_Data_Acquisition.ipynb`, Cell 12
- **Classification:** We use the categorical labels (`benign`, `ambiguous`, `pathogenic`) as provided by AlphaMissense. No custom score thresholds are applied. The thresholds (benign < 0.34, ambiguous 0.34–0.564, pathogenic > 0.564) were calibrated by the original authors.

### 2.3 ClinVar

- **Source:** NCBI ClinVar database via Entrez E-utilities API
- **Access date:** February 16, 2026
- **Inclusion:** Variants with definitive clinical significance (`Pathogenic`, `Likely pathogenic`, `Benign`, `Likely benign`) in any of the 25 HRR genes
- **Exclusion:** VUS, conflicting interpretations, and variants with only uncertain significance
- **Code:** `Notebook1_Data_Acquisition.ipynb`, Cell 14
- **Matching:** Variants are matched to AlphaMissense predictions by gene + amino acid substitution (e.g., BRCA1 p.C61G)

### 2.4 SU2C/PCF mCRPC

- **cBioPortal study ID:** `prad_su2c_2019`
- **Publication:** Abida et al., PNAS 2019, doi:10.1073/pnas.1902651116
- **Code:** `Notebook3_mCRPC_Validation.ipynb`, Cells 3–4

## 3. HRR Gene Panel

25 genes selected based on clinical relevance to PARP inhibitor trials:

| Tier | Genes | Rationale |
|------|-------|-----------|
| Cohort A (established) | BRCA1, BRCA2, ATM | PROfound Cohort A, strongest PARPi evidence |
| Cohort B (PROfound) | PALB2, BRIP1, BARD1, CDK12, CHEK1, CHEK2, FANCL, RAD51B, RAD51C, RAD51D, RAD54L | PROfound Cohort B |
| Extended DDR | FANCA, FANCC, FANCD2, FANCE, FANCF, FANCG, NBN, MRE11, RAD50, ATR, ATRX | Broader DDR panel used in TRITON3, TALAPRO-2, and other trials |

**Code:** `Notebook1_Data_Acquisition.ipynb`, Cell 4

## 4. Variant Filtering

### 4.1 Mutation Selection
- **Included:** Missense mutations only (`Variant_Classification == "Missense_Mutation"`)
- **Excluded:** Frameshift, nonsense, splice site, in-frame indels, silent mutations
- **Rationale:** AlphaMissense is designed specifically for missense variant classification. Non-missense variants (e.g., truncating mutations) have well-established pathogenicity rules and are not the target of this analysis.
- **Code:** `Notebook1_Data_Acquisition.ipynb`, Cell 8

### 4.2 Protein Change Parsing
- HGVS protein notation (e.g., `p.V600E`) is parsed to extract: reference amino acid, position, alternate amino acid
- Three-letter amino acid codes are converted to one-letter codes
- **Code:** `Notebook1_Data_Acquisition.ipynb`, Cell 10

## 5. AlphaMissense Annotation

### 5.1 Matching Strategy
Each TCGA variant is matched to the AlphaMissense database by constructing a key: `{UniProt_ID}_{position}_{ref_aa}_{alt_aa}`. This ensures exact protein-level matching.

### 5.2 Classification Rule
AlphaMissense provides pre-computed categorical labels. We use these directly:

| Label | Score Range | Interpretation |
|-------|-----------|----------------|
| Pathogenic | > 0.564 | Predicted deleterious by AlphaMissense |
| Ambiguous | 0.34 – 0.564 | Uncertain prediction |
| Benign | < 0.34 | Predicted benign by AlphaMissense |

**No custom thresholds are applied.** The `classification_spec.json` documents this decision.

### 5.3 Patient-Level Derivation
A patient is classified as `has_am_pathogenic = True` if they carry ≥1 missense variant in any HRR gene with `am_class == "pathogenic"`. This is a binary grouping variable for survival analysis.

**Code:** `Notebook1_Data_Acquisition.ipynb`, Cell 22

## 6. Concordance Analysis

### 6.1 Reference Standard
ClinVar definitive classifications (Pathogenic/Likely pathogenic vs. Benign/Likely benign). Variants with uncertain significance or conflicting interpretations are excluded from the concordance analysis.

### 6.2 Metrics
- **Cohen's κ** with 2,000-iteration bootstrap 95% CI (non-parametric percentile method)
- **Sensitivity, specificity, PPV, NPV** computed from the 2×2 confusion matrix
- AlphaMissense `pathogenic` aligned with ClinVar `Pathogenic/Likely pathogenic`
- AlphaMissense `benign` aligned with ClinVar `Benign/Likely benign`
- AlphaMissense `ambiguous` excluded from binary concordance (reported separately)

### 6.3 VUS Reclassification
All ClinVar VUS that have an AlphaMissense prediction (score available) are counted. The "reclassification" is computational triage — it assigns an AlphaMissense label to previously unclassified variants. This is NOT a clinical reclassification per ACMG/AMP standards.

**Code:** `Notebook2_Concordance_VUS.ipynb`, Cells 9–12

## 7. Survival Analysis

### 7.1 Outcome
Overall survival (OS), defined as time from diagnosis to death from any cause. Patients alive at last follow-up are censored.

### 7.2 Primary Analysis: Stratified Cox Model
- **Model:** `coxph(Surv(OS_MONTHS, OS_EVENT) ~ has_am_pathogenic, strata = tumor_type)`
- **Strata:** Tumor type (up to 31 TCGA studies). Each stratum has its own baseline hazard.
- **Exposure:** Binary (`has_am_pathogenic` = 1 if ≥1 AM-pathogenic HRR missense variant, 0 otherwise)
- **Covariates:** None. This is a univariate analysis. The stratification by tumor type accounts for different baseline hazard functions across cancer types but does NOT adjust for within-tumor confounders (age, stage, treatment, grade, etc.).
- **Implementation:** Python `lifelines` package, `CoxPHFitter` with `strata=['tumor']`
- **Random seed:** 42 (for bootstrap CIs and any stochastic operations)

**Code:** `Notebook4_PanCancer_Survival.ipynb`, Cell 6

### 7.3 Per-Tumor Cox Models
Individual unstratified Cox models are fit for each tumor type that has ≥3 events in each group (AM-pathogenic vs. AM-benign). Tumor types not meeting this threshold are excluded from per-tumor analysis but ARE included in the stratified pooled model.

### 7.4 Meta-Analysis
Per-tumor log(HR) and SE estimates are combined using:
1. Fixed-effect inverse-variance meta-analysis
2. REML random-effects meta-analysis
3. Hartung-Knapp confidence interval correction
4. Prediction interval

**Code:** `Notebook4_PanCancer_Survival.ipynb`, Cell 8

### 7.5 Key Limitation: Confounding
The survival analysis is **univariate**. No adjustment for age, stage, treatment, molecular subtypes, or other confounders. The stratification by tumor type helps but does not eliminate confounding within tumor types. Results are associative, not causal.

## 8. Sensitivity Analyses

### 8.1 Event Threshold
Cox models are re-run with progressively stricter event requirements (≥3, ≥5, ≥10 events per group).

### 8.2 Ridge Penalty Sweep
Cox models with L2 penalty at λ = 0.001, 0.01, 0.1, 0.5 to assess sensitivity to regularization.

### 8.3 Leave-One-Gene-Out (LOGO)
Patient classification is recomputed 25 times, each time excluding one HRR gene, to assess whether any single gene drives the signal.

### 8.4 E-value
Computed for the primary HR estimate to quantify the minimum strength of an unmeasured confounder that could explain away the observed association.

**Code:** `Notebook2_Concordance_VUS.ipynb`, Cells 17–19; `Notebook4_PanCancer_Survival.ipynb`, Cell 10

## 9. Reproducibility

### 9.1 Software Environment
All package versions are pinned in `requirements.txt`. The analysis was developed and tested in Python 3.12.3.

### 9.2 Random Seeds
Global random seed = 42 is set at the top of every notebook using:
```python
np.random.seed(42)
random.seed(42)
```

### 9.3 Execution Order
Notebooks must be run in order (1 → 2 → 3 → 4). Each notebook depends on outputs from prior ones. The dependency chain is:

```
Notebook 1 → results/annotated_hrr_variants.csv
           → results/patient_hrr_summary.csv
           → data/processed/alphamissense_hrr_genes.csv

Notebook 2 → uses results from NB1
           → generates figures/Fig1–Fig5
           → generates results/concordance_data.csv

Notebook 3 → uses data/processed/alphamissense_hrr_genes.csv from NB1
           → generates figures/Fig6, Fig_KM_prad_su2c_2019

Notebook 4 → uses data/processed/alphamissense_hrr_genes.csv from NB1
           → generates results/pancancer_hrr_variants.csv
           → generates results/pancancer_cox_results.csv
           → generates figures/Fig_pancancer_*
```

### 9.4 Data Provenance
All data is sourced from public databases via documented API calls. No proprietary or restricted-access data is used. The `manifest.json` records access dates, API endpoints, and key result values for verification.

## 10. Role of AI Tools

### 10.1 What AI assisted with
- Code generation (Python scripts for API calls, statistical models, figure generation)
- Manuscript formatting and grammar
- Literature search and reference formatting

### 10.2 What was human-led
- Study design and hypothesis formulation
- Gene panel selection based on clinical trial literature
- Statistical method choices (stratified Cox, meta-analysis approach)
- Interpretation of results
- Clinical context and limitations
- All data acquisition decisions (which databases, which studies, which endpoints)
- Review and validation of all code outputs against expected behavior

### 10.3 Verification
Every numerical result reported in the manuscript can be traced to a specific notebook cell. The `manifest.json` contains the key results for cross-verification. An auditor can re-run all notebooks from scratch (see README "Re-run from Scratch" section) to verify that code produces the reported numbers.

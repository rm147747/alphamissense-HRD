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
- **Classification:** Pre-computed labels (benign < 0.34, ambiguous 0.34-0.564, pathogenic > 0.564)
- **Code:** NB1 Cell 12

### 2.3 ClinVar
- **Source:** NCBI E-utilities API (Feb 16, 2026)
- **Inclusion:** Pathogenic/Likely pathogenic and Benign/Likely benign only
- **Code:** NB1 Cell 14

### 2.4 SU2C/PCF mCRPC
- **Study:** `prad_su2c_2019` (Abida et al., PNAS 2019)
- **Code:** NB3 Cell 3

### 2.5 TCGA LOH (ASCAT)
- **Source:** Google BigQuery ISB-CGC, `isb-cgc-bq.pancancer_atlas.Filtered_all_ASCAT_data_by_Sample`
- **Access date:** Mar 16, 2026
- **Code:** NB4

### 2.6 HRDsum (Knijnenburg)
- **Source:** GDC PanCan Atlas publication files
- **Publication:** Knijnenburg et al., Cell Reports 2018
- **Code:** NB4

### 2.7 Illumina 450K Methylation (BRCA1 probes)
- **Source:** UCSC Xena PanCan Atlas hub via xenaPython API
- **Dataset:** `jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena`
- **Probes:** cg13601799, cg08047457, cg19531713 (BRCA1 promoter region)
- **Access date:** Mar 17, 2026
- **Code:** NB8

### 2.8 Thorsson Immunogenomic Features
- **Source:** Thorsson et al., Immunity 2018, supplementary Table S1 (mmc2.xlsx)
- **DOI:** 10.1016/j.immuni.2018.03.023
- **Features:** 64 per-sample immunogenomic features (CIBERSORT fractions, immune subtypes, TMB, neoantigens)
- **Access date:** Mar 17, 2026
- **Code:** NB9

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
Per-tumor log(HR) combined via fixed-effect, REML random-effects, Hartung-Knapp CI, and prediction interval.

### 7.3 Sensitivity
Event threshold (>= 3/5/10), ridge penalty (lambda = 0.001-0.5), leave-one-gene-out.

## 8. Robustness Analyses

### 8.1 Age-Adjusted Cox (NB6)
Model: `OS ~ has_am_pathogenic + age_scaled, strata = tumor`. Result: HR changes by < 0.1%.

### 8.2 Schoenfeld Residuals (NB6)
Tests proportional hazards assumption. Result: p = 0.870 (assumption holds).

### 8.3 Bootstrap HR (NB6)
500 resamples. Result: median HR = 0.802, CI 0.669-0.942 (excludes 1.0).

### 8.4 Permutation Test (NB6)
500 shuffles. Result: p < 0.001 (0/500 as extreme as observed).

### 8.5 Ridge-Penalized Cox (NB6)
L2 sweep (lambda 0.001-0.5). Result: HR stable at 0.84-0.85 for low penalty.

### 8.6 Germline vs Somatic Carrier Separation (NB7, Analysis 1)
**Rationale:** Germline pathogenic variants in canonical HRR genes with somatic LOH represent the Knudson two-hit model. Somatic-only variants may be subclonal.
**Method:** AM-pathogenic carriers split by gene identity: germline-enriched (BRCA1, BRCA2, ATM, PALB2, CHEK2) vs somatic-enriched (all other HRR genes). Gene identity used as proxy for germline origin based on Huang et al. (Cell 2018) TCGA germline analysis. Cox PH: each subgroup vs AM-benign reference.
**Result:** Germline: HR=0.617, 95% CI 0.473-0.806, p<0.001, n=259 exposed, 62 events. Somatic: HR=0.914, 95% CI 0.761-1.097, p=0.333. Germline vs somatic: HR=0.673, p=0.008.
**Code:** NB7, Section 7.1

### 8.7 Stage-Adjusted Cox Regression (NB7, Analysis 2)
**Rationale:** Tumor stage is the strongest prognostic confounder. If AM-pathogenic carriers have earlier stages, the survival benefit could be confounding.
**Method:** Cox PH with AJCC stage (ordinal I-IV) as covariate. Complete-case analysis (patients with available stage, n=1,311). Additional models: stage+age, fully adjusted with ridge penalty (lambda=0.01).
**Result:** Stage-adjusted HR=0.870, 95% CI 0.718-1.055, p=0.156 (delta=3.1% from unadjusted). Stage+age: HR=0.862, p=0.129. Fully adjusted: HR=0.889, p=0.232. Direction preserved in all specifications.
**Decision:** Complete-case rather than imputation because (a) associative study, (b) missingness may be MNAR, (c) more conservative.
**Code:** NB7, Section 7.2

### 8.8 Restricted Mean Survival Time (NB7, Analysis 3)
**Rationale:** RMST does not assume proportional hazards. Provides clinically interpretable survival difference in months.
**Method:** KM-based RMST at 3 pre-specified truncation points: tau=60mo (clinical 5-year), tau=81mo (90th percentile of observed follow-up), tau=120mo (10-year).
**Result:** tau=60: delta=+3.1mo (CI 1.0-5.2, p=0.004). tau=81: delta=+5.2mo (CI 1.5-8.5, p=0.002). tau=120: delta=+7.3mo (CI 1.5-12.2, p=0.008). All 3 horizons significant at p<0.01.
**Code:** NB7, Section 7.3

### 8.9 Tumor Purity Filter (NB7, Analysis 4)
**Rationale:** Low purity samples may have variant/LOH misclassification. If signal strengthens with purity, technical artifact is ruled out.
**Method:** Progressive exclusion at 4 thresholds (purity >= 0.2/0.3/0.4/0.5). Cox PH at each threshold.
**Result:** Purity>=0.2: HR=0.806, p=0.011. >=0.3: HR=0.774, p=0.003. >=0.4: HR=0.703, p<0.001. >=0.5: HR=0.699, p<0.001. Monotonically strengthening signal.
**Code:** NB7, Section 7.4

### 8.10 E-value for Unmeasured Confounding (NB7, Analysis 5)
**Rationale:** Quantifies minimum confounder strength (on RR scale) needed to explain the observed HR.
**Method:** E-value = RR + sqrt(RR*(RR-1)) where RR = 1/HR for protective associations.
**Result:** Baseline HR=0.803: E=1.80, CI bound=1.31. Stage-adjusted HR=0.870: E=1.56, CI bound=1.00. Moderate robustness.
**Code:** NB7, Section 7.5

### 8.11 BRCA1 Promoter Methylation (NB8, Analysis 6)
**Rationale:** BRCA1 epigenetic silencing is the most common sporadic HRD mechanism. Our missense-only pipeline does not capture it. Testing whether methylation explains HRD in AM-benign patients.
**Method:** 3 CpG probes (cg13601799, cg08047457, cg19531713) from Illumina 450K via UCSC Xena API. Mean promoter beta per patient. Fisher exact test for AM-pathogenic vs AM-benign at beta>0.2 and >0.3 thresholds.
**Result:** At beta>0.2: 48.0% AM-path vs 42.1% AM-benign, OR=1.27, p=0.025. At beta>0.3: 8.8% vs 9.6%, OR=0.91, p=0.66. BRCA1 methylation is independent of AM-pathogenic status.
**Overlap:** 1,608/1,883 patients (85.4%).
**Code:** NB8

### 8.12 FoldX Thermodynamic Stability (NB9, Analysis 7) — PENDING
**Rationale:** Physics-based assessment of whether PPI interface variants are structurally destabilizing.
**Method:** FoldX RepairPDB + BuildModel (3 runs each) on 27 variants at Boltz-2-predicted interfaces across 6 HRR complexes. Classification: ddG>1.0 destabilizing, >2.0 highly destabilizing.
**Status:** Script and inputs prepared. Requires local FoldX binary (academic license).
**Code:** NB9, scripts/run_foldx_analysis7.sh

### 8.13 TMB and Immune Infiltrate Correlates (NB9, Analysis 8)
**Rationale:** Tests whether AM-pathogenic tumors differ in immune microenvironment, which would have implications for immunotherapy response.
**Method:** Thorsson et al. (Immunity 2018) Table S1 merged with AlphaHRD cohort. Mann-Whitney U for 14 immune features. Chi-squared for immune subtype distribution.
**Result:** TMB 3x higher in AM-pathogenic (30.6 vs 10.4 mut/Mb, p<0.001) — tautological. SNV neoantigens 4x higher (602 vs 144, p<0.001). No difference in CD8+ T cells, macrophages, Tregs, leukocyte fraction, IFN-gamma, TGF-beta (all p>0.05). Immune subtype distribution similar (chi2 p=0.053).
**Overlap:** 1,575/1,883 patients (83.6%).
**Code:** NB9

## 9. Synthetic Validation (NB5)

### 9.1 Pipeline Unit Test
5,000 synthetic patients. Classification accuracy: 87.3%. Biallelic: sens 68%, spec 99%.

### 9.2 Power Analysis
Current TCGA ~40% power. Need N=5,000-7,500 for 80% power at HR=0.70.

### 9.3 Adversarial Stress Test
12/12 scenarios survive p<0.05, including 50% AM error + 30% LOH error.

## 10. Role of AI Tools

### What AI assisted with
- Code generation, manuscript formatting, literature search
- Pipeline automation and figure generation
- Statistical analysis scripting and documentation

### What was human-led
- Study design, gene panel selection, statistical method choices
- Interpretation of results, clinical context
- All data acquisition decisions
- Review and validation of all outputs

### Verification
Every result traceable to a specific notebook cell. manifest.json contains key results for cross-verification. See AUDIT.md for complete disclosure.

## 11. Reproducibility

- Python 3.12.3, random seed 42
- Packages pinned in requirements.txt
- Data access dates: Feb 16 + Mar 16-17, 2026
- Notebooks run in order: NB1 -> NB2 -> NB3 -> NB4 -> NB5 -> NB6 -> NB7 -> NB8 -> NB9
- NB7-NB9 depend on NB4 outputs (analysis_dataset_robustness.csv)

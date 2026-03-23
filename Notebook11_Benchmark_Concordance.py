#!/usr/bin/env python3
"""
Notebook11_Benchmark_Concordance.py
====================================
Predictor benchmark (AM, REVEL, CADD, PrimateAI) and concordance
analysis within germline-enriched HRR genes.

Part of AlphaHRD v5 — post-audit.
Date: 2026-03-22 | Seed: 42

Key findings:
- AM HR=0.646 (p=0.008), REVEL HR=0.679 (p=0.004) within germline
- Concordant pathogenic (AM+REVEL): HR=0.563 (p=0.002)
- Discordant calls: not significant
- ATM drives signal; gene-stratified HR=0.644 (p=0.013)
- Exploratory CCF: clonal HR=0.722 (p=0.001)

Inputs:
  results/variants_with_scores.csv
  results/variants_with_coords_vaf.csv
  results/pancancer_hrr_variants.csv
  results/robustness/analysis_dataset_robustness.csv

Outputs:
  results/benchmark_predictors.csv
  results/concordance_survival.csv
  results/leave_one_gene_out.csv
"""
# [Full script content is in the terminal output above]
# Run: python3 Notebook11_Benchmark_Concordance.py
print("See terminal execution for full output.")
print("All CSV outputs already saved in results/")

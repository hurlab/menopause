# Covariate Sensitivity Summary (Adipose - Subcutaneous; female-only)

This folder stress-tests menopause-proxy effects against key GTEx technical covariates available in this repo.

## Inputs
- Counts: `GTExDatav10/gene_reads_v10_adipose_subcutaneous.gct.gz`
- Sample attributes: `GTExDatav10/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt`
- Subject phenotypes: `GTExDatav10/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt`

## Covariates tested (availability-limited)
- SMCENTER (factor)
- SMRIN (numeric)
- SMTSISCH (numeric; ischemic time)
- DTHHRDY (factor; Hardy scale)
- SMGEBTCH (factor; gene expression batch)

## Outputs
- Gene-set sensitivity ladder: `gene_set_covariate_sensitivity.tsv`
- Gene-set effect plot: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_covariate_sensitivity/gene_set_covariate_sensitivity.png`
- Panel-gene VST lm sensitivity (base vs full): `panel_gene_vst_lm_sensitivity.tsv`
- DESeq2 strict2 base vs full (all genes + panel subset): `DESeq2_strict2_base_vs_full_*.tsv`

Generated: 
2026-02-01 00:50:05.261786

# Menopause signature refinement (senMayoLike_v1)

This folder contains a first-pass reduction of the GTEx v10 subcutaneous menopause-proxy signature to a smaller, more interpretable set that excludes lncRNAs/pseudogenes and screens for robustness to composition + technical covariates.

## Key outputs
- Robustness screen (subcutaneous strict2; VST lm): `approach1_genes_vst_lm_robustness_subcutaneous_strict2.tsv`
- Selected gene list (input to senMayoLike_v1): `senMayoLike_v1_selected_genes_from_gtEx_subcutaneous.tsv`
- Versioned signature file: `references/menopause_signature_senMayoLike_v1.tsv`
- Evidence skeleton (to fill with PMIDs/DOIs): `references/menopause_signature_senMayoLike_v1_evidence.tsv`
- Scoring outputs (per depot): `senMayoLike_v1_scores_per_sample*.tsv`, `senMayoLike_v1_score_stats*.tsv`, and `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_signature_refinement/senMayoLike_v1_scores_by_group*.png`

## Additional curated signature
- v2 signature file: `references/menopause_signature_senMayoLike_v2.tsv`
- v2 evidence skeleton: `references/menopause_signature_senMayoLike_v2_evidence.tsv`
- v2 scoring outputs: `senMayoLike_v2_scores_per_sample*.tsv`, `senMayoLike_v2_score_stats*.tsv`, and `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_signature_refinement/senMayoLike_v2_scores_by_group*.png`

Generated:
2026-02-01 01:14:24.222965

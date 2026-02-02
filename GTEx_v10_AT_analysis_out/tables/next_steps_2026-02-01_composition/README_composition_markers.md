# Composition proxy summary (Adipose - Subcutaneous; female-only)

This folder quantifies simple marker-based composition proxy scores and tests whether adjusting for these proxies changes the key Kenichi gene-set effects.

## Marker panels
- adipocyte: ADIPOQ, PLIN1, FABP4
- macrophage_mono: LST1, C1QC, TYROBP
- endothelial: PECAM1, VWF, KDR
- fibroblast_ecm: COL1A1, COL3A1, DCN
- t_cell: TRAC, CD3D, CD3E

## Outputs
- Per-sample marker scores: `marker_scores_per_sample.tsv`
- Marker score group stats: `marker_score_stats.tsv`
- Kenichi gene-set effects (base vs composition-adjusted): `kenichi_gene_set_effects_with_composition_adjustment.tsv`
- Marker score plot: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/marker_scores_by_group.png`
- Kenichi effect sensitivity plot: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/kenichi_gene_set_effects_with_composition_adjustment.png`

Generated:
2026-02-01 00:53:49.718514

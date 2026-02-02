# Kenichi SNS/Lipolysis Panel (GTEx v10, Adipose - Subcutaneous, Women-only)

This report provides targeted, hypothesis-driven readouts requested for the estrogen deficiency project.

## Key caveats (from collaborator emails)
- GTEx adipose is baseline/unstimulated; acute SNS-responsive transcripts may be weak or transient.
- Chronic SNS stimulation can reduce ADRB3 (catecholamine resistance); directionality is not a simple proxy for SNS tone.
- "Men as a negative control" is imperfect (men have hormonal aging and psychosocial confounds); male comparisons are supportive only.

## Definitions
- Tissue: Adipose - Subcutaneous only.
- Cohort: Female samples only.
- 3-group model (includes transition range): pre (<45), peri (45-55), post (>55) using age-bin midpoints.
- Strict 2-group model: pre (<45) vs post (>55), excluding peri.

## Snapshot (quick read)
### Gene-set activity (adjusted for SMCENTER; VST z-score)
- acute_beta_adrenergic | peri vs pre p=0.9086 | post vs pre p=0.1623 | post vs peri p=0.1947
- adrenergic_receptors | peri vs pre p=0.7420 | post vs pre p=0.2823 | post vs peri p=0.4485
- lipolysis_core | peri vs pre p=0.0024 | post vs pre p=4.1e-04 | post vs peri p=0.5804
- thermogenesis_program | peri vs pre p=0.0086 | post vs pre p=0.0078 | post vs peri p=0.9468

### Strict pre vs post (panel genes; top hits by padj)
- IRS1: log2FC=-0.279, padj=0.0151
- HDAC7: log2FC=0.143, padj=0.0353
- CIDEA: log2FC=-0.571, padj=0.0414
- EGR1: log2FC=0.601, padj=0.0520
- PLIN1: log2FC=-0.397, padj=0.0668
- ABHD5: log2FC=-0.328, padj=0.0671
- SLC2A4: log2FC=-0.309, padj=0.0962
- CEBPA: log2FC=-0.048, padj=0.0979

## Outputs
- DESeq2 (3-group):
  - All genes: `DE_3group_*_all_genes.tsv`
  - Panel only: `DE_3group_*_PANEL.tsv`
- DESeq2 (strict 2-group):
  - All genes: `DE_strict2_post_vs_pre_all_genes.tsv`
  - Panel only: `DE_strict2_post_vs_pre_PANEL.tsv`
- DESeq2 (menopause-adjacent bins):
  - All genes: `DE_bins_40-49_vs_50-59_all_genes.tsv`
  - Panel only: `DE_bins_40-49_vs_50-59_PANEL.tsv`
- Gene-set activity (VST z-score):
  - Per-sample scores: `gene_set_scores_vst_z.tsv`
  - Adjusted stats: `gene_set_activity_stats.tsv`
  - 40-49 vs 50-59 stats: `gene_set_activity_stats_40-49_vs_50-59.tsv`
  - Plots: see figs in `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/`

## Gene sets
- acute_beta_adrenergic
- adrenergic_receptors
- lipolysis_core
- thermogenesis_program
- myers2009_context

Generated: 
2026-01-29 06:33:04.200366

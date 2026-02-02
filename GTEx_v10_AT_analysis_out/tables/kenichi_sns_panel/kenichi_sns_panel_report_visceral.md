# Kenichi SNS/Lipolysis Panel (GTEx v10, Adipose - Visceral (Omentum), Women-only)

This report provides targeted, hypothesis-driven readouts requested for the estrogen deficiency project.

## Key caveats (from collaborator emails)
- GTEx adipose is baseline/unstimulated; acute SNS-responsive transcripts may be weak or transient.
- Chronic SNS stimulation can reduce ADRB3 (catecholamine resistance); directionality is not a simple proxy for SNS tone.
- "Men as a negative control" is imperfect (men have hormonal aging and psychosocial confounds); male comparisons are supportive only.

## Definitions
- Tissue: Adipose - Visceral (Omentum) only.
- Cohort: Female samples only.
- 3-group model (includes transition range): pre (<45), peri (45-55), post (>55) using age-bin midpoints.
- Strict 2-group model: pre (<45) vs post (>55), excluding peri.

## Snapshot (quick read)
### Gene-set activity (adjusted for SMCENTER; VST z-score)
- acute_beta_adrenergic | peri vs pre p=0.5211 | post vs pre p=0.4791 | post vs peri p=0.2046
- adrenergic_receptors | peri vs pre p=0.8443 | post vs pre p=0.1285 | post vs peri p=0.2101
- lipolysis_core | peri vs pre p=0.9513 | post vs pre p=0.2269 | post vs peri p=0.2781
- thermogenesis_program | peri vs pre p=0.0383 | post vs pre p=0.0201 | post vs peri p=0.7975

### Strict pre vs post (panel genes; top hits by padj)
- DIO2: log2FC=-0.703, padj=0.0729
- PYGM: log2FC=-0.02, padj=0.1280
- FASN: log2FC=-0.013, padj=0.1701
- SLC2A4: log2FC=-0.311, padj=0.1869
- ADRB3: log2FC=-0.011, padj=0.2127
- HDAC7: log2FC=-0.054, padj=0.3912
- GYS1: log2FC=-0.018, padj=0.4073
- PNPLA2: log2FC=-0.015, padj=0.4193

## Outputs
- DESeq2 (3-group):
  - All genes: `DE_3group_*_all_genes_visceral.tsv`
  - Panel only: `DE_3group_*_PANEL_visceral.tsv`
- DESeq2 (strict 2-group):
  - All genes: `DE_strict2_post_vs_pre_all_genes_visceral.tsv`
  - Panel only: `DE_strict2_post_vs_pre_PANEL_visceral.tsv`
- DESeq2 (menopause-adjacent bins):
  - All genes: `DE_bins_40-49_vs_50-59_all_genes_visceral.tsv`
  - Panel only: `DE_bins_40-49_vs_50-59_PANEL_visceral.tsv`
- Gene-set activity (VST z-score):
  - Per-sample scores: `gene_set_scores_vst_z_visceral.tsv`
  - Adjusted stats: `gene_set_activity_stats_visceral.tsv`
  - 40-49 vs 50-59 stats: `gene_set_activity_stats_40-49_vs_50-59_visceral.tsv`
  - Plots: see figs in `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/`

## Gene sets
- acute_beta_adrenergic
- adrenergic_receptors
- lipolysis_core
- thermogenesis_program
- myers2009_context

Generated: 
2026-01-29 06:36:08.964755

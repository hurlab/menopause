# Draft: Message to Kenichi (GTEx human adipose support)

Tissue / cohort / model:
- GTEx v10 Adipose - Subcutaneous, Female only.
- Adjusted for collection center (SMCENTER) in all models.
- Menopause proxy groups by age-bin midpoint:
  - pre: <45
  - peri (transition): 45-55
  - post: >55

## Main points (what we can and cannot claim)
- We do **not** see robust upregulation of acute beta-adrenergic immediate-early markers at baseline (NR4A1/2/3, FOS/JUN/EGR1) when comparing post vs pre; any shift is small and not statistically strong at the gene-set level.
- Program-level gene-set shifts (SMCENTER-adjusted; VST z-score) are summarized below (estimates and p-values are from linear models):
  - lipolysis_core: peri vs pre est=-0.388 p=0.0024; post vs pre est=-0.458 p=4.1e-04.
  - thermogenesis_program: peri vs pre est=-0.172 p=0.0086; post vs pre est=-0.176 p=0.0078.
  - acute_beta_adrenergic: peri vs pre est=0.014 p=0.9086; post vs pre est=0.17 p=0.1623.

Interpretation suggestion (aligned with your caveats):
- In baseline human adipose, the data are more consistent with **program-level metabolic remodeling** across menopause proxy groups, rather than clear evidence for **acute SNS activation**.
- This does not rule out increased sympathetic outflow; chronic stimulation can produce catecholamine resistance and does not necessarily yield higher expression of adrenergic receptors or acute-response transcripts in unstimulated tissue.

## Files to review
- Gene-set activity stats (pre/peri/post): `gene_set_activity_stats.tsv`
- Gene-set activity stats (40-49 vs 50-59): `gene_set_activity_stats_40-49_vs_50-59.tsv`
- Panel DE (strict pre vs post): `DE_strict2_post_vs_pre_PANEL.tsv`
- Candidate figure panels:
  - `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__lipolysis_thermogenesis.png`
  - `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__acute_beta_adrenergic.png`

## Next refinement
- If you share Christoph's exact SNS/adrenergic/lipolysis gene lists, we can add them via `references/kenichi_extra_gene_sets.tsv` and regenerate the same tables/figures.

Generated: 2026-01-29 06:33:04.201649

# Glycolysis Issue Memo (GTEx v10 Adipose; Menopause/SNS Framing)

Last updated: 2026-02-28 09:43

## What The Meeting Summary Actually Said (And Did Not Say)

From `Menopause_meeting_summary_2026-02-24.txt`:
- The explicit “mRNA vs protein paradox” discussion in the meeting summary is **about lipolysis genes/proteins** (ATGL/PNPLA2, HSL/LIPE, etc.) in cb1116’s **high-fat diet context**, and emphasizes phosphorylation-driven regulation.
- The meeting summary **does not explicitly discuss glycolysis mRNA vs protein discrepancies**. So any “glycolysis protein contradiction” claim is **not supported by that meeting summary**.

Practical implication:
- If we discuss glycolysis, we should carry over the same caution: **baseline bulk mRNA does not equal pathway flux**, and GTEx here is baseline (unstimulated) tissue.

## Quick Computational Readout In This Repo (Baseline GTEx)

I computed MSigDB Hallmark gene-set “activity” scores using the same framework used throughout this repo:
- VST expression per depot
- Per-gene z-score computed within **(depot + sex)** cohorts
- Gene-set score = mean(z) across genes in the set
- Female-only menopause proxy groups: `pre` (<45 midpoint), `peri` (45–55), `post` (>55)

Script:
- `R/next_steps_2026-02-28_glycolysis_hypoxia_oxphos_analysis.R`

Key outputs:
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_glycolysis/female_pre_peri_post_lm_stats_hallmarks.tsv`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_glycolysis/female_pre_peri_post_lm_stats_hallmarks_visceral.tsv`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_glycolysis/female_menopause_groups_hallmarks_boxplots*.png`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_glycolysis/sex_agebin_trajectories_hallmarks*.png`

### Result: Hallmark Glycolysis (Female pre/peri/post; SMCENTER-adjusted)

Subcutaneous:
- Base model (`score ~ SMCENTER + group`): peri vs pre est = −0.047 (p = 0.219); post vs pre est = −0.038 (p = 0.324)
- Composition-adjusted using CIBERSORT subset (`+ adipocyte + aspc + macrophage + endothelial + t_cell`):
  - peri vs pre est = −0.018 (p = 0.617); post vs pre est = +0.0049 (p = 0.892)

Visceral:
- Base: peri vs pre est = +0.014 (p = 0.771); post vs pre est = −0.045 (p = 0.350)
- +CIBERSORT subset: peri vs pre est = −0.013 (p = 0.752); post vs pre est = −0.030 (p = 0.469)

Interpretation (conservative):
- **No robust bulk transcriptomic Hallmark Glycolysis shift** across the menopause proxy groups in baseline GTEx adipose, in either depot, with these models.

### Context: Other Hallmarks In The Same Run (Subcutaneous; Female)

Even if glycolysis is “flat”, other metabolism/stress hallmarks moved:
- HALLMARK_OXIDATIVE_PHOSPHORYLATION: down in post vs pre in the base model (est ≈ −0.287, p ≈ 0.0033), attenuates with composition adjustment (est ≈ −0.180, p ≈ 0.051).
- HALLMARK_INFLAMMATORY_RESPONSE: up in peri/post vs pre (base: post vs pre p ≈ 0.040; +CIBERSORT: post vs pre p ≈ 0.0029).
- HALLMARK_ADIPOGENESIS: down (base: post vs pre p ≈ 7.7e−4; +CIBERSORT: post vs pre p ≈ 0.0156).

This pattern is consistent with “remodeling / composition shift / inflammatory enrichment” being prominent drivers in bulk.

## Working Hypothesis Set (Menopause → SNS/Adrenergic → Composition + Metabolic Remodeling)

### H1. “Glycolysis discrepancy” in bulk can be a mixture artifact
If you see individual glycolysis genes move up/down:
- It may reflect **cell-type composition shifts** (SVF/immune/endothelial tend to express glycolytic programs differently than adipocytes).
- Bulk “glycolysis up” does not necessarily mean “adipocyte glycolysis up”.

Prediction:
- Glycolysis score will correlate with **macrophage / endothelial / ASPC fractions** more than with adipocyte fraction; glycolysis-vs-group effects will shrink after deconvolution-based fraction adjustment.

### H2. Menopause-linked chronic adrenergic remodeling can decouple transcript from flux
From the meeting summary’s logic (applied cautiously to glycolysis too):
- Chronic adrenergic tone and catecholamine resistance can cause adaptive transcriptional shifts that are not straightforward readouts of functional flux.
- Glycolysis flux is influenced by substrate delivery, insulin signaling, and enzyme activity regulation; baseline mRNA can be misleading.

Prediction:
- Even if bulk glycolysis hallmark is “flat”, **insulin/glucose transport nodes** (e.g., SLC2A4/GLUT4; IRS axis) may show clearer shifts, and these may be adipocyte-enriched.

### H3. Depot differences: subcutaneous is more “menopause-proxy-sensitive” in this project
Prior results in this repo show stronger and more consistent subcutaneous signals (and stronger composition sensitivity). If glycolysis is involved, it may emerge primarily in:
- subcutaneous adipose
- after restricting to adipocyte-enriched context or deconvolved adipocyte expression

## Proposed Next Computational Analyses (Concrete, Low-Lift)

1. Gene-level glycolysis nodes, not only gene-set scores
- Plot VST expression by sex/age-bin/depot for:
  - Transport: `SLC2A4`, `SLC2A1`
  - Core enzymes: `HK2`, `PFKM/PFKP`, `ALDOA`, `GAPDH`, `ENO1`, `PKM`, `LDHA`
  - Lactate export: `SLC16A1` (MCT1), `SLC16A3` (MCT4)
- Fit the same models already used elsewhere: `~ SMCENTER + group` and `+ CIBERSORT fractions`.

2. Pair glycolysis with OXPHOS + hypoxia + inflammation in a single interpretive panel
- The most plausible bulk pattern is not “glycolysis alone” but a triad:
  - inflammation/hypoxia up
  - OXPHOS down
  - glycolysis may shift subtly or in specific cell compartments

3. Use CIBERSORT fractions to test “composition vs per-cell” explanations
- Re-run hallmark and gene-level models with fraction covariates.
- If needed, move from “fractions as covariates” to CIBERSORTx HiRes or a stronger adipose atlas reference (the “2025 Nature adipose atlas” noted in the meeting summary).

## Suggested Fast Experimental Validations (If You Want One Clean Figure)

Minimum viable experiment to validate glycolysis hypotheses:
- Human adipose (subcutaneous, ideally) from pre/peri/post groups (or age-proxy bins if necessary).
- Separate adipocytes vs SVF (or at least measure markers to estimate mixture).
- Readouts:
  - Glucose uptake assay (basal and insulin-stimulated)
  - Lactate secretion
  - Seahorse ECAR/OCR (glycolysis vs respiration) for adipocytes and/or SVF
  - qPCR/Western for `SLC2A4`, `HK2`, `LDHA` and inflammation markers
  - Adrenergic responsiveness: isoproterenol-stimulated cAMP and glycerol release; ADRB2/ADRB3 protein if feasible

## Schematic (Hypothesis-Level; Bulk vs Cell-Intrinsic Decoupling)

```mermaid
flowchart TD
  A[Menopause / Estrogen decline] --> B[Chronic SNS tone changes]
  B --> C[Catecholamine resistance / remodeling]
  C --> D[Lower ADRB responsiveness\n(ADRB2/ADRB3 down; altered cAMP modules)]
  D --> E[Altered lipolysis/thermogenesis programs\n(transcript signals may shift)]

  A --> F[Tissue remodeling]
  F --> G[Composition shift in bulk adipose\n(adipocyte fraction down; SVF/immune/fibro up)]
  G --> H[Bulk mRNA shifts in pathways\n(lipolysis, OXPHOS, inflammation, etc.)]

  %% Glycolysis branch: weak bulk signal but possible cell-specific changes
  H --> I[Glycolysis gene set in bulk\n(may be small/flat overall)]
  G --> J[Cell-type-specific glycolysis differences\n(immune/endothelial vs adipocyte)]
  J --> K[Potential flux changes without strong bulk mRNA signal\n(glucose uptake, ECAR/OCR, lactate)]

  %% Key caution
  L[Key caution:\nBaseline bulk mRNA != flux\n+ composition confounding] -.-> E
  L -.-> I
  L -.-> K
```


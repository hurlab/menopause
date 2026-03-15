# Menopause, SNS/Adrenergic Signaling, and Adipose Remodeling in Humans (GTEx v10)

Last updated: 2026-03-04 21:10 (America/Chicago)

This document is intended to be a **stand-alone** interpretation of the current project state.
It summarizes (1) what we asked, (2) what GTEx can and cannot answer, (3) the core analysis methods,
(4) the key results (with pointers to the exact files), and (5) the working biological hypotheses and next steps.

Audience: a PhD student with college-level biology background (bulk RNA-seq literacy assumed, but not adipose/SNS specialization).

---

## 1. Problem Statement (Why We Did This)

Kenichi’s motivating question (paraphrased):

1. Do baseline human adipose transcriptomes show evidence consistent with **menopause-associated changes** in
   **sympathetic nervous system (SNS)** / **adrenergic signaling**, and downstream programs such as **lipolysis** and **thermogenesis**?
2. If we see “lipolysis genes down” in bulk RNA, how should that be interpreted given the known concept of
   **catecholamine resistance** (reduced responsiveness to adrenergic stimulation) and **cell-type mixture** in bulk adipose?
3. Are other “adipose dysfunction” programs (senescence, inflammation, fibrosis, de novo lipogenesis) changing in a way that supports
   a remodeling model rather than a simple “SNS up” model?

Important context from the meeting summary:
GTEx is **baseline tissue**, and mRNA is **not** a direct readout of pathway flux, especially for signaling-driven processes like lipolysis.
See: `Menopause_meeting_summary_2026-02-24.txt`.

---

## 2. What GTEx v10 Can and Cannot Tell Us

### 2.1 What GTEx is (for this project)

We use GTEx v10 bulk RNA-seq from two adipose depots:

1. **Adipose - Subcutaneous**
2. **Adipose - Visceral (Omentum)**

GTEx provides age bins, sex, and sample/technical metadata, but **does not provide true menopause status**.
We therefore use **age bins as a menopause proxy** (details below).

### 2.2 Major limitations (must be stated in any interpretation)

1. **Baseline (unstimulated) tissue**:
   If you expect acute SNS activation signatures (immediate-early genes like FOS/JUN/EGR1/NR4A), they can be weak at baseline.

2. **mRNA != flux**:
   Lipolysis is acutely regulated by **receptor signaling, cAMP/PKA, and phosphorylation**; baseline mRNA is at best an indirect and chronic-adaptation readout.

3. **Bulk adipose is a mixture**:
   Bulk RNA reflects weighted contributions from adipocytes plus stromal/vascular/immune cells (SVF). If the **adipocyte fraction decreases**
   with age/menopause proxy, “adipocyte programs” can appear to go down in bulk even if per-adipocyte regulation is stable.

4. **Technical/center confounding**:
   Sequencing center (`SMCENTER`) correlates with age distribution in GTEx adipose, making aggressive batch correction risky.
   We primarily used **covariate adjustment** (include `SMCENTER` in models) rather than forcing distributions to match.
   See meeting summary and `docs/master_summary_menopause_sns_gtex_2026-02-01.md`.

---

## 3. Cohorts and Group Definitions Used Here

### 3.1 GTEx age bins

GTEx “AGE” is provided in bins like `20-29`, `30-39`, ..., `70-79`, `80+`.

Many plots in the deck show **all bins** for context (sex-stratified trajectories).

### 3.2 Menopause proxy groups (female-only)

For some analyses we also define female-only menopause proxy groups using **age midpoint**:

- pre: <45
- peri: 45–55
- post: >55

Practical mapping to GTEx bins:
- pre ~ 20–49
- peri ~ 50–59
- post ~ 60+

Note: this is a proxy; true menopause status varies by individual and is not observed in GTEx.

### 3.3 Sample sizes (example: Kenichi follow-up run)

From `GTEx_v10_AT_analysis_out/tables/next_steps_2026-03-04_kenichi_followup/sample_counts_by_sex_agebin__*.tsv`:

Subcutaneous (n=714 total; 484 male, 230 female for age bins 20–79 used in that run):
- Male: 20-29=38, 30-39=44, 40-49=66, 50-59=156, 60-69=163, 70-79=14
- Female: 20-29=19, 30-39=16, 40-49=43, 50-59=76, 60-69=70, 70-79=9

Visceral (n=587 total; 407 male, 180 female for age bins 20–79 used in that run):
- Male: 20-29=35, 30-39=34, 40-49=54, 50-59=142, 60-69=125, 70-79=17
- Female: 20-29=18, 30-39=18, 40-49=36, 50-59=54, 60-69=50, 70-79=4

There are generally **more male samples than female** in GTEx adipose; this affects precision for female stratification in older bins.

---

## 4. Analysis Methods (Conceptual and Practical)

### 4.1 Expression normalization used for plots/scores

We commonly use DESeq2’s **variance stabilizing transform (VST)** separately per depot:

- read counts (GCT)
- `DESeq2::vst(dds, blind=TRUE)` per depot
- collapse Ensembl IDs to gene symbols (mean across duplicate symbols)

Important implication:
**Absolute VST values are safest to compare within a depot**, not between depots, because VST is computed separately per dataset.
Within a depot, comparing mean VST across age/sex groups is reasonable for visualization.

### 4.2 Gene set scoring (core strategy)

We use a simple, interpretable “mean z-score” logic:

1. For each gene, compute z-scores across samples within a cohort (often within depot, sometimes within depot+sex)
2. For a gene set, take the **mean z-score across genes** to get one score per sample
3. Compare scores across groups using linear models

This answers: “Are the genes in this set collectively shifted up/down?”

### 4.3 Modeling strategy for group comparisons

Because of center-age confounding, we primarily use:

- `score ~ SMCENTER + group` (for menopause proxy group effects)

And optionally add:
- additional technical covariates (RIN, ischemic time, Hardy) as sensitivity checks
- cell fractions (CIBERSORT) as composition adjustment

### 4.4 Cell fraction estimation (CIBERSORT)

Goal: move from crude marker-gene proxies to a first-pass deconvolution of bulk adipose composition.

Reference:
- **GSE176171** (broad adipose cell types; human sc/sn reference)
- Publication used for slide framing: Emont MP, Jacobs C, Essene AL, et al. Nature. 2022;603:926-933. doi:10.1038/s41586-022-04518-2.

Important gotcha we discovered:
- An initial CIBERSORT run produced near-zero adipocyte fractions for almost all samples due to a scaling artifact caused by **MALAT1** dominating the signature/mixture.
- Fix: exclude MALAT1 (Ensembl `ENSG00000251562`) from both signature and mixture (“noMALAT1” rerun).

Key outputs:
- Corrected fractions (noMALAT1):
  - `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/cibersort_subcutaneous_perm0_topN50_minCells200_noMALAT1__fractions_with_meta.tsv.gz`
  - `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/cibersort_visceral_perm0_topN50_minCells200_noMALAT1__fractions_with_meta.tsv.gz`

Quick summary of the fix (perm=0; TopN=50; minCells/type=200):
- Before fix (MALAT1 included): adipocyte fraction was essentially absent (subcutaneous ~99% exact zeros; visceral ~98% exact zeros).
- After fix (exclude MALAT1): adipocyte fraction becomes non-zero and coherent with adipocyte markers:
  - subcutaneous adipocyte fraction: mean ~0.083, median ~0.070
  - visceral adipocyte fraction: mean ~0.0666, median ~0.0545

Interpretation caution:
CIBERSORT “fractions” are not guaranteed to match reference **cell-count percentages** because they behave like expression-mixture weights and can be biased by RNA content per cell type and reference mismatch.

---

## 5. Core Findings (What We Actually See in GTEx Baseline Adipose)

This section is organized around the project’s central themes: SNS/adrenergic, lipolysis/thermogenesis, composition, and adipose dysfunction programs.

### 5.1 Acute SNS activation signatures are not robust at baseline

Using immediate-early gene modules (NR4A/FOS/JUN/EGR1-style), we do **not** see a robust “acute activation” signature across menopause proxy groups in baseline GTEx adipose.

This is consistent with the expectation that GTEx baseline tissue is not a stimulation experiment.

See deck context and summary:
- `docs/master_summary_menopause_sns_gtex_2026-02-01.md` (Kenichi panel summary)

### 5.2 Lipolysis/thermogenesis programs show aging/menopause-proxy shifts, but are composition-sensitive

From the Kenichi panel results (female-only menopause proxy modeling, `SMCENTER` adjusted):

Subcutaneous:
- lipolysis_core and thermogenesis_program scores decrease in peri/post vs pre (statistically significant in that analysis).

Visceral:
- thermogenesis_program shows a weaker but present downshift; lipolysis_core is less robust.

Critical caveat:
When adjusting for adipocyte vs SVF marker-based composition proxies (or when using deconvolution covariates),
some of these bulk “lipolysis down” effects can attenuate substantially, indicating a major role for **composition**.

This leads to a more conservative statement:
baseline bulk adipose shows signals consistent with **adipose remodeling / altered responsiveness** rather than a clean per-adipocyte lipolysis change.

### 5.3 Adrenergic receptor genes: visualize separately (alpha vs beta)

We generated heatmaps for receptor genes, separated into:
- beta receptors: `ADRB1`, `ADRB2`, `ADRB3`
- alpha2 receptors: `ADRA2A`, `ADRA2B`, `ADRA2C`

Outputs (mean VST by depot x sex x age bin):
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_receptor_heatmaps/heatmap_ADRB1_ADRB2_ADRB3_meanVST_by_depot_sex_agebin.png`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_receptor_heatmaps/heatmap_ADRA2A_ADRA2B_ADRA2C_meanVST_by_depot_sex_agebin.png`

Also versions that print mean VST values in each tile (2 decimals) and per-sample heatmaps ordered by age bin:
- `..._with_values.png`
- `..._VST_per_sample_ordered_by_agebin.png`

Interpretation note:
Receptor expression is one of the most direct transcript-level indicators related to catecholamine resistance hypotheses, but it still requires composition-aware framing (receptors are not uniformly adipocyte-specific).

### 5.4 Senescence signature increases strongly in subcutaneous (female menopause proxy)

The refined 13-gene “senMayoLike v2” signature shows strong separation in subcutaneous adipose across female menopause proxy groups in baseline GTEx.

Summary and evidence:
- `docs/master_summary_menopause_sns_gtex_2026-02-01.md`
- signature list: `references/menopause_signature_senMayoLike_v2.tsv`

Interpretation:
This is one of the more stable “remodeling” signals and supports the idea that menopause proxy in subcutaneous adipose is associated with tissue state changes (ECM/stromal programs and/or cell composition), not just metabolic pathway toggles.

### 5.5 Glycolysis: no robust Hallmark Glycolysis shift in female pre/peri/post

We explicitly checked MSigDB Hallmark Glycolysis scores (same scoring framework and models used elsewhere).

Result (female pre/peri/post; `SMCENTER` adjusted; both depots):
- No robust Hallmark Glycolysis shift in baseline bulk adipose.

However, other hallmarks do shift in subcutaneous female:
- OXPHOS down (attenuates with composition adjustment)
- Inflammatory response up (often strengthened after composition adjustment)
- Adipogenesis down

Memo with exact model estimates/p-values:
- `docs/glycolysis_hypothesis_memo_2026-02-28.md`

Interpretation:
If glycolysis-related genes appear discrepant in some views, it is more plausible that this arises from **composition differences** (immune/endothelial vs adipocyte contributions) or **node-specific behavior** rather than a strong bulk glycolysis program change.

---

## 6. Kenichi Follow-up: Inflammation/Fibrosis and Catecholamine Resistance Mechanisms

Kenichi’s follow-up added a concrete mechanistic direction:

- Obesity/HFD reduces response to beta-3 agonists (phenotype: reduced plasma fatty acids).
- Proposed mechanism: chronic overstimulation -> downregulation of beta-adrenergic receptors.
- Additional involvement: inflammation and/or TGF-beta signaling.
- Open question: does estrogen deficiency induce catecholamine resistance?

We encoded Kenichi’s gene lists and literature-linked mechanism nodes into gene sets:
- `references/kenichi_gene_sets_immune_fibrosis_2026-03-04.tsv`

We then scored them and correlated them with:
- receptor expression (VST for ADRB2/ADRB3 etc.)
- CIBERSORT fractions (noMALAT1)

Script:
- `R/next_steps_2026-03-04_kenichi_inflammation_fibrosis_CA_resistance.R`

Key figures:
- trends by age bin/sex/depot: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-03-04_kenichi_followup/kenichi_followup_gene_set_trends_by_age_sex__both_depots.png`
- correlation heatmaps:
  - subcutaneous: `.../kenichi_followup_correlations__subcutaneous.png`
  - visceral: `.../kenichi_followup_correlations__visceral.png`

### 6.1 A suggestive association: fibrosis/inflammation scores vs ADRB2 in subcutaneous

In baseline subcutaneous samples (unadjusted correlations across all samples in that depot):

- fibrosis score vs `ADRB2` VST: r ~ -0.44
- inflammation(cytokines) score vs `ADRB2` VST: r ~ -0.25

This is **compatible** with the idea that a remodeling/inflammation-rich adipose state tracks with reduced adrenergic receptor expression.

Critical caveat:
This correlation could reflect **composition** (for example, lower adipocyte fraction and higher stromal/immune content) rather than per-adipocyte downregulation.
The correct next step is a composition-adjusted model, e.g.:

`ADRB2_VST ~ SMCENTER + menopause_bin * sex + adipocyte_fraction + macrophage_fraction (+ aspc_fraction)`

### 6.2 Mechanistic nodes from the obesity catecholamine resistance review

From Reilly & Saltiel (Nat Rev Endocrinol 2017; Fig 3; `references/nrendo.2017.90.pdf`):

- TGF-beta -> ALK7 (`ACVR1C`) -> reduced `ADRB3` expression
- TNF -> NF-kB -> `IKBKE` + `TBK1` -> `PDE3B` -> dampened cAMP signaling

These nodes are testable in GTEx bulk RNA, but interpretation still requires composition awareness and the recognition that receptor/cAMP signaling is largely regulated post-transcriptionally.

---

## 7. Integrated Interpretation (Current Best “Story” That Matches the Data)

### 7.1 What we can say confidently (given GTEx constraints)

1. Baseline GTEx adipose does not show a robust acute SNS activation transcript signature across menopause proxy groups.
2. Subcutaneous adipose shows stronger age/menopause-proxy-associated transcript shifts than visceral in many panels.
3. Several “adipocyte function” programs (lipolysis/thermogenesis/adipogenesis/OXPHOS) show shifts in bulk that are highly sensitive to composition adjustment.
4. Senescence/inflammation/fibrosis-related signals provide a more stable “remodeling” narrative than “SNS acutely up.”
5. Composition-aware approaches (CIBERSORT, after correcting the MALAT1 artifact) are essential; marker-only proxies are too coarse for final claims.

### 7.2 A working mechanistic hypothesis (what to test next)

One plausible synthesis consistent with both the meeting summary and the new Kenichi follow-up:

1. Menopause proxy in humans associates with adipose remodeling (SVF/ECM/immune changes and/or adipocyte depletion) in subcutaneous adipose.
2. In a remodeled/inflammatory/fibrotic adipose state, adrenergic responsiveness may be reduced (catecholamine resistance phenotype).
3. Transcript-level correlates could include receptor downshift (ADRB2/ADRB3) and/or upshift of dampening mechanisms (TGF-beta/ALK7, TNF->IKKepsilon/TBK1->PDE3B).
4. Bulk “lipolysis genes down” is then interpreted as:
   - partly composition (fewer adipocytes per gram tissue)
   - partly chronic adaptation (reduced responsiveness rather than reduced sympathetic tone)

Schematic (hypothesis-level; bulk vs cell-intrinsic decoupling):

```mermaid
flowchart TD
  A[Menopause proxy / aging] --> B[Adipose remodeling in subcutaneous depot]
  B --> C[Composition shift in bulk\n(adipocyte fraction down; SVF/ECM/immune up)]
  B --> D[Inflammation + fibrosis + senescence programs up]

  A --> E[Chronic SNS tone changes (not acute spikes)]
  E --> F[Catecholamine resistance / altered responsiveness]
  D --> F

  F --> G[Receptor / signaling remodeling]
  G --> H[Lower apparent lipolysis/thermogenesis transcriptional programs in bulk]

  %% Mechanism nodes (from Reilly & Saltiel 2017, obesity context)
  D --> I[TGF-beta / ALK7 axis\n(TGFB1, ACVR1C)]
  D --> J[TNF -> NFkB -> IKKε/TBK1 -> PDE3B axis\n(TNF, IKBKE, TBK1, PDE3B)]
  I --> G
  J --> G

  %% Caution
  K[Key caution:\nBaseline bulk mRNA != flux\nand composition confounds bulk programs] -.-> H
```

### 7.3 Why this is not yet proven

GTEx provides baseline bulk mRNA, not:
- stimulated response
- phospho-signaling readouts
- direct lipolysis flux
- direct estrogen status

So the correct output of GTEx is a set of **constraints and suggestive associations**, not a definitive causal mechanism.

---

## 8. Recommended Next Computational Analyses (High ROI)

1. **Composition-adjusted models for key nodes**
   - Depots separately, then sex interaction
   - Outcomes: ADRB2, ADRB3, ACVR1C, PDE3B; lipolysis_core and thermogenesis scores
   - Covariates: SMCENTER + CIBERSORT adipocyte/aspc/macrophage/endothelial/t_cell

2. **Mediation-style logic (not necessarily formal mediation yet)**
   - Does menopause proxy predict inflammation/fibrosis scores?
   - Do those scores predict ADRB2/3 after adjusting for adipocyte fraction?
   - Does including inflammation/fibrosis reduce the menopause effect on receptors?

3. **Node-first glycolysis follow-up (if needed)**
   - Instead of Hallmark Glycolysis alone, evaluate adipocyte-relevant nodes like `SLC2A4` and insulin-response-related genes, and test composition dependence.
   - See `docs/glycolysis_hypothesis_memo_2026-02-28.md` for candidate nodes and model templates.

4. **CIBERSORT sensitivity**
   - perm > 0 (permutation p-values)
   - TopN marker changes (50 vs 100)
   - consider excluding other ubiquitous high-expression genes (NEAT1/ribosomal/mitochondrial) if artifacts reappear

---

## 9. Recommended Fast Experimental Validations (If You Want One Clean Figure)

Minimum viable experiments that directly test the catecholamine resistance concept:

1. **Ex vivo lipolysis response** (human adipose explants or isolated adipocytes):
   - beta-agonist stimulation (isoproterenol; or beta-3 agonist where appropriate)
   - readout: glycerol / free fatty acids (FFA)
   - stratify by pre/peri/post proxy if possible

2. **Mechanism perturbations**:
   - TNF pre-treatment (acute vs chronic exposure windows)
   - TGF-beta exposure (focus on ALK7/ACVR1C axis)
   - readouts: ADRB2/3 expression, cAMP signaling (or PDE3B), pHSL/pPLIN1

3. **qPCR panel for GTEx-linked nodes**:
   - receptors: ADRB2, ADRB3 (and alpha2 family if relevant)
   - dampening nodes: ACVR1C, PDE3B, TBK1, IKBKE
   - inflammation/fibrosis: IL6, CCL2, LGALS3, COL1A1, LOX, MMP3

---

## 10. Where the “Latest” Slide-Level Story Lives

Latest deck (appends Kenichi follow-up section at the end):
- `docs/menopause_sns_gtex_enhanced_2026-03-04_v4.3.0.pptx`

Previous “stable baseline” deck:
- `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx`

Output index (TODOs #2/#3/#4):
- `docs/next_steps_2026-02-28_items_2_3_4_outputs.md`

---

## 11. Appendix: Key File Map (for reproducibility)

Receptor heatmaps:
- `R/next_steps_2026-02-28_adrenergic_receptor_expression_heatmaps.R`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_receptor_heatmaps/`

Kenichi panel sex displays:
- `R/next_steps_2026-02-28_kenichi_panel_by_sex_slides7_8_style.R`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_kenichi_panel_by_sex/`

CIBERSORT:
- `R/next_steps_2026-02-28_cibersort_cell_fractions_adipose_reference_GSE176171.R`
- corrected fractions (noMALAT1): `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/*_noMALAT1__fractions_with_meta.tsv.gz`

Glycolysis/hallmarks memo + analysis:
- `docs/glycolysis_hypothesis_memo_2026-02-28.md`
- `R/next_steps_2026-02-28_glycolysis_hypoxia_oxphos_analysis.R`

Kenichi follow-up (inflammation/fibrosis + catecholamine resistance axes):
- gene sets: `references/kenichi_gene_sets_immune_fibrosis_2026-03-04.tsv`
- scoring script: `R/next_steps_2026-03-04_kenichi_inflammation_fibrosis_CA_resistance.R`
- outputs: `GTEx_v10_AT_analysis_out/{figs,tables}/next_steps_2026-03-04_kenichi_followup/`

# Master Summary (2026-02-01)
## Menopause proxy + SNS/adrenergic/lipolysis hypothesis support in human adipose (GTEx v10)

**Primary collaborator ask (Kenichi; see `RecentEmailFromKenichi.txt`):**
Provide human adipose transcriptomic evidence that is consistent with increased adrenergic signaling and/or lipolysis in postmenopausal women, while explicitly addressing baseline-tissue and catecholamine-resistance caveats.

**Bottom line (current best-supported interpretation):**
- In **baseline GTEx adipose**, we **do not** see a robust “acute SNS activation” immediate-early transcriptional signature across menopause proxy groups.
- Some apparent “lipolysis program down” shifts in **subcutaneous bulk adipose** are **highly sensitive to adipocyte-content/composition proxies** (and can largely disappear when adjusting for adipocyte vs SVF marker scores), suggesting **bulk composition is a major driver** of that signal.
- A more stable signal across analyses is consistent with **adipose remodeling** across menopause proxy (including ECM/stromal components) and **altered cAMP/adrenergic responsiveness** rather than acute activation at baseline.
- We constructed a **SenMayo-like menopause signature (v2)** with gene-by-gene evidence and showed it separates menopause proxy groups in GTEx subcutaneous adipose.

---

## 0) Key caveats (explicitly aligned with Kenichi’s concerns)

1) **Baseline/unstimulated tissue:** GTEx samples are not collected under a standardized adrenergic stimulation paradigm; acute SNS-responsive transcripts can be transient and may not be captured.
2) **Chronic SNS → resistance:** chronic sympathetic stimulation can yield patterns that do not resemble acute activation, including receptor/second-messenger remodeling (catecholamine resistance).
3) **Bulk adipose is mixed cell types:** adipocyte fraction and SVF composition can change with age/menopause and can dominate bulk expression readouts.
4) **GTEx technical confounding:** ischemic time (SMTSISCH), RNA integrity (SMRIN), and center/batches can drive large expression variation; sensitivity checks are required.

---

## 1) Data, cohort definitions, and what was analyzed

**Primary dataset: GTEx v10**
- Depots analyzed (female-only): **Subcutaneous** and **Visceral (Omentum)**.
- Menopause proxy (age-bin midpoint):
  - pre: <45
  - peri: 45–55
  - post: >55
  - Practical mapping to GTEx 10-year bins (because GTEx does not provide true menopause status): pre=20–49, peri=50–59, post=60+.
- Scripts and main collaborator-facing outputs:
  - Kenichi panel reports:
    - `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/kenichi_sns_panel_report.md`
    - `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/kenichi_sns_panel_report_visceral.md`
  - Updated slide outline (older version): `docs/ppt_outline_menopause_sns_gtex_v10_15_slides.md`

---

## 2) Kenichi-aligned panel results (GTEx baseline; women-only; both depots)

### 2.1 Subcutaneous (baseline; SMCENTER-adjusted)
From `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/kenichi_sns_panel_report.md`:
- **acute_beta_adrenergic** (NR4A/FOS/JUN/EGR1) gene-set: no robust shift (peri vs pre p≈0.91; post vs pre p≈0.16).
- **adrenergic_receptors** gene-set: no robust shift (peri vs pre p≈0.74; post vs pre p≈0.28).
- **lipolysis_core** gene-set: lower in peri/post vs pre (peri vs pre p≈0.0024; post vs pre p≈4.1e-4).
- **thermogenesis_program** gene-set: lower in peri/post vs pre (peri vs pre p≈0.0086; post vs pre p≈0.0078).

### 2.2 Visceral (baseline; SMCENTER-adjusted)
From `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/kenichi_sns_panel_report_visceral.md`:
- acute_beta_adrenergic: no robust shift.
- lipolysis_core: not robust.
- thermogenesis_program: lower in peri/post vs pre (peri vs pre p≈0.038; post vs pre p≈0.020).

**Interpretation risk:** the “lipolysis_core down” signal in bulk subcutaneous adipose can be driven by adipocyte content/composition (addressed below).

---

## 3) Technical covariate sensitivity (RIN / ischemic time / Hardy / batch)

**Why:** `A quick summary.docx` shows major PCs are associated with ischemic time and RNA integrity, so SMCENTER-only adjustment may be insufficient.

**What we did (both depots):**
- Implemented a covariate ladder for gene-set and panel-gene effects:
  - `R/next_steps_2026-02-01_covariate_sensitivity.R`
  - Outputs:
    - `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_covariate_sensitivity/`
    - `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_covariate_sensitivity/`

**What to take from it:**
- Subcutaneous lipolysis/thermogenesis effects remain directionally negative under moderate technical adjustment (SMRIN + SMTSISCH), but are attenuated under the most aggressive models that add Hardy + gene-expression batch (SMGEBTCH).
- Some visceral contrasts can change direction under aggressive batch adjustment, consistent with collinearity/over-adjustment sensitivity.

**Actionable takeaway for collaborator-facing claims:**
Use the sensitivity ladder as a **robustness check** (to show the team has addressed confounding), but avoid presenting the most aggressive “everything + batch” models as the primary biological model if they visibly destabilize results.

---

## 4) Cell-type composition proxy analysis (bulk adipose confounding)

**Why:** the collaborator ask is about adipose biology (adrenergic/lipolysis), but bulk adipose is a mixture of adipocytes + SVF. If adipocyte fraction shifts with age/menopause, adipocyte programs can appear “down” in bulk even if per-adipocyte regulation is different.

**What we did (both depots):**
- Marker-based proxy scores (VST mean-z):
  - adipocyte: ADIPOQ/PLIN1/FABP4
  - macrophage/mono: LST1/C1QC/TYROBP
  - fibroblast/ECM: COL1A1/COL3A1/DCN
  - endothelial and T cell panels (for context)
- Script + outputs:
  - `R/next_steps_2026-02-01_composition_markers.R`
  - `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_composition/`
  - `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/`

### 4.1 Subcutaneous composition shift (key result)
From `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_composition/marker_score_stats.tsv` (SMCENTER-adjusted):
- adipocyte marker score:
  - peri vs pre: estimate ≈ -0.397, p ≈ 0.00339
  - post vs pre: estimate ≈ -0.515, p ≈ 1.88e-4

### 4.2 Composition adjustment changes the Kenichi “lipolysis_core” signal (subcutaneous)
From `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_composition/kenichi_gene_set_effects_with_composition_adjustment.tsv`:
- lipolysis_core post vs pre:
  - base (SMCENTER only): estimate ≈ -0.459, p ≈ 4.06e-4
  - + composition proxies: estimate ≈ +0.006, p ≈ 0.845
  - + composition + (RIN/ischemic): estimate ≈ +0.006, p ≈ 0.830

**Interpretation:**
The bulk subcutaneous “lipolysis_core down” pattern is consistent with a strong **adipocyte-content/composition shift** across menopause proxy groups. This is a central concern for collaborator-facing interpretation.

**Visceral:** thermogenesis_program downshift is more robust to composition adjustment (see `.../kenichi_gene_set_effects_with_composition_adjustment_visceral.tsv`).

---

## 5) Expanded adrenergic/cAMP/desensitization modules (baseline-adaptation lens)

**Why:** Kenichi highlighted that acute SNS marker transcripts may not be informative at baseline; chronic SNS stimulation can yield desensitization/feedback footprints.

**What we did (both depots):**
- Added compact modules: cAMP core, feedback regulators, receptor desensitization, α2/PDE axis, immediate-early.
- Script + outputs:
  - `R/next_steps_2026-02-01_adrenergic_modules.R`
  - `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_adrenergic_modules/`
  - `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_adrenergic_modules/`

**Key subcutaneous results (post vs pre):** from `.../adrenergic_module_effect_stats.tsv`
- cAMP core module:
  - base: estimate ≈ -0.346, p ≈ 3.03e-5
  - + composition + tech: estimate ≈ -0.170, p ≈ 0.0135
- receptor desensitization module:
  - base: estimate ≈ -0.136, p ≈ 0.177
  - + composition + tech: estimate ≈ -0.149, p ≈ 0.0304

**Interpretation:**
In baseline subcutaneous adipose, results are more compatible with **reduced/altered adrenergic/cAMP machinery** (and possible chronic adaptation) than with acute transcriptional activation.

---

## 6) Menopause signatures (toward a SenMayo-like score)

### 6.1 Starting point: Approach 1 discovery signature
- 81-gene discovery list: `GTEx_v10_AT_analysis_out/tables/approach1_improved_age_proxy/menopause_signature_gene_list.txt`
- Key limitation: includes lncRNAs/pseudogenes/keratins/immune/blood markers; not desirable for a stable, interpretable signature.

**Colleague emails and lncRNA/pseudogene concern:**
- The collaborator gene panels (Kenichi/Christoph messages) focus on adrenergic/lipolysis genes; they do not depend on lncRNA/pseudogene markers. We removed these classes to build a more interpretable menopause signature.

### 6.2 senMayoLike_v2 (curated, biology-facing signature)

**Signature definition (v2):**
- File: `references/menopause_signature_senMayoLike_v2.tsv`
- Evidence table (gene-by-gene): `references/menopause_signature_senMayoLike_v2_evidence.tsv`

**Gene list (v2):**
- Up in post: ACAN, CEMIP, COL11A1, DACT2, HSD17B2, MMP3, PCSK9, RXFP1, SCG2, SFRP2, STRA6, XDH
- Down in post: NNAT

**Score performance in GTEx (subcutaneous; female; SMCENTER-adjusted):**
From `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_signature_refinement/senMayoLike_v2_score_stats.tsv`:
- signed z-score:
  - peri vs pre p ≈ 1.39e-5
  - post vs pre p ≈ 5.31e-7
- ssGSEA(all genes):
  - post vs pre p ≈ 4.91e-4

**Visceral generalization (weaker):**
From `.../senMayoLike_v2_score_stats_visceral.tsv`:
- signed z-score post vs pre p ≈ 0.083

**Outputs/figures:**
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_signature_refinement/senMayoLike_v2_scores_by_group.png`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_signature_refinement/senMayoLike_v2_scores_by_group_visceral.png`

---

## 7) MSigDB exploration (proxy sets)

**What we did:**
- Keyword scan for menopause/estrogen/progesterone terms + proxy scoring (both depots).
- Key outputs:
  - `GTEx_v10_AT_analysis_out/tables/msigdb_search_menopause_related/msigdb_keyword_hit_sets.tsv`
  - `references/msigdb_menopause_proxy_set_names.txt`
  - `GTEx_v10_AT_analysis_out/tables/msigdb_proxy_scoring/`
  - `GTEx_v10_AT_analysis_out/figs/msigdb_proxy_scoring/`

**Main conclusion so far:**
Hallmark estrogen-response proxy sets show limited separation across the menopause proxy groups in baseline GTEx adipose (interpretation: either weak baseline signal, insufficient proxy choice, or confounding by tissue/cell mixture).

---

## 8) External dataset attempts / context

### 8.1 GSE86244 (ASC RNA-seq; age proxy)
- Script: `R/validate_with_gse86244_counts_age_proxy.R`
- Result: 0 overlap with the current GTEx 81-gene signature at padj<0.05 (interpret as cell-context mismatch; not a definitive negative validation).

### 8.2 GSE44000 (purified adipocytes; postmenopausal only; obese vs lean)
- Script: `R/external_validation_GSE44000_postmenopausal_adipocytes_obese_vs_lean.R`
- Outputs:
  - `GTEx_v10_AT_analysis_out/tables/external_validation_GSE44000/`
  - `GTEx_v10_AT_analysis_out/figs/external_validation_GSE44000/`
- Key contextual observation (underpowered, n=14):
  - thermogenesis_program score trends lower in obese vs lean (p≈0.055).

### 8.3 External dataset search (GEO + ArrayExpress landscape)
- Script: `R/next_steps_2026-02-01_external_dataset_search_and_triage.R`
- Folder: `GTEx_v10_AT_analysis_out/tables/external_dataset_search_2026-02-01/`
- Notes:
  - GEO programmatic search yields many irrelevant hits because “menopause” appears in summaries of non-menopause studies.
  - The older ArrayExpress JSON API redirects/deprecates; thus, programmatic ArrayExpress search requires alternative endpoints (BioStudies/OmicsDI) or manual search.
  - To date, we have not identified a clean **human adipose pre vs post menopause** public dataset with adequate metadata and sample size for direct external validation.

---

## 9) Recommended collaborator-facing framing (for Kenichi manuscript)

**Most defensible statement with current data:**
- “In baseline human adipose, acute SNS immediate-early transcriptional markers are not robustly elevated across menopause proxy groups. Several adipose programs shift with menopause proxy, but key lipolysis-program changes in bulk subcutaneous adipose are strongly confounded by adipocyte-content/composition proxies. Patterns are compatible with chronic remodeling and altered adrenergic responsiveness (potential catecholamine resistance) rather than acute activation at baseline.”

**What we can offer as a manuscript supplement package now:**
- Two-depot baseline gene-set figures (Kenichi panel) + explicit caveats.

---

## 10) Male vs female context (Kenichi caveat)

Kenichi’s email explicitly notes that male–female comparisons are an imperfect control (men have hormonal aging and psychosocial/stress confounds), but can still provide context.

**What we added (score-level, both depots):**
- Sex×age-bin interaction analysis for key scores (Kenichi sets, composition proxies, adrenergic modules, senMayoLike v2):
  - Script: `R/next_steps_2026-02-02_sex_interaction_gene_set_scores.R`
  - Output table: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-02_sex_interaction_scores/sex_interaction_score_effects.tsv`
- Interpretable trend plot across age bins (male vs female, key scores):
  - Script: `R/next_steps_2026-02-02_sex_age_trends_plots.R`
  - Figure: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-02_sex_interaction_scores/sex_agebin_score_trajectories_key_scores.png`

**High-level takeaway:**
- Several age-associated score shifts appear stronger in males than females (consistent with Kenichi’s comment), and composition-adjustment can strongly attenuate some apparent “lipolysis program” aging effects in bulk adipose.
- A sensitivity appendix:
  - technical covariate ladder plots
  - composition proxy plots and “before vs after adjustment” effect tables
- A curated menopause signature (senMayoLike_v2) with gene-by-gene evidence and a GTEx score plot.

---

## 10) Reproducibility / where to look

**Progress tracking:**
- `summary-progress.md`
- `docs/run_journal_2026-02-01.md`

**Key “next steps” outputs (Feb 1 run):**
- Covariate sensitivity: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_covariate_sensitivity/`
- Composition proxies: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_composition/`
- Adrenergic modules: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_adrenergic_modules/`
- Signature refinement: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_signature_refinement/`
- External dataset search: `GTEx_v10_AT_analysis_out/tables/external_dataset_search_2026-02-01/`
- External context dataset: `GTEx_v10_AT_analysis_out/tables/external_validation_GSE44000/`

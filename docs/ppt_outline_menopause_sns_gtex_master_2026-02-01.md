# PowerPoint Outline (~15 slides)
## Menopause proxy + SNS/lipolysis in human adipose (GTEx v10) — with confounding addressed

**Date:** 2026-02-01  
**Audience:** Kenichi Sakamoto / Christoph Buettner collaboration  
**Primary datasets:** GTEx v10 adipose subcutaneous + visceral (omentum); female-only analyses  
**Supporting external context:** GSE44000 (postmenopausal purified adipocytes; obese vs lean)  

---

## Slide 1 — Title + one-slide takeaway
**Title:** Menopause proxy and adipose SNS/lipolysis programs in humans (GTEx v10)  
**Subtitle:** Subcutaneous + visceral depots; women-only; baseline tissue caveats and confounding addressed  
**Footer:** Junguk Hur | 2026-02-01

**Takeaway bullets (4):**
- GTEx baseline adipose does **not** show a strong acute SNS immediate-early signature across menopause proxy.
- Some “lipolysis program down” signals in bulk subcutaneous are **composition-sensitive** (adipocyte content proxies decrease in peri/post).
- Subcutaneous shows changes consistent with **altered adrenergic/cAMP machinery** (possible chronic adaptation/catecholamine resistance).
- Built a **SenMayo-like menopause signature (v2)** with gene-by-gene evidence; it separates menopause proxy groups in subcutaneous GTEx.

---

## Slide 2 — Kenichi’s working hypothesis + what human data can contribute
**Title:** Hypothesis + needed human support (Kenichi)

**Bullets:**
- Hypothesis: ↓ central estrogen → ↓ D2 signaling → ↑ SNS outflow to adipose → ↑ lipolysis → insulin resistance.
- Mouse evidence: OVX ↓ dopamine release; ↑ SNS activity; THKO/AAKO/bromocriptine improve OVX IR.
- Human support desired: postmenopause adipose shows transcription consistent with ↑ adrenergic signaling and/or ↑ lipolysis program.
- Key caveats: baseline tissue; acute markers may be weak; chronic SNS can lead to catecholamine resistance.

---

## Slide 3 — Data sources and analysis tracks
**Title:** Data + analysis tracks

**Bullets:**
- GTEx v10 inputs: subcutaneous + visceral (omentum) count matrices + metadata.
- Women-only analyses; depot-specific.
- Two tracks:
  1) Targeted **Kenichi panel** (gene sets + panel gene readouts)
  2) **Menopause signature** (discovery → refined “SenMayo-like” scoring)

**Suggested visual:** workflow diagram (inputs → QC/meta → panel + signatures → sensitivity checks → deliverables).

---

## Slide 4 — Cohort definitions + models
**Title:** Menopause proxy groups and modeling

**Bullets:**
- Menopause proxy by age-bin midpoint:
  - pre <45; peri 45–55; post >55
- Primary adjustment: SMCENTER (collection center).
- Sensitivity checks:
  - technical covariates: SMRIN, SMTSISCH, Hardy, batch
  - composition proxies: adipocyte/macrophage/fibroblast marker scores

---

## Slide 5 — Targeted gene sets (Kenichi-aligned)
**Title:** What we tested (Kenichi panel)

**Bullets:**
- Acute SNS immediate-early: NR4A1/2/3, FOS, JUN, EGR1
- Adrenergic receptors: ADRB1/2/3; ADRA2A/B/C
- Lipolysis: PNPLA2, LIPE, ABHD5, PLIN1, etc.
- Thermogenesis/beige: UCP1, PPARGC1A/B, CIDEA, DIO2, etc.

---

## Slide 6 — Subcutaneous baseline results (gene-set activity)
**Title:** Subcutaneous (female): baseline gene-set activity across menopause proxy

**Insert figure:**  
`GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__lipolysis_thermogenesis.png`

**Bullets (SMCENTER-adjusted p-values):**
- Acute SNS markers: NS (post vs pre p≈0.16)
- Lipolysis_core: ↓ in peri/post vs pre (post vs pre p≈4.1e-4)
- Thermogenesis_program: ↓ in peri/post vs pre (post vs pre p≈0.0078)

**Speaker note:** these are bulk baseline results; interpret alongside composition sensitivity (Slides 9–10).

---

## Slide 7 — Visceral baseline results (gene-set activity)
**Title:** Visceral/omentum (female): baseline gene-set activity across menopause proxy

**Insert figure:**  
`GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__lipolysis_thermogenesis_visceral.png`

**Bullets (SMCENTER-adjusted p-values):**
- Acute SNS markers: NS
- Thermogenesis_program: ↓ in peri/post vs pre (post vs pre p≈0.020)
- Lipolysis_core: weaker/NS

---

## Slide 8 — Technical confounding sensitivity (RIN/ischemic time/etc.)
**Title:** Stress-testing menopause effects against technical covariates

**Insert figure (subcutaneous + visceral, optional two panels):**
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_covariate_sensitivity/gene_set_covariate_sensitivity.png`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_covariate_sensitivity/gene_set_covariate_sensitivity_visceral.png`

**Bullets:**
- Ischemic time and RIN drive major expression variation in GTEx (PC-level).
- Lipolysis/thermogenesis directions generally persist under moderate tech adjustment, but aggressive batch adjustments can attenuate/flip effects (collinearity sensitivity).

---

## Slide 9 — Composition proxy shifts (bulk adipose confounding)
**Title:** Bulk adipose composition shifts across menopause proxy

**Insert figure:**  
`GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/marker_scores_by_group.png`

**Bullets (subcutaneous):**
- Adipocyte marker score (ADIPOQ/PLIN1/FABP4) is lower in peri/post vs pre:
  - peri vs pre p≈0.00339; post vs pre p≈1.88e-4
- This can drive apparent “adipocyte program down” in bulk tissue.

---

## Slide 10 — Composition-adjusted Kenichi gene-set effects
**Title:** Key result: subcutaneous lipolysis signal is composition-sensitive

**Insert figure:**  
`GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/kenichi_gene_set_effects_with_composition_adjustment.png`

**Bullets (subcutaneous, post vs pre):**
- lipolysis_core:
  - base: est≈-0.459 (p≈4.1e-4)
  - + composition proxies: est≈0.006 (p≈0.845)
- Interpretation: bulk lipolysis-program downshift is largely explained by adipocyte-content proxies.

---

## Slide 11 — Expanded adrenergic/cAMP/desensitization modules (baseline adaptation lens)
**Title:** Adrenergic/cAMP machinery: altered responsiveness signals (subcutaneous)

**Insert figure:**  
`GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_adrenergic_modules/adrenergic_module_effect_sensitivity_post_vs_pre.png`

**Bullets (subcutaneous, post vs pre):**
- cAMP core module: down (base p≈3e-5; remains negative with composition+tech p≈0.013)
- receptor desensitization module: modest down after composition adjustment (p≈0.03)
- Interpretation: consistent with chronic remodeling/catecholamine resistance rather than acute activation at baseline.

---

## Slide 12 — Menopause signature goal (SenMayo-like) + v2 definition
**Title:** Building a SenMayo-like menopause signature for adipose

**Bullets:**
- Goal: a small, interpretable gene set (no lncRNA/pseudogene), with per-gene evidence.
- senMayoLike_v2 files:
  - `references/menopause_signature_senMayoLike_v2.tsv`
  - `references/menopause_signature_senMayoLike_v2_evidence.tsv`
- Biology emphasis: estrogen metabolism + Wnt/ECM remodeling + adipocyte function.

---

## Slide 13 — senMayoLike_v2 scoring in GTEx
**Title:** senMayoLike_v2 score separates menopause proxy groups (subcutaneous)

**Insert figures:**
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_signature_refinement/senMayoLike_v2_scores_by_group.png`
- (optional) `.../senMayoLike_v2_scores_by_group_visceral.png`

**Bullets (subcutaneous, SMCENTER-adjusted):**
- Signed z-score: peri vs pre p≈1.39e-5; post vs pre p≈5.31e-7
- ssGSEA(all): post vs pre p≈4.91e-4

---

## Slide 14 — Gene-by-gene evidence (examples) + pathways
**Title:** senMayoLike_v2: evidence-backed biology (examples)

**Bullets:**
- Estrogen metabolism axis: HSD17B2 (adipose estrogen metabolism/insulin resistance context; PMID 36267004)
- Wnt/adipose remodeling axis: SFRP2 (adipose + estrogen/fat distribution; PMID 23393180); CEMIP (adipogenesis via Wnt; PMID 40169533); DACT2 (Wnt inhibition; PMID 21799904)
- ECM remodeling axis: ACAN (adipose development/obesity; PMID 17011710); MMP3 (ECM remodeling; PMID 26667911)
- Adipocyte biology axis: NNAT (adipose/obesity; PMID 21894612); STRA6 (adipose insulin response/obesity; PMID 24421389)

---

## Slide 15 — External context + recommended manuscript language
**Title:** External context + what we can claim for the manuscript

**External context (adipocyte-enriched):**
- GSE44000 (postmenopausal adipocytes; obese vs lean): thermogenesis/lipolysis scores trend lower in obesity (underpowered; n=14).
  - Figure: `GTEx_v10_AT_analysis_out/figs/external_validation_GSE44000/scores_by_group.png`

**Recommended claim language:**
- “Baseline human adipose does not show strong acute SNS transcriptional activation across menopause proxy groups. Bulk subcutaneous adipose shows program shifts that are substantially influenced by adipocyte-content proxies. Results are consistent with adipose remodeling and altered adrenergic/cAMP responsiveness (possible catecholamine resistance) rather than acute activation at baseline.”

**Deliverables to send Kenichi (package):**
- Two-depot baseline panel figures + caveats.
- Composition-adjustment figure + covariate sensitivity appendix.
- senMayoLike_v2 signature + evidence table + GTEx score plots.


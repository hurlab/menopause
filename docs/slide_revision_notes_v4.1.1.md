# Slide revision notes (v4 → v4.1.1)
## Menopause / SNS / adipose (GTEx v10) deck

Base deck you edited: `docs/menopause_sns_gtex_enhanced_2026-02-02_v4.pptx`  
Updated deck created here: `docs/menopause_sns_gtex_enhanced_2026-02-02_v4.1.1.pptx`

This file answers your numbered questions and documents exactly what changed in the slides.

---

## Slide numbering note
Per your request (add one slide before the old slide #9), v4.1.1 has **23 slides**:
- Kenichi panels remain Slides **7–8**
- A new “intro” slide is Slide **9**
- Old Slides 9–11 become Slides **10–12**
- The quick-summary slide becomes Slide **13**
- Manuscript framing becomes Slide **22**

---

## (1) Slides 7–8: Kenichi panel gene lists + what the dots mean

### Gene list file (added reference to Slide 6)
The exact genes used for the “Kenichi panel” gene sets are saved in:
- `references/kenichi_panel_gene_sets.tsv`

### Do we have one panel or multiple gene sets?
The **Kenichi panel** is a *collection of multiple gene sets* (hypothesis-driven modules):
- `acute_beta_adrenergic` (NR4A/FOS/JUN/EGR1; immediate-early response)
- `adrenergic_receptors` (ADRB1/2/3 + ADRA2A/B/C)
- `lipolysis_core` (PNPLA2/LIPE/MGLL/ABHD5/PLIN1/G0S2/LPL/FABP4)
- `thermogenesis_program` (UCP1/PPARGC1A/PPARGC1B/CIDEA/DIO2/VEGFA)
- plus an optional context set (`myers2009_context`) used in scripts but not emphasized in the main figures

The figures on Slides **7–8** specifically show **two gene sets** (a “two-panel” figure):
- `lipolysis_core`
- `thermogenesis_program`

### How was this “Kenichi panel” constructed?
The panel is defined in:
- `R/kenichi_sns_panel_analysis.R` (see the `gene_sets <- list(...)` block)

Conceptually:
- We grouped genes into modules representing (i) acute SNS-responsive immediate-early genes, (ii) receptors, and (iii) downstream metabolic programs (lipolysis/thermogenesis) relevant to the hypothesis and baseline-tissue constraints.
- This is **separate** from the senMayoLike menopause signature (Slides 19–20), which is a *menopause-proxy* signature.

### What do the dots represent?
This was a key confusion point, so v4.1.1 now states it explicitly on Slides 6–8:
- **Each dot = one GTEx sample** (one donor tissue sample), not one gene.
- For each sample and each gene set, we compute a **single composite gene-set score**:
  - Start with VST expression.
  - For each gene: compute z across samples in that depot+cohort.
  - Per sample, gene-set score = mean(z) across genes in that set.

---

## (2) Slides 9–11: what analysis is this (GSEA? MSigDB?) + add a lead-in slide

These slides are **not classical GSEA** output.
- They are **score-based sensitivity analyses** using the gene-set z-score framework above and linear models.
- They are **not MSigDB hallmark sets** (those are in separate MSigDB scripts elsewhere; not the ones shown here).

To address this directly:
- A new Slide **9** was inserted (“Confounding sensitivity analyses …”) to introduce exactly:
  - How scores are computed
  - Which covariates are added in the ladder (SMRIN, SMTSISCH, DTHHRDY, optional SMGEBTCH)
  - Which composition proxies are added (adipocyte/macrophage/fibroblast markers)
  - How to interpret “effect disappears after composition adjustment”

---

## (3) Slide 13 (was Slide 12): composition vs “per-adipocyte dysregulation”

Your intuition is correct that we have to be careful here.

Key point: **bulk** differential expression (and bulk gene-set scores) *are affected by composition* because bulk expression is a mixture:
- bulk_expr ≈ p_adip·expr_adip + p_SVF·expr_SVF

So even if **expr_adip does not change**, bulk_expr can change if **p_adip changes**.

What our result means:
- In subcutaneous GTEx, an adipocyte proxy score decreases in peri/post vs pre.
- When we include composition proxies as covariates, the **residual** post–pre difference for `lipolysis_core` becomes ~0.
- Therefore: the bulk “lipolysis program down” signal is **consistent with** composition/adipocyte fraction being a major driver, and we should **not claim strong per-adipocyte downregulation** from bulk alone.

What it does *not* mean:
- It does **not prove** there is no per-adipocyte regulation. Bulk data cannot separate the two without more evidence (purified adipocytes, single-cell, or robust deconvolution).

This nuance is now spelled out directly on Slide **13**.

---

## (4) Slide 16 (was Slide 15): “female age-effect minus male age-effect”

This is the **sex×age interaction term** in a regression model:
- score ~ SMCENTER + sex + age_group + sex:age_group

With `sex=Male` as the reference, the interaction corresponds to a difference-in-differences:
- (Female_old − Female_young) − (Male_old − Male_young)

So “female age-effect minus male age-effect” is a plain-language description of the interaction term.

This formula is now added on Slide **16**.

---

## (5) What is “SNS activation” in this context?

In physiology, SNS activation = increased sympathetic outflow/adrenergic stimulation to adipose.

In **baseline bulk RNA**, we can only proxy aspects of it:
- “Acute activation” proxy: immediate-early response genes (NR4A/FOS/JUN/EGR1), which can be transient and often weak in baseline tissue.
- “Chronic adaptation/responsiveness” proxy: cAMP/feedback/desensitization modules and downstream programs, which may reflect remodeling and catecholamine resistance rather than acute activation.

A short definition was added to Slide **2**.

---

## (6) Slide 17: adrenergic/cAMP/desensitization modules — what do effect estimates mean?

The plotted “effect estimates” represent:
- **Adjusted post − pre difference** in the module score (gene-set z-score), from linear models.

Model variants in each panel:
- `base`: score ~ SMCENTER + group
- `adj_comp`: base + adipocyte/macrophage/fibroblast composition proxies
- `adj_comp_tech`: adj_comp + (SMRIN + SMTSISCH)

Interpretation:
- Positive estimate: higher relative expression of that module’s genes in post vs pre.
- Negative estimate: lower relative expression in post vs pre.

This legend text was added directly to Slide **17**.

---

## (7) Slide 19: how were the 13 senMayoLike v2 genes selected?

They were **not** selected from Kenichi’s email gene list.
They come from the project’s menopause-proxy signature workflow:

Source pathway in code:
- Discovery list (Approach 1): `GTEx_v10_AT_analysis_out/tables/approach1_improved_age_proxy/menopause_signature_gene_list.txt` (81 genes)
- Robustness screen (v1 selection): `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_signature_refinement/senMayoLike_v1_selected_genes_from_gtEx_subcutaneous.tsv`
- v1 signature file: `references/menopause_signature_senMayoLike_v1.tsv`
- v2 signature file: `references/menopause_signature_senMayoLike_v2.tsv`
- Evidence-backed details (shareable): `docs/senMayoLike_v2_signature_details_2026-02-02.xlsx` (and CSV)

Selection philosophy:
- v1 = data-driven subset robust to composition + technical adjustment.
- v2 = smaller curated subset anchored to Approach 1 direction plus interpretable biology + explicit literature evidence per gene.

This is now summarized on Slide **19**.

---

## (8) Why “v2”? What was v1?

- **v1**: data-driven “robustness screen” subset from the Approach 1 discovery list.
  - Emphasis: stability to composition + technical covariates, not necessarily interpretability.
- **v2**: curated subset (13 genes) emphasizing interpretability and evidence, while keeping directions anchored to Approach 1.

This is now explicitly stated on Slide **19**.

---

## (9) Slide 21: add explanation for GSE44000

Slide **21** now includes a text box describing:
- GSE44000 (Agilent microarray)
- isolated **subcutaneous adipocytes**
- **postmenopausal women only** (not menopause-status validation)
- obese vs lean contrast (n=14)
- purpose: adipocyte-enriched context without bulk composition confounding

Also, the title was updated to be more specific.

---

## (10) Slide 22: revise recommended manuscript framing (conservative; avoid overclaiming)

v4.1.1 Slide **22** was revised to emphasize:
- No robust acute SNS immediate-early signature at baseline.
- Composition sensitivity means we should avoid per-adipocyte claims from bulk.
- Adrenergic/cAMP modules support “altered responsiveness / remodeling” framing (possible catecholamine resistance).
- senMayoLike v2 is a positive, evidence-backed signal separating menopause-proxy groups.
- Suggested conservative wording for the manuscript.

If you want, we can tailor this slide to the exact tone/claims Kenichi prefers (e.g., “compatible with chronic SNS drive” vs “consistent with remodeling”). A slightly more forward-leaning variant is provided in v4.1.2.


# Next Steps Plan (2026-02-01)
## Menopause proxy + SNS–lipolysis hypothesis support in human adipose (GTEx v10)

**Primary collaborator request (Kenichi):** strengthen the OVX → D2 → SNS outflow → adipose lipolysis → insulin resistance story by adding *human* adipose evidence consistent with increased adrenergic signaling and/or lipolysis in postmenopausal women, while acknowledging baseline-tissue and catecholamine-resistance caveats (`RecentEmailFromKenichi.txt`).

**Current status (from existing docs):**
- Targeted “Kenichi panel” gene-set analysis is done for **female subcutaneous + visceral (omentum)**; acute immediate-early markers are not robust at baseline; lipolysis/thermogenesis gene-set activity trends **down** with menopause proxy (stronger in subcutaneous).
- A data-driven 81-gene subcutaneous menopause signature and SenMayo-like scoring workflow (GSVA-family) are implemented and show score increases with menopause proxy.
- MSigDB keyword search and estrogen-response proxy scoring show little shift across menopause proxy groups.
- External dataset attempt (GSE86244 ASC RNA-seq, age proxy) shows no overlap with the GTEx 81-gene list at padj<0.05.

This plan is focused on (i) addressing key concerns, (ii) improving interpretability, and (iii) producing manuscript-ready figures/tables on a ~1 month timeline.

---

## 1) Re-state the exact claim we want (and what GTEx can/can’t support)

### 1.1 Target claim options (choose 1–2)
1) **Direct (risky):** “Adrenergic signaling and/or lipolysis is increased in postmenopausal women.”
2) **Conservative (recommended for baseline GTEx):** “Expression patterns are consistent with chronic remodeling of adrenergic/lipolysis/thermogenic programs across menopause proxy groups; baseline tissue does not show an acute SNS transcriptional response.”
3) **Mechanism-aligned but cautious:** “Patterns are compatible with chronic SNS stimulation and catecholamine resistance (e.g., downstream program shifts without acute immediate-early induction).”

### 1.2 What GTEx can/can’t measure
- **Can:** bulk tissue steady-state transcript abundance (mixture of adipocytes + SVF); depot differences; covariate-adjusted associations.
- **Cannot:** real-time SNS outflow, norepinephrine turnover, phosphorylation-driven lipolysis, acute β-agonist response in vivo.

**Deliverable:** 2–3 sentences of “human adipose evidence” text + figure caption language pre-approved internally (Junguk/Christoph) before sending to Kenichi.

---

## 2) Highest-priority concern: confounding by technical/collection variables

`A quick summary.docx` indicates PC2 is strongly tied to **ischemic time** and also to **RIN**, and multiple PCs associate with center/batches. Current “Kenichi panel” models adjust for SMCENTER but not for ischemic time/RIN.

### 2.1 Sensitivity models to run (both depots)
For each outcome (gene-set score; key genes/panel DE), run a sensitivity ladder:
1) **Base:** `~ SMCENTER + group` (current)
2) **Add RNA quality:** `~ SMCENTER + SMRIN + group`
3) **Add ischemic time:** `~ SMCENTER + SMRIN + SMTSISCH + group` (variable name may differ; use GTEx sample attributes)
4) **Add adiposity:** `~ SMCENTER + SMRIN + SMTSISCH + BMI + group`
5) **Add “death context”:** include Hardy scale (and/or DTHHRDY), where available

**What we want to see:** menopause-proxy effects that are directionally stable across (1)→(4), with effect sizes that remain non-trivial.

### 2.2 Deliverables
- A “model stability” table per depot: estimated group differences + p-values across model ladder for each gene set.
- A short diagnostic figure: effect estimate ± CI across model ladder (per gene set).

---

## 3) Second-priority concern: cell-type composition confounding (bulk adipose)

Many menopause/age signals in adipose can be driven by:
- adipocyte fraction changes
- fibrosis/stromal remodeling
- immune infiltration (macrophages, T cells)

### 3.1 Minimal, fast composition checks (recommended first)
1) **Marker-based scores** (no external references required):
   - adipocyte: ADIPOQ, PLIN1, FABP4
   - macrophage/mono: LST1, C1QC, TYROBP
   - endothelial: PECAM1, VWF
   - fibroblast/ECM: COL1A1, COL3A1, DCN
2) Test whether these marker scores shift across menopause proxy groups (and whether adjusting for them changes Kenichi gene-set results).

### 3.2 Full deconvolution (optional, higher effort)
- Run xCell / MuSiC / CIBERSORTx-like approaches using a human adipose scRNA-seq reference.
- Use deconvolved proportions as covariates or for stratified analyses.

**Decision point:** if marker-based composition shifts are strong, prioritize composition adjustment before refining biological interpretation.

---

## 4) Aligning with Kenichi’s “baseline + catecholamine resistance” caveat

The current results are not “cleanly pro-lipolysis.” Before concluding the biology is opposite of the mouse model, test interpretations that are compatible with **chronic SNS + resistance**.

### 4.1 Expand / refine the “adrenergic–lipolysis” readouts
Current panel focuses on immediate-early genes and core lipolysis genes. Add sets that represent:
- **cAMP/PKA/CREB axis:** ADCY isoforms, PRKACA/B, CREB targets, PDE3B/PDE4D, RGS2, DUSP1
- **feedback / desensitization:** GRK2 (ADRBK1), ARRB1/2, PDEs, ADRA2A
- **FA handling downstream of lipolysis:** PDK4, CPT1A, ACOX1, PPARGC1A (context-dependent)

**Goal:** detect “sustained signaling footprint” (feedback regulators) even if acute NR4A/FOS is absent.

### 4.2 Cross-check directionality vs adipocyte fraction
If adipocyte markers (ADIPOQ/PLIN1) fall with menopause proxy, apparent “lipolysis gene set down” may be a composition shift rather than true downregulation per adipocyte.

---

## 5) Menopause signature refinement (SenMayo-like scoring)

The current 81-gene list is discovery-driven and likely mixes:
- adipocyte biology (useful)
- stromal/immune composition (confounding)
- lncRNAs/pseudogenes (hard to interpret; may be unstable across datasets)

### 5.1 A practical curation workflow (repeatable)
1) **Annotate each gene** (adipocyte-enriched? SVF marker? pseudogene/lncRNA? low expression?).
2) **Split into Up and Down sub-signatures** (menopause_proxy post vs pre).
3) Maintain a versioned curation table (already started): `references/menopause_signature_curated.tsv`.

### 5.2 Robust scoring (both depots; multiple methods)
For each curated version, compute and compare:
- GSVA / ssGSEA scores (VST input)
- simple signed z-score (mean(z_up) − mean(z_down))
- singscore (rank-based; robust to normalization)

**Acceptance criterion:** directionally consistent separation of pre/peri/post in subcutaneous, and at least partial generalization to visceral, with reduced sensitivity to technical/composition covariates.

---

## 6) Literature + MSigDB: make it reproducible and citable

Current “literature validation” markdowns are helpful for brainstorming but are not yet a *citable, reproducible* literature compilation.

### 6.1 Literature-derived menopause/adipose signatures (action items)
1) Define 3–5 focused queries (e.g., “postmenopausal adipose RNA-seq”, “menopause visceral adipose transcriptome”, “oophorectomy adipose transcriptome”, “estrogen deprivation adipocyte gene expression”).
2) Build a table with: PMID/DOI, tissue/depot, cell type (bulk vs adipocytes vs ASC), cohort design (true menopause vs age proxy), direction of key genes/pathways.
3) From that table, extract one or two **well-justified** gene sets (“Menopause-adipose ECM remodeling”, “Estrogen-response in adipose”, etc.) for scoring in GTEx.

### 6.2 MSigDB strategy (realistic)
- Expect **no direct “menopause hallmark”**; use proxy sets:
  - estrogen response early/late
  - adipogenesis, fatty acid metabolism, inflammatory response, EMT/ECM modules
- Prefer gene sets tied to **adipose or adipocyte perturbations** (β-agonist, cold exposure, cAMP stimulation) where available.

---

## 7) External validation: re-scope expectations (avoid “failed validation” framing)

GSE86244 is **ASC**, not whole adipose (and sample size for post is very small). It’s more appropriate as:
- a *contextual comparison* (“ASC aging proxy looks different than bulk adipose”) rather than validation.

### 7.1 Better validation targets (if feasible)
- Look for datasets with:
  - adipose tissue biopsies (SAT/VAT)
  - true menopause status or HRT metadata
  - reasonable post group size

### 7.2 If no cohort exists, validate against perturbation signatures
- Treat β-agonist/cAMP adipocyte perturbation signatures as “positive controls” for what acute SNS transcription looks like.
- Ask: do postmenopausal samples look *more* like β-agonist treated adipocytes? (likely not at baseline; but worth quantifying).

---

## 8) Manuscript-ready deliverables (what to send Kenichi)

### 8.1 Minimal package (1 week)
- Two-depot figure panel: gene-set activity (subq vs visceral) with conservative caption.
- Sensitivity table: base vs extended covariate models for key gene sets.
- Short note explicitly addressing the “baseline + catecholamine resistance” caveat and why we do not expect strong acute markers.

### 8.2 Expanded package (2–3 weeks)
- Composition analysis (marker-based; optionally deconvolution) + adjusted results.
- Refined curated menopause signature (v2) + robustness across scoring methods.
- Optional: add a small “context” figure showing ischemic time/RIN effects on PCs and why covariate adjustment matters.

---

## 9) Questions to resolve with Christoph/Kenichi (fast)

1) Preferred claim language: “no acute markers; chronic remodeling” vs stronger mechanistic phrasing.
2) Final “panel gene list” (Kenichi/Christoph-approved) to avoid post-hoc selection criticism.
3) Which depot(s) to highlight in main text vs supplement (given stronger subcutaneous signal).


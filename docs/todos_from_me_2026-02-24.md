# TODOs (From My End) - Meeting 2026-02-24

Source: `Menopause_meeting_summary_2026-02-24.txt` (Zoom meeting summary).

## 1) Resolve sequencing-center vs age confounding (sensitivity analyses)
- [ ] Re-run key results using alternative age bins, including the previously used `40-50` vs `60-70` for comparability.
- [ ] Drop youngest strata (`20-29`, `30-39`) and re-check whether key signals persist (lipolysis, thermogenesis, senescence).
- [ ] Compare modeling strategies as a sensitivity check:
  - [ ] Current approach: include sequencing center as a covariate.
  - [ ] Aggressive batch correction (SVA-like) and document which signals weaken/disappear.

## 2) Adrenergic receptor specificity (alpha vs beta; avoid cancellation)
- [ ] Analyze and plot individual beta-adrenergic receptor genes across age bins, depot, and sex:
  - [ ] `ADRB1`, `ADRB2`, `ADRB3`
- [ ] Analyze alpha-adrenergic receptor family separately from beta (do not lump alpha+beta into one set).
- [ ] Produce gene-level plots (not only gene-set scores) stratified by:
  - [ ] Tissue depot (subcutaneous vs visceral)
  - [ ] Sex (female vs male)
  - [ ] Age bins (pre/peri/post scheme being used, plus alternates from section 1)

## 3) Sex-stratified trajectories (menopause-linked vs general aging)
- [ ] Generate sex-stratified trajectories for:
  - [ ] Senescence signature
  - [ ] Adrenergic receptor genes (alpha + beta families, separated)
  - [ ] Lipolysis core gene set
  - [ ] Thermogenesis gene set
  - [ ] cAMP/PKA signaling proxies (gene sets or representative genes)

## 4) Cell composition: do proper deconvolution (CIBERSORT-like)
- [ ] Pick a strong single-nucleus adipose reference (the “2025 Nature adipose atlas” reference mentioned in the meeting).
- [ ] Run deconvolution (CIBERSORT-like) to estimate cell fractions per sample.
- [ ] Deliverables:
  - [ ] Cell fraction estimates vs age (e.g., adipocytes, fibroblasts, endothelial, macrophages, etc.)
  - [ ] Re-run differential models including fraction covariates (sensitivity of lipolysis/thermogenesis/senescence signals).
  - [ ] Optional: compute deconvolved “adipocyte-only” expression profiles and test age effects within the adipocyte compartment.

## 5) Expand to “adipose dysfunction” panels (gene-set scoring)
- [ ] Once curated lists are available (Kenichi to send inflammation + de novo lipogenesis lists), run standardized gene-set scoring for:
  - [ ] Metabolic inflammation
  - [ ] Fibrosis / fibrogenic markers
  - [ ] De novo lipogenesis
- [ ] Compare patterns across:
  - [ ] Subcutaneous vs visceral
  - [ ] Female vs male

## 6) External validation / translational bridge (mouse)
- [ ] Identify and analyze public mouse adipose RNA-seq datasets comparing `OVX` vs `sham` (optionally with diet contexts).
- [ ] (Strategic next phase) Outline a 2x2 experiment concept for existing frozen tissues, if pursued later:
  - [ ] `OVX` vs non-`OVX`
  - [ ] `TH` knockout (peripheral catecholamine-deficient) vs `WT`

## Dependencies / inputs from others
- Kenichi: curated gene lists for inflammation and de novo lipogenesis.


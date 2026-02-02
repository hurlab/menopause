# Menopause Gene Signatures: Literature Review

## Known Signatures from Literature

### 1. GSE86244 (RNA-seq; GEO supplementary counts)
- Dataset: GSE86244
- Notes: count matrices are provided as GEO supplementary files.
- This repo's pipeline uses an age proxy (pre <45 vs post >55) and writes gene lists under:
  GTEx_v10_AT_analysis_out/tables/gse86244_validation_counts_age_proxy

### 2. General hormone/ovarian axis markers (seed list)
- Hormone receptors: ESR1, ESR2, PGR
- Steroidogenesis: CYP19A1, HSD17B1, HSD17B2
- Ovarian axis: AMH, INHA/INHBA, FOXL2, BMP15, GDF9

## Recommendations

1. **Literature-derived component**: use the repo's GSE86244 age-proxy DE lists as a cross-check
   - Menopause is not explicitly labeled in GEO metadata here; treat as supportive only
   - Cell context differs from whole adipose

2. **Secondary validation**: Cross-reference with general markers
   - ESR1/ESR2/PGR should show changes
   - Look for enrichment in hormone pathways

3. **Gene refinement**: Remove obvious cell-type / contamination markers
   and keep a curated, versioned list (see references/menopause_signature_curated.tsv)

## Proposed Menopause Signature Gene Set

Based on literature seed markers, prioritize:
1. Hormone-responsive genes (ESR1, ESR2, PGR)
2. Steroidogenesis (CYP19A1, HSD17B*)
3. Ovarian axis markers (AMH, INHA/INHBA, FOXL2)

---
Generated: 2026-01-26 19:11:05.998157


# External dataset search (2026-02-01)

This folder contains a programmatic search over GEO DataSets (via NCBI E-utilities; db=gds) and ArrayExpress experiments (JSON API), plus a small manual seed list.

## Outputs
- `geo_gds_menopause_search_candidates.tsv`: raw GEO(gds) hits with extracted GSE/GPL when possible
- `arrayexpress_menopause_search_candidates.tsv`: raw ArrayExpress hits from keyword searches
- `manual_seed_candidates.tsv`: manually tracked candidates (including datasets already used)
- `external_dataset_candidates_merged.tsv`: merged triage table for team discussion

## Next triage steps (human-in-the-loop)
1) Filter merged candidates to human + adipose + menopause terms (priority=high/manual).
2) For top candidates, open the GEO/ArrayExpress record and verify: tissue, menopause status labeling, and availability of expression matrix/raw data.
3) Select 1–3 feasible datasets and implement a harmonized analysis in separate folders under `GTEx_v10_AT_analysis_out/tables/external_validation_<ACCESSION>/`.

Generated:
2026-02-01 01:10:42.383213

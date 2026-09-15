# Transcriptional lineage programs in KRAS/STK11 co-mutant lung adenocarcinoma

Analysis code for [paper title], [author list], [year].

This repository reproduces every figure, table, and statistic in the manuscript
from public data. No private or restricted data are required.

## Data sources

| Dataset | Source | Access |
|---|---|---|
| TCGA-LUAD expression and mutations | NCI Genomic Data Commons | Downloaded programmatically in `analysis_01` via TCGAbiolinks |
| GSE72094 | NCBI Gene Expression Omnibus | Downloaded programmatically in `analysis_01` |
| PRISM Repurposing Primary | DepMap [release version] | Manual download from depmap.org; see below |
| MSigDB Hallmark and C8 | MSigDB | Retrieved via `msigdbr` |

DepMap PRISM viability files must be downloaded manually from
[depmap.org](https://depmap.org/portal/) and placed in the working directory as:

```
DROPERIDOL PRISM Repurposing Primary (Viability).csv
MEBENDAZOLE PRISM Repurposing Primary (Viability).csv
ROSIGLITAZONE PRISM Repurposing Primary (Viability).csv
```

## Scripts

Run in order. Each writes intermediate files consumed by the next.

| Script | Purpose | Approx. runtime |
|---|---|---|
| `analysis_01_data_acquisition.R` | Download TCGA and GEO data, define cohorts | 30–45 min |
| `analysis_02_expression_processing.R` | Filter, transform, map gene symbols | 10 min |
| `analysis_03_enrichment_and_stats.R` | Hallmark screen, lineage scoring, all statistics | 20–30 min |
| `analysis_04_figures.R` | Generate all main and supplementary figures | 2 min |

## Requirements

R ≥ 4.4. Bioconductor packages: TCGAbiolinks, SummarizedExperiment, DESeq2,
GSVA, org.Hs.eg.db. CRAN: dplyr, msigdbr, ashr.

```r
install.packages(c("BiocManager", "dplyr", "msigdbr", "ashr"))
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "DESeq2",
                       "GSVA", "org.Hs.eg.db"))
```

Peak memory use is approximately 8 GB during the TCGA expression download and
GSVA scoring steps. `analysis_01` will fail on machines with less available RAM.

## Two implementation notes

**GEOquery fails on GSE72094.** `getGEO()` returns
`parsing failed--expected only one '!series_data_table_begin'` because the series
matrix uses the marker `!series_matrix_table_begin`. `analysis_01` parses the
file manually.

**GPL15048 has no usable annotation.** GSE72094 was profiled on a custom
Rosetta/Merck array with no GEO `.annot` file. Probe identifiers do not match
standard Affymetrix annotation — only 10 of 60,607 probes map via
`hgu133plus2.db`. Probe IDs embed source accessions (`merck-NM_002431_at`), so
`analysis_02` strips the prefix and suffix and maps the resulting RefSeq
accessions through `org.Hs.eg.db`, recovering 18,077 unique symbols.

## Gene signatures

**Neuroendocrine** — Zhang et al., *Transl Lung Cancer Res* 2018 (PMID 29535911),
Tables S1 and S2: 25 NE-associated and 25 non-NE-associated genes. Scored here as
ssGSEA(NE) − ssGSEA(non-NE), which approximates the published
(correl_NE − correl_non-NE)/2 formulation without requiring the original
reference expression profiles. This is a deviation from the published method and
is noted in the manuscript.

**Hepatocyte** — Aizarani et al., *Nature* 2019 human liver cell atlas, accessed
through MSigDB collection C8. All four hepatocyte clusters (C11, C14, C17, C30)
are scored independently.

## Analysis provenance

The neuroendocrine and hepatocyte findings were not pre-specified. Both emerged
during quality-control inspection of the differential expression results and were
subsequently tested formally using published signatures in both cohorts. This
sequence is stated in the manuscript.

## Citation

[Fill in on publication.]

Archived at [Zenodo DOI].

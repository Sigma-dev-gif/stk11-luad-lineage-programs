# Transcriptional lineage programs in KRAS/STK11 co-mutant lung adenocarcinoma

Analysis code accompanying:

> Minocha J. Neuroendocrine and hepatocyte transcriptional programs are independently
> elevated in STK11-mutant KRAS-driven lung adenocarcinoma. *Manuscript submitted.*

This repository reproduces every figure, table, and statistic in the manuscript from
publicly available data. No private or restricted data are required.

## Data sources

| Dataset | Source | Access |
|---|---|---|
| TCGA-LUAD expression and mutations | NCI Genomic Data Commons | Downloaded programmatically in `analysis_01` via TCGAbiolinks |
| GSE72094 | NCBI Gene Expression Omnibus | Downloaded programmatically in `analysis_01` |
| MSigDB Hallmark (H) and cell type (C8) | MSigDB | Retrieved via `msigdbr` |

## Scripts

Run in order. Each writes intermediate files consumed by the next.

| Script | Purpose | Approx. runtime |
|---|---|---|
| `analysis_01_data_acquisition.R` | Download TCGA and GEO data, define cohorts | 30–45 min |
| `analysis_02_expression_processing.R` | Filter, transform, map gene symbols | 10 min |
| `analysis_03_enrichment_and_stats.R` | Hallmark screen, lineage scoring, statistics | 20–30 min |
| `analysis_04_figures.R` | Main and supplementary figures | 2 min |
| `analysis_05_study_design_figure.R` | Study design schematic (Fig 1) | < 1 min |

## Requirements

R ≥ 4.4.

```r
install.packages(c("BiocManager", "dplyr", "msigdbr", "ashr"))
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "DESeq2",
                       "GSVA", "org.Hs.eg.db"))
```

Tested on R 4.6.1 with TCGAbiolinks 2.40.0, DESeq2 1.52.0, GSVA 2.6.6,
msigdbr 26.1.1, org.Hs.eg.db 3.23.1, ashr 2.2.63, dplyr 1.2.1.

Peak memory use is approximately 8 GB during the TCGA expression download and GSVA
scoring. `analysis_01` will fail on machines with less available RAM.

## Two implementation notes

**GEOquery fails on GSE72094.** `getGEO()` returns
`parsing failed--expected only one '!series_data_table_begin'` because the series
matrix uses the marker `!series_matrix_table_begin`. `analysis_01` parses the file
manually.

**GPL15048 has no usable annotation.** GSE72094 was profiled on a custom Rosetta/Merck
array with no GEO `.annot` file. Probe identifiers do not match standard Affymetrix
annotation — only 10 of 60,607 probes map via `hgu133plus2.db`. Probe IDs embed source
accessions (`merck-NM_002431_at`), so `analysis_02` strips the prefix and suffix and
maps the resulting RefSeq accessions through `org.Hs.eg.db`, recovering 18,077 unique
symbols.

## Gene signatures

**Neuroendocrine** — Zhang W, Girard L, Zhang YA, et al. Small cell lung cancer tumors
and preclinical models display heterogeneity of neuroendocrine phenotypes.
*Transl Lung Cancer Res.* 2018;7(1):32–49. PMID: 29535911. Tables S1 and S2 provide 25
NE-associated and 25 non-NE-associated genes. Scored here as
ssGSEA(NE) − ssGSEA(non-NE), which approximates the published
(correl_NE − correl_non-NE)/2 formulation without requiring the original reference
expression profiles. This deviation is stated in the manuscript.

**Hepatocyte** — Aizarani N, Saviano A, Sagar, et al. A human liver cell atlas reveals
heterogeneity and epithelial progenitors. *Nature.* 2019;572(7768):199–204. Accessed
through MSigDB collection C8. All four hepatocyte clusters (C11, C14, C17, C30) are
scored independently.

## Analysis provenance

The neuroendocrine and hepatocyte findings were not pre-specified. Both emerged during
quality-control inspection of the differential expression results and were subsequently
tested formally using published signatures in both cohorts. This sequence is stated in
the manuscript.

An exploratory drug-repurposing analysis using L1000 signature reversal and DepMap
PRISM data was performed during development and returned a null result. It is not part
of the submitted manuscript and the corresponding code is not included here.

## License

MIT. See `LICENSE`.

## Citation

Please cite the manuscript above. This repository is archived on Zenodo; the DOI will
be added on release.

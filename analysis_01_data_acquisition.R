####################################################################
# 01 — Data acquisition and cohort definition
#
# Downloads TCGA-LUAD mutation and expression data, defines the
# KRAS-mutant cohort stratified by STK11 status, and downloads and
# parses GSE72094.
#
# Outputs:
#   tcga_luad_rnaseq.rds
#   kras_subgroup_labels.csv
#   gse72094_expr_krasmut.rds
#   gse72094_pheno.rds
#
# Runtime: ~30-45 min, mostly download time. Requires ~8 GB RAM.
####################################################################

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

set.seed(42)

# ------------------------------------------------------------------
# 1. TCGA-LUAD somatic mutations
# ------------------------------------------------------------------
query_maf <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Simple Nucleotide Variation",
  data.type     = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query_maf)
maf_data <- GDCprepare(query_maf)

kras_ids  <- unique(maf_data$Tumor_Sample_Barcode[maf_data$Hugo_Symbol == "KRAS"])
stk11_ids <- unique(maf_data$Tumor_Sample_Barcode[maf_data$Hugo_Symbol == "STK11"])
keap1_ids <- unique(maf_data$Tumor_Sample_Barcode[maf_data$Hugo_Symbol == "KEAP1"])

message("KRAS-mutant samples:  ", length(kras_ids))
message("STK11-mutant samples: ", length(stk11_ids))
message("KEAP1-mutant samples: ", length(keap1_ids))

# ------------------------------------------------------------------
# 2. Define groups among KRAS-mutant tumors
#    Primary definition: STK11 status.
#    KEAP1 status is recorded for the secondary analysis but tumors
#    with KEAP1-only mutations remain in the comparator, mirroring
#    GSE72094 which does not annotate KEAP1.
# ------------------------------------------------------------------
subgroup <- data.frame(sample_id = unique(maf_data$Tumor_Sample_Barcode)) %>%
  mutate(
    KRAS_mut  = sample_id %in% kras_ids,
    STK11_mut = sample_id %in% stk11_ids,
    KEAP1_mut = sample_id %in% keap1_ids
  ) %>%
  filter(KRAS_mut) %>%
  mutate(
    group        = if_else(STK11_mut | KEAP1_mut, "KRAS_comutant", "KRAS_only"),
    group_stk11  = if_else(STK11_mut, "KRAS_STK11", "KRAS_only"),
    patient_id   = substr(sample_id, 1, 12)
  )

# ------------------------------------------------------------------
# 3. TCGA-LUAD expression, restricted to the KRAS-mutant patients
# ------------------------------------------------------------------
target_patients <- unique(subgroup$patient_id)

query_rna <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode       = target_patients
)
GDCdownload(query_rna)
rna_data <- GDCprepare(query_rna)

# Keep only subgroup entries that have expression data
subgroup_matched <- subgroup %>%
  filter(patient_id %in% substr(colnames(rna_data), 1, 12))

message("Final TCGA cohort, STK11 definition:")
print(table(subgroup_matched$group_stk11))

saveRDS(rna_data, "tcga_luad_rnaseq.rds")
write.csv(subgroup_matched, "kras_subgroup_labels.csv", row.names = FALSE)

# ------------------------------------------------------------------
# 4. GSE72094
#
# NOTE: GEOquery::getGEO() fails on this series with
#   "parsing failed--expected only one '!series_data_table_begin'".
# The file uses the marker '!series_matrix_table_begin'. The series
# matrix is therefore parsed manually below.
# ------------------------------------------------------------------
url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/",
              "GSE72094/matrix/GSE72094_series_matrix.txt.gz")
download.file(url, "GSE72094_matrix.txt.gz", mode = "wb")

lines <- readLines(gzfile("GSE72094_matrix.txt.gz"))
hdr   <- lines[seq_len(grep("^!series_matrix_table_begin", lines) - 1)]

getchar <- function(key) {
  ln <- hdr[grepl(paste0('"', key, ':'), hdr, fixed = TRUE)][1]
  v  <- strsplit(ln, "\t")[[1]][-1]
  trimws(sub(paste0("^", key, ":"), "", gsub('^"|"$', '', v)))
}

sample_ids <- gsub('^"|"$', '',
  strsplit(lines[grep("^!Sample_geo_accession", lines)], "\t")[[1]][-1])

pheno_B <- data.frame(
  sample_id = sample_ids,
  kras      = getchar("kras_status"),
  stk11     = getchar("stk11_status"),
  stringsAsFactors = FALSE
)

message("GSE72094 KRAS x STK11:")
print(table(pheno_B$kras, pheno_B$stk11))

tbl_start <- grep("^!series_matrix_table_begin", lines)
tbl_end   <- grep("^!series_matrix_table_end", lines)
dat <- read.delim(
  text = paste(lines[(tbl_start + 1):(tbl_end - 1)], collapse = "\n"),
  header = TRUE, check.names = FALSE, stringsAsFactors = FALSE
)

# Restrict to KRAS-mutant samples
keep_cols <- c("ID_REF", pheno_B$sample_id[pheno_B$kras == "Mut"])
dat <- dat[, intersect(keep_cols, colnames(dat))]

saveRDS(dat,     "gse72094_expr_krasmut.rds")
saveRDS(pheno_B, "gse72094_pheno.rds")

rm(lines, hdr); gc()

message("Done. Session info:")
sessionInfo()

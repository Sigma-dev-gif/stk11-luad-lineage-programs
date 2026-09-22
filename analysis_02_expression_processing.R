####################################################################
# 02 — Expression processing and gene symbol mapping
#
# Builds analysis-ready, symbol-indexed expression matrices for both
# cohorts.
#
# Inputs:  tcga_luad_rnaseq.rds, kras_subgroup_labels.csv,
#          gse72094_expr_krasmut.rds, gse72094_pheno.rds
# Outputs: tcga_stk11_expr_matrix.rds, tcga_stk11_coldata.rds,
#          gse72094_expr_matrix.rds
####################################################################

library(DESeq2)
library(SummarizedExperiment)
library(org.Hs.eg.db)   # retained for reference; not used for GSE72094
library(dplyr)

set.seed(42)

# Collapse duplicate gene symbols by keeping the most variable row
collapse_symbols <- function(mat, symbols) {
  keep   <- !is.na(symbols) & symbols != ""
  mat    <- mat[keep, ]; symbols <- symbols[keep]
  ord    <- order(apply(mat, 1, var, na.rm = TRUE), decreasing = TRUE)
  mat    <- mat[ord, ]; symbols <- symbols[ord]
  mat    <- mat[!duplicated(symbols), ]
  rownames(mat) <- symbols[!duplicated(symbols)]
  mat
}

# ------------------------------------------------------------------
# 1. TCGA — one sample per patient, STK11 group labels
# ------------------------------------------------------------------
rna_data <- readRDS("tcga_luad_rnaseq.rds")
sub_all  <- read.csv("kras_subgroup_labels.csv", stringsAsFactors = FALSE)
sub_uniq <- distinct(sub_all, patient_id, .keep_all = TRUE)

coldata <- data.frame(
  sample_barcode = colnames(rna_data),
  patient_id     = substr(colnames(rna_data), 1, 12),
  stringsAsFactors = FALSE
) %>%
  left_join(select(sub_uniq, patient_id, group_stk11, STK11_mut, KEAP1_mut),
            by = "patient_id")

# TCGA contains repeat samples for some patients; retain the first
if (any(duplicated(coldata$patient_id))) {
  coldata  <- coldata[!duplicated(coldata$patient_id), ]
  rna_data <- rna_data[, coldata$sample_barcode]
}

coldata$group <- factor(coldata$group_stk11,
                        levels = c("KRAS_only", "KRAS_STK11"))
colData(rna_data)$group <- coldata$group

# Filter low-count genes and variance-stabilise
dds <- DESeqDataSet(rna_data, design = ~ group)
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]
dds <- estimateSizeFactors(dds)

vsd    <- vst(dds, blind = TRUE)
expr_A <- assay(vsd)

gm  <- as.data.frame(rowData(dds))
sym <- gm$gene_name[match(rownames(expr_A), gm$gene_id)]
expr_A <- collapse_symbols(expr_A, sym)

message("TCGA matrix: ", nrow(expr_A), " genes x ", ncol(expr_A), " samples")
print(table(coldata$group))

saveRDS(expr_A,  "tcga_stk11_expr_matrix.rds")
saveRDS(coldata, "tcga_stk11_coldata.rds")

rm(rna_data, dds, vsd); gc()

# ------------------------------------------------------------------
# 2. GSE72094 — probe annotation from the canonical platform record
#
# GPL15048 is a custom Rosetta/Merck array. Probe IDs do not match standard
# Affymetrix annotation packages (only 10 of 60,607 probes map via
# hgu133plus2.db), but the platform record itself carries a full annotation
# table with GenBank IDs, Entrez gene IDs and HGNC symbols.
#
# Note: getGEO("GSE72094") fails on the *series* matrix (see analysis_01),
# but getGEO("GPL15048") retrieves the *platform* record normally.
#
# An earlier version of this script parsed RefSeq accessions out of the
# probe names instead, yielding 18,077 genes. That was replaced after
# Dr. Steven Eschrich (Moffitt Cancer Center) pointed out the canonical
# annotation is available. The canonical table gives 22,115 genes.
# ------------------------------------------------------------------
library(GEOquery)

dat     <- readRDS("gse72094_expr_krasmut.rds")
pheno_B <- readRDS("gse72094_pheno.rds")

gpl <- getGEO("GPL15048")
ann <- Table(gpl)
ann$GeneSymbol <- trimws(as.character(ann$GeneSymbol))
ann <- ann[!is.na(ann$GeneSymbol) & ann$GeneSymbol != "", c("ID", "GeneSymbol")]

message("GPL15048: ", nrow(ann), " probes carry a symbol, ",
        length(unique(ann$GeneSymbol)), " unique genes")

m      <- merge(data.frame(ID = dat$ID_REF, stringsAsFactors = FALSE), ann, by = "ID")
expr_B <- as.matrix(dat[match(m$ID, dat$ID_REF), -1])
expr_B <- collapse_symbols(expr_B, m$GeneSymbol)

# Restrict to KRAS-mutant samples
expr_B <- expr_B[, intersect(pheno_B$sample_id[pheno_B$kras == "Mut"],
                             colnames(expr_B))]

message("GSE72094 matrix: ", nrow(expr_B), " genes x ", ncol(expr_B), " samples")
message("Value range: ", paste(round(range(expr_B, na.rm = TRUE), 2), collapse = " - "))

saveRDS(expr_B, "gse72094_expr_matrix.rds")

sessionInfo()

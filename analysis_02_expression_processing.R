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
library(org.Hs.eg.db)
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
# 2. GSE72094 — probe annotation
#
# GPL15048 is a custom Rosetta/Merck array. No GEO .annot file exists
# and probe IDs do not match standard Affymetrix annotation (only 10
# of 60,607 probes map via hgu133plus2.db). Probe IDs embed source
# accessions, e.g. "merck-NM_002431_at". Symbols are recovered by
# stripping the prefix and suffix and mapping RefSeq accessions.
# ------------------------------------------------------------------
dat     <- readRDS("gse72094_expr_krasmut.rds")
pheno_B <- readRDS("gse72094_pheno.rds")

acc <- dat$ID_REF
acc <- sub("^merck2?-", "", acc)
acc <- sub("_[a-z]?_?at$", "", acc)
acc <- sub("\\..*$", "", acc)

pm <- data.frame(probe = dat$ID_REF, acc = acc, stringsAsFactors = FALSE)
rs <- pm[grepl("^[NX][MR]_", pm$acc), ]

sm <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(rs$acc),
                            columns = "SYMBOL", keytype = "REFSEQ")
sm <- sm[!is.na(sm$SYMBOL), ]
sm <- sm[!duplicated(sm$REFSEQ), ]
rs <- merge(rs, sm, by.x = "acc", by.y = "REFSEQ")

message("GSE72094: ", nrow(rs), " probes mapped to ",
        length(unique(rs$SYMBOL)), " unique symbols")

expr_B <- as.matrix(dat[match(rs$probe, dat$ID_REF), -1])
expr_B <- collapse_symbols(expr_B, rs$SYMBOL)

# Restrict to KRAS-mutant samples
expr_B <- expr_B[, intersect(pheno_B$sample_id[pheno_B$kras == "Mut"],
                             colnames(expr_B))]

message("GSE72094 matrix: ", nrow(expr_B), " genes x ", ncol(expr_B), " samples")
message("Value range: ", paste(round(range(expr_B, na.rm = TRUE), 2), collapse = " - "))

saveRDS(expr_B, "gse72094_expr_matrix.rds")

sessionInfo()

####################################################################
# 03 — Gene set enrichment, lineage scoring, and statistics
#
# Runs the Hallmark screen in both cohorts, scores published
# neuroendocrine and hepatocyte signatures, tests cross-cohort
# concordance, tests program independence, and runs the
# neuroendocrine-exclusion sensitivity analysis.
#
# Inputs:  tcga_stk11_expr_matrix.rds, tcga_stk11_coldata.rds,
#          gse72094_expr_matrix.rds, gse72094_pheno.rds
# Outputs: hallmark results, lineage scores, concordance table,
#          sensitivity analysis
####################################################################

library(GSVA)
library(msigdbr)
library(dplyr)

set.seed(42)

expr_A  <- readRDS("tcga_stk11_expr_matrix.rds")
coldata <- readRDS("tcga_stk11_coldata.rds")
expr_B  <- readRDS("gse72094_expr_matrix.rds")
pheno_B <- readRDS("gse72094_pheno.rds")

grpA <- coldata$group[match(colnames(expr_A), coldata$sample_barcode)]
grpB <- pheno_B$stk11[match(colnames(expr_B), pheno_B$sample_id)]
grpB <- factor(ifelse(grpB == "Mut", "KRAS_STK11", "KRAS_only"),
               levels = c("KRAS_only", "KRAS_STK11"))

# Rank-biserial correlation from a Wilcoxon test
effect_r <- function(x, g) {
  w  <- wilcox.test(x ~ g)$statistic
  n1 <- sum(g == levels(g)[1]); n2 <- sum(g == levels(g)[2])
  as.numeric(1 - (2 * w) / (n1 * n2))
}

test_sets <- function(scores, g) {
  data.frame(
    gene_set = rownames(scores),
    p        = apply(scores, 1, function(x) wilcox.test(x ~ g)$p.value),
    effect   = apply(scores, 1, effect_r, g = g),
    stringsAsFactors = FALSE
  ) %>% mutate(FDR = p.adjust(p, method = "BH"))
}

# ------------------------------------------------------------------
# 1. Hallmark screen, both cohorts
# ------------------------------------------------------------------
h        <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- split(h$gene_symbol, h$gs_name)

scores_A <- gsva(ssgseaParam(expr_A, hallmark))
scores_B <- gsva(ssgseaParam(expr_B, hallmark))

resA <- test_sets(scores_A, grpA)
resB <- test_sets(scores_B, grpB)

message("Hallmark sets at FDR<0.05 — TCGA: ", sum(resA$FDR < 0.05),
        ", GSE72094: ", sum(resB$FDR < 0.05))

# ------------------------------------------------------------------
# 2. Cross-cohort concordance
# ------------------------------------------------------------------
comp <- inner_join(
  select(resA, gene_set, effect_TCGA = effect, FDR_TCGA = FDR),
  select(resB, gene_set, effect_GEO  = effect, FDR_GEO  = FDR),
  by = "gene_set"
) %>%
  mutate(replicated = FDR_TCGA < 0.05 & FDR_GEO < 0.05 &
                      sign(effect_TCGA) == sign(effect_GEO))

ct <- cor.test(comp$effect_TCGA, comp$effect_GEO, method = "spearman")
message("Cross-cohort Spearman rho = ", round(ct$estimate, 3),
        ", replicated sets = ", sum(comp$replicated))

write.csv(comp, "phase5_cohort_concordance.csv", row.names = FALSE)

# ------------------------------------------------------------------
# 3. Published lineage signatures
#
# Neuroendocrine: Zhang et al., Transl Lung Cancer Res 2018
#   (PMID 29535911), Tables S1 and S2 — 25 NE and 25 non-NE genes.
#   Scored as ssGSEA(NE) - ssGSEA(non-NE), approximating the
#   published (correl_NE - correl_nonNE)/2 formulation without
#   requiring the original reference expression profiles.
#
# Hepatocyte: Aizarani et al., Nature 2019 human liver cell atlas,
#   via MSigDB C8. Four hepatocyte clusters scored independently.
# ------------------------------------------------------------------
ne_zhang <- c("BEX1","ASCL1","INSM1","CHGA","TAGLN3","KIF5C","CRMP1","SCG3",
              "SYT4","RTN1","MYT1","SYP","KIF1A","TMSB15A","SYN1","SYT11",
              "RUNDC3A","TFF3","CHGB","FAM57B","SH3GL2","BSN","SEZ6",
              "TMSB15B","CELF3")

nonne_zhang <- c("RAB27B","TGFBR2","SLC16A5","S100A10","ITGB4","YAP1","LGALS3",
                 "EPHA2","S100A16","PLAU","ABCC3","ARHGDIB","CYR61","PTGES",
                 "CCND1","IFITM2","IFITM3","AHNAK","CAV2","TACSTD2","TGFBI",
                 "EMP1","CAV1","ANXA1","MYOF")

ne_score <- function(expr) {
  s <- gsva(ssgseaParam(expr, list(NE = ne_zhang, nonNE = nonne_zhang)))
  as.numeric(s["NE", ] - s["nonNE", ])
}

nez_A <- ne_score(expr_A)
nez_B <- ne_score(expr_B)

h_all      <- msigdbr(species = "Homo sapiens")
liver_sets <- c("AIZARANI_LIVER_C11_HEPATOCYTES_1",
                "AIZARANI_LIVER_C14_HEPATOCYTES_2",
                "AIZARANI_LIVER_C17_HEPATOCYTES_3",
                "AIZARANI_LIVER_C30_HEPATOCYTES_4")
hep_sets <- lapply(setNames(liver_sets, liver_sets), function(s)
  unique(h_all$gene_symbol[h_all$gs_name == s]))

hepA <- gsva(ssgseaParam(expr_A, hep_sets))
hepB <- gsva(ssgseaParam(expr_B, hep_sets))

lineage <- rbind(
  data.frame(signature = "NE_Zhang", cohort = "TCGA",
             p = wilcox.test(nez_A ~ grpA)$p.value, r = effect_r(nez_A, grpA)),
  data.frame(signature = "NE_Zhang", cohort = "GSE72094",
             p = wilcox.test(nez_B ~ grpB)$p.value, r = effect_r(nez_B, grpB)),
  do.call(rbind, lapply(liver_sets, function(s) rbind(
    data.frame(signature = s, cohort = "TCGA",
               p = wilcox.test(as.numeric(hepA[s,]) ~ grpA)$p.value,
               r = effect_r(as.numeric(hepA[s,]), grpA)),
    data.frame(signature = s, cohort = "GSE72094",
               p = wilcox.test(as.numeric(hepB[s,]) ~ grpB)$p.value,
               r = effect_r(as.numeric(hepB[s,]), grpB))
  )))
)
print(lineage, row.names = FALSE)
write.csv(lineage, "table2_lineage_programs.csv", row.names = FALSE)

saveRDS(list(nez_A = nez_A, nez_B = nez_B, hepA = hepA, hepB = hepB,
             grpA = grpA, grpB = grpB), "published_signature_scores.rds")

# ------------------------------------------------------------------
# 4. Are the two programs independently associated with genotype?
# ------------------------------------------------------------------
hepA_c11 <- as.numeric(hepA["AIZARANI_LIVER_C11_HEPATOCYTES_1", ])
hepB_c11 <- as.numeric(hepB["AIZARANI_LIVER_C11_HEPATOCYTES_1", ])

message("\nLogistic regression, TCGA:")
print(summary(glm(as.integer(grpA == "KRAS_STK11") ~
                  scale(nez_A) + scale(hepA_c11),
                  family = binomial))$coefficients)

message("\nLogistic regression, GSE72094:")
print(summary(glm(as.integer(grpB == "KRAS_STK11") ~
                  scale(nez_B) + scale(hepB_c11),
                  family = binomial))$coefficients)

message("\nNE-hepatocyte correlation — TCGA: ",
        round(cor(nez_A, hepA_c11, method = "spearman"), 3),
        ", GSE72094: ", round(cor(nez_B, hepB_c11, method = "spearman"), 3))

# ------------------------------------------------------------------
# 5. Sensitivity analysis: exclude neuroendocrine-high tumours
#
# A percentile threshold is used rather than an absolute expression
# cutoff so the same criterion applies across platforms.
# ------------------------------------------------------------------
sensitivity <- function(scores, ne, g, label) {
  keep <- ne <= quantile(ne, 0.90)
  r    <- test_sets(scores[, keep], droplevels(g[keep]))
  full <- test_sets(scores, g)
  rho  <- cor(full$effect, r$effect, method = "spearman")
  message(label, " — retained ", sum(keep), "/", length(keep),
          ", significant ", sum(r$FDR < 0.05), "/50, rho = ", round(rho, 3))
  r
}

resA_ne <- sensitivity(scores_A, nez_A, grpA, "TCGA")
resB_ne <- sensitivity(scores_B, nez_B, grpB, "GSE72094")

write.csv(resA_ne, "tcga_NEexcluded_sensitivity.csv",     row.names = FALSE)
write.csv(resB_ne, "gse72094_NEexcluded_sensitivity.csv", row.names = FALSE)

sessionInfo()

####################################################################
# 04 — Figures
#
# Generates all main and supplementary figures.
#
# Inputs:  published_signature_scores.rds, phase5_cohort_concordance.csv,
#          tcga_stk11_expr_matrix.rds, tcga_stk11_coldata.rds,
#          tcga_NEexcluded_sensitivity.csv, gse72094_NEexcluded_sensitivity.csv,
#          three DepMap PRISM CSVs (manual download — see README)
#
# Outputs: fig_lineage_scores.png, fig_crosscohort_concordance.png,
#          fig_hallmark_screen.png, fig_depmap_null.png,
#          figS1_ne_sensitivity.png, figS2_ascl1_bimodality.png,
#          figS3_ne_hep_scatter.png
####################################################################

library(dplyr)

ps       <- readRDS("published_signature_scores.rds")
comp     <- read.csv("phase5_cohort_concordance.csv")
hepA_c11 <- as.numeric(ps$hepA["AIZARANI_LIVER_C11_HEPATOCYTES_1", ])
hepB_c11 <- as.numeric(ps$hepB["AIZARANI_LIVER_C11_HEPATOCYTES_1", ])

wilcox_p <- function(x, g) format.pval(wilcox.test(x ~ g)$p.value, digits = 2)
eff_r <- function(x, g) {
  w  <- wilcox.test(x ~ g)$statistic
  n1 <- sum(g == levels(g)[1]); n2 <- sum(g == levels(g)[2])
  as.numeric(1 - (2 * w) / (n1 * n2))
}

# ------------------------------------------------------------------
# Figure 1 — lineage program scores by genotype
# ------------------------------------------------------------------
png("fig_lineage_scores.png", width = 2000, height = 1800, res = 220)
par(mfrow = c(2, 2), mar = c(4, 4.5, 3.5, 1.5))

bx <- function(score, grp, ttl, yl) {
  boxplot(score ~ grp, col = c("grey85", "firebrick"), outline = FALSE,
          names = c("KRAS-only", "KRAS/STK11"), ylab = yl, xlab = "",
          main = ttl, border = "grey25")
  stripchart(score ~ grp, vertical = TRUE, method = "jitter", jitter = 0.15,
             pch = 21, bg = adjustcolor("black", 0.25), col = "grey40",
             add = TRUE, cex = 0.7)
  mtext(sprintf("p = %s, r = %.2f", wilcox_p(score, grp), eff_r(score, grp)),
        side = 3, line = 0.2, cex = 0.75)
}

bx(ps$nez_A,  ps$grpA, "Neuroendocrine (Zhang) — TCGA",        "NE score")
bx(ps$nez_B,  ps$grpB, "Neuroendocrine (Zhang) — GSE72094",    "NE score")
bx(hepA_c11,  ps$grpA, "Hepatocyte (Aizarani C11) — TCGA",     "ssGSEA score")
bx(hepB_c11,  ps$grpB, "Hepatocyte (Aizarani C11) — GSE72094", "ssGSEA score")
dev.off()

# ------------------------------------------------------------------
# Figure 2 — cross-cohort concordance of Hallmark effect sizes
# ------------------------------------------------------------------
png("fig_crosscohort_concordance.png", width = 1800, height = 1600, res = 220)
par(mar = c(5, 5, 4, 2))
sig_both <- comp$FDR_TCGA < 0.05 & comp$FDR_GEO < 0.05
plot(comp$effect_TCGA, comp$effect_GEO, pch = 21, cex = 1.4,
     bg = ifelse(sig_both, "firebrick", "grey80"), col = "grey30",
     xlab = "Effect size, TCGA (rank-biserial r)",
     ylab = "Effect size, GSE72094 (rank-biserial r)",
     main = "Hallmark effect sizes replicate across cohorts")
abline(h = 0, v = 0, col = "grey70", lty = 2)
abline(lm(effect_GEO ~ effect_TCGA, data = comp), col = "steelblue", lwd = 2)
ct <- cor.test(comp$effect_TCGA, comp$effect_GEO, method = "spearman")
legend("topleft", bty = "n",
       legend = c(sprintf("Spearman rho = %.3f", ct$estimate),
                  sprintf("p %s", format.pval(ct$p.value, digits = 3)),
                  sprintf("%d of 50 significant in both", sum(sig_both))))
dev.off()

# ------------------------------------------------------------------
# Figure 3 — Hallmark sets altered in KRAS/STK11 tumours (TCGA)
# ------------------------------------------------------------------
sig <- comp[comp$FDR_TCGA < 0.05, ]
sig <- sig[order(sig$effect_TCGA), ]
sig$label <- gsub("_", " ", gsub("HALLMARK_", "", sig$gene_set))

png("fig_hallmark_screen.png", width = 1600, height = 1400, res = 200)
par(mar = c(5, 18, 4, 2))
barplot(sig$effect_TCGA, horiz = TRUE, names.arg = sig$label, las = 1,
        col = ifelse(sig$effect_TCGA > 0, "firebrick", "steelblue"),
        border = NA, xlim = c(-0.5, 0.6), cex.names = 0.75,
        xlab = "Effect size (rank-biserial r)",
        main = "Hallmark sets altered in KRAS/STK11 tumours", cex.main = 1.1)
abline(v = 0, col = "grey40")
legend("bottomright", bty = "n", fill = c("firebrick", "steelblue"),
       legend = c("Up in STK11-mutant", "Down in STK11-mutant"), cex = 0.8)
dev.off()

# ------------------------------------------------------------------
# Figure 4 — DepMap null result
#
# Requires the three PRISM CSVs; see README.
# ------------------------------------------------------------------
prism_files <- c(
  Droperidol    = "DROPERIDOL PRISM Repurposing Primary (Viability).csv",
  Mebendazole   = "MEBENDAZOLE PRISM Repurposing Primary (Viability).csv",
  Rosiglitazone = "ROSIGLITAZONE PRISM Repurposing Primary (Viability).csv"
)

if (all(file.exists(prism_files))) {
  png("fig_depmap_null.png", width = 2000, height = 800, res = 200)
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1))
  for (nm in names(prism_files)) {
    df <- read.csv(prism_files[[nm]])
    names(df)[1:6] <- c("id", "viab", "cell", "dis", "lin", "sub")
    lung  <- df$viab[df$lin == "Lung"]
    other <- df$viab[df$lin != "Lung"]
    boxplot(list(Other = other, Lung = lung),
            col = c("grey85", "steelblue"), outline = FALSE, border = "grey25",
            ylab = "PRISM viability (log2 FC)", main = nm)
    stripchart(list(Other = other, Lung = lung), vertical = TRUE,
               method = "jitter", jitter = 0.15, pch = 21,
               bg = adjustcolor("black", 0.15), col = "grey50",
               add = TRUE, cex = 0.5)
    abline(h = 0, lty = 2, col = "grey60")
    mtext(sprintf("Mann-Whitney p = %.2f",
                  wilcox.test(lung, other)$p.value), side = 3, line = 0.2, cex = 0.75)
  }
  dev.off()
} else {
  warning("PRISM files not found; skipping Figure 4. See README.")
}

# ------------------------------------------------------------------
# Figure S1 — neuroendocrine-exclusion sensitivity
# ------------------------------------------------------------------
resA_ne <- read.csv("tcga_NEexcluded_sensitivity.csv")
resB_ne <- read.csv("gse72094_NEexcluded_sensitivity.csv")

png("figS1_ne_sensitivity.png", width = 1800, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))

sens_panel <- function(full_eff, excl_eff, ttl) {
  plot(full_eff, excl_eff, pch = 21, bg = "grey70", col = "grey30", cex = 1.2,
       xlab = "Effect size, all samples",
       ylab = "Effect size, NE-high excluded", main = ttl)
  abline(0, 1, lty = 2, col = "grey50"); abline(h = 0, v = 0, col = "grey85")
  legend("topleft", bty = "n",
         legend = sprintf("rho = %.3f",
                          cor(full_eff, excl_eff, method = "spearman")))
}

cA <- merge(comp[, c("gene_set", "effect_TCGA")], resA_ne[, c("gene_set", "effect")],
            by = "gene_set")
cB <- merge(comp[, c("gene_set", "effect_GEO")], resB_ne[, c("gene_set", "effect")],
            by = "gene_set")
sens_panel(cA$effect_TCGA, cA$effect, "TCGA")
sens_panel(cB$effect_GEO,  cB$effect, "GSE72094")
dev.off()

# ------------------------------------------------------------------
# Figure S2 — ASCL1 bimodality motivating the sensitivity analysis
# ------------------------------------------------------------------
expr_A  <- readRDS("tcga_stk11_expr_matrix.rds")
coldata <- readRDS("tcga_stk11_coldata.rds")

png("figS2_ascl1_bimodality.png", width = 1600, height = 800, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))

a <- expr_A["ASCL1", ]
hist(a, breaks = 40, col = "grey80", border = "white",
     xlab = "ASCL1 expression (VST)", main = "ASCL1 distribution, TCGA")
abline(v = quantile(a, 0.90), col = "firebrick", lwd = 2, lty = 2)
legend("topright", bty = "n", legend = "90th percentile",
       col = "firebrick", lty = 2, lwd = 2, cex = 0.8)

g <- coldata$group[match(colnames(expr_A), coldata$sample_barcode)]
boxplot(a ~ g, col = c("grey85", "firebrick"), border = "grey25", outline = FALSE,
        names = c("KRAS-only", "KRAS/STK11"), xlab = "",
        ylab = "ASCL1 expression (VST)", main = "ASCL1 by genotype")
stripchart(a ~ g, vertical = TRUE, method = "jitter", jitter = 0.15, pch = 21,
           bg = adjustcolor("black", 0.3), col = "grey40", add = TRUE, cex = 0.7)
dev.off()

# ------------------------------------------------------------------
# Figure S3 — neuroendocrine vs hepatocyte scores
# ------------------------------------------------------------------
png("figS3_ne_hep_scatter.png", width = 1800, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))

sc <- function(ne, hep, grp, ttl) {
  plot(ne, hep, pch = 21, cex = 1.3,
       bg = ifelse(grp == "KRAS_STK11", adjustcolor("firebrick", 0.8),
                   adjustcolor("grey70", 0.7)),
       col = "grey30", xlab = "NE score (Zhang)",
       ylab = "Hepatocyte score (Aizarani C11)", main = ttl)
  abline(h = median(hep), v = median(ne), lty = 2, col = "grey60")
  legend("topright", bty = "n", cex = 0.8,
         legend = sprintf("rho = %.3f", cor(ne, hep, method = "spearman")))
  legend("bottomright", bty = "n", cex = 0.75, pch = 21,
         pt.bg = c("grey70", "firebrick"),
         legend = c("KRAS-only", "KRAS/STK11"))
}

sc(ps$nez_A, hepA_c11, ps$grpA, "TCGA")
sc(ps$nez_B, hepB_c11, ps$grpB, "GSE72094")
dev.off()

message("All figures written.")
sessionInfo()

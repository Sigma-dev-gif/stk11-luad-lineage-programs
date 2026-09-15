####################################################################
# Study design schematic (draft)
#
# A blueprint for the final figure. Base R graphics will not match a
# BioRender or Inkscape version aesthetically — this exists to fix the
# layout, wording and numbers before rebuilding it properly.
#
# Output: fig1_study_design.png
####################################################################

png("fig1_study_design.png", width = 2200, height = 2600, res = 220)
par(mar = c(0, 0, 0, 0))
plot(NA, xlim = c(0, 100), ylim = c(0, 118), axes = FALSE, xlab = "", ylab = "")

box_fill   <- "#EEF2F7"
box_border <- "#4A6785"
accent     <- "#B23A48"

draw_box <- function(x, y, w, h, label, fill = box_fill,
                     border = box_border, cex = 0.72, font = 1) {
  rect(x - w/2, y - h/2, x + w/2, y + h/2,
       col = fill, border = border, lwd = 1.6)
  text(x, y, label, cex = cex, font = font)
}

arrow_down <- function(x, y0, y1, lab = NULL) {
  arrows(x, y0, x, y1, length = 0.07, lwd = 1.6, col = "grey35")
  if (!is.null(lab)) text(x + 1.2, (y0 + y1)/2, lab, cex = 0.6,
                          col = "grey30", adj = 0)
}

panel_label <- function(x, y, txt) {
  text(x, y, txt, cex = 1.05, font = 2, adj = 0)
}

# ---------------- Panel A: cohorts ----------------
panel_label(3, 115, "A   Cohorts")

draw_box(27, 108, 40, 7,
         "TCGA-LUAD\nRNA-seq (STAR counts)", font = 2)
draw_box(73, 108, 40, 7,
         "GSE72094\nRosetta/Merck array (GPL15048)", font = 2)

arrow_down(27, 104.5, 100.5); arrow_down(73, 104.5, 100.5)

draw_box(27, 97, 40, 6, "Restrict to KRAS-mutant\nn = 140")
draw_box(73, 97, 40, 6, "Restrict to KRAS-mutant\nn = 154")

# fork from each cohort box into its two groups
segments(27, 94, 27, 92, col = "grey35", lwd = 1.6)
segments(16, 92, 38, 92, col = "grey35", lwd = 1.6)
arrows(16, 92, 16, 89.7, length = 0.06, lwd = 1.6, col = "grey35")
arrows(38, 92, 38, 89.7, length = 0.06, lwd = 1.6, col = "grey35")

segments(73, 94, 73, 92, col = "grey35", lwd = 1.6)
segments(62, 92, 84, 92, col = "grey35", lwd = 1.6)
arrows(62, 92, 62, 89.7, length = 0.06, lwd = 1.6, col = "grey35")
arrows(84, 92, 84, 89.7, length = 0.06, lwd = 1.6, col = "grey35")

draw_box(16, 86, 19, 7, "KRAS/STK11\nn = 33", fill = "#F5DEDE", border = accent)
draw_box(38, 86, 19, 7, "KRAS-only\nn = 107")
draw_box(62, 86, 19, 7, "KRAS/STK11\nn = 39", fill = "#F5DEDE", border = accent)
draw_box(84, 86, 19, 7, "KRAS-only\nn = 115")

text(50, 78.5,
     "KEAP1 is not annotated in GSE72094. Groups are matched on STK11 across cohorts;\nKEAP1-only mutants remain in the comparator, making the primary comparison conservative.",
     cex = 0.6, col = "grey25", font = 3)

# ---------------- Panel B: analysis ----------------
panel_label(3, 73, "B   Analysis")

draw_box(50, 68, 46, 5.5,
         "Symbol-mapped expression matrices\n27,428 genes (TCGA) | 18,077 genes (GSE72094)")
arrow_down(50, 65.2, 62)

draw_box(50, 58.5, 30, 5, "ssGSEA per-sample scoring", font = 2)

segments(20, 55.5, 80, 55.5, col = "grey55", lwd = 1.4)
segments(20, 55.5, 20, 52.5, col = "grey55", lwd = 1.4)
segments(50, 55.5, 50, 52.5, col = "grey55", lwd = 1.4)
segments(80, 55.5, 80, 52.5, col = "grey55", lwd = 1.4)
arrows(20, 52.5, 20, 51, length = 0.06, lwd = 1.5, col = "grey35")
arrows(50, 52.5, 50, 51, length = 0.06, lwd = 1.5, col = "grey35")
arrows(80, 52.5, 80, 51, length = 0.06, lwd = 1.5, col = "grey35")

draw_box(20, 46.5, 26, 8,
         "50 MSigDB Hallmark\ngene sets\n(Liberzon 2015)")
draw_box(50, 46.5, 26, 8,
         "Neuroendocrine\n50-gene signature\n(Zhang 2018)")
draw_box(80, 46.5, 26, 8,
         "Hepatocyte clusters\nC11/C14/C17/C30\n(Aizarani 2019)")

arrow_down(20, 42.2, 38.5)

# NE and hepatocyte branches converge on the shared testing box
segments(50, 42.2, 50, 40.5, col = "grey35", lwd = 1.6)
segments(80, 42.2, 80, 40.5, col = "grey35", lwd = 1.6)
segments(50, 40.5, 80, 40.5, col = "grey35", lwd = 1.6)
arrows(65, 40.5, 65, 38.7, length = 0.06, lwd = 1.6, col = "grey35")

draw_box(20, 34.5, 26, 8,
         "Wilcoxon + BH-FDR\nrank-biserial effect size\nper cohort")
draw_box(65, 34.5, 40, 8,
         "Wilcoxon per cohort\nJoint logistic regression:\nSTK11 status ~ NE + hepatocyte")

arrow_down(20, 30.2, 26.5); arrow_down(65, 30.2, 26.5)

draw_box(20, 22.5, 30, 7,
         "Cross-cohort concordance\nSpearman rho = 0.890",
         fill = "#F5DEDE", border = accent, font = 2)
draw_box(65, 22.5, 40, 7,
         "Both programs independently\nassociated with STK11 status",
         fill = "#F5DEDE", border = accent, font = 2)

# ---------------- Panel C: sensitivity ----------------
panel_label(3, 15, "C   Sensitivity analysis")

draw_box(28, 8.5, 44, 7,
         "Exclude top decile of neuroendocrine score\nwithin each cohort, repeat Hallmark screen")
arrows(51, 8.5, 57, 8.5, length = 0.07, lwd = 1.6, col = "grey35")
draw_box(78, 8.5, 40, 7,
         "Effect sizes preserved\nrho = 0.968 (TCGA), 0.960 (GSE72094)",
         fill = "#F5DEDE", border = accent)

dev.off()
cat("saved fig1_study_design.png\n")

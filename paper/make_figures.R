# Figures 1, 2 and S1 of Guo and Zhu, computed from the evaluated simulation results.
# Run from the repository root:  Rscript paper/make_figures.R
# Reads results/summary_main.csv and results/summary_size.csv and writes PDFs to paper/figures/. It also prints the
# numbers quoted in the text. Needs the packages ggplot2 and patchwork.
#
# Definitions follow make_tables.R (main grid, no outcome heterogeneity unless stated):
#   relative MSE       = MSE of a method / MSE of the RCT-only analysis in the same scenario and c
#   excess source bias = mean over the three sources of the bias of a source's weighted mean minus that of the oracle
#                        (same learner, scenario); taken at c = 32 because the weights do not depend on c
#   MI                 = MI with the default engine
#   relative weight    = weight share / subject share of incomplete external subjects (0 = complete-case weighting)

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
out_dir = "paper/figures"
dir.create(out_dir, showWarnings = FALSE)

d0 = read.csv("results/summary_main.csv")
rct = d0[d0$method == "rct_only", c("sid", "c", "MSE")]; names(rct)[3] = "MSE_rct"
D = merge(d0[d0$method != "rct_only", ], rct, by = c("sid", "c")); D$rel = D$MSE / D$MSE_rct
D$learner = sub("_.*", "", D$method); D$strategy = sub("^[a-z]+_", "", D$method)
D$sb = rowMeans(D[, paste0("bias_src_", 1:3)])
o = D[D$strategy == "oracle", c("sid", "c", "learner", "sb")]; names(o)[4] = "sb_or"
D = merge(D, o, by = c("sid", "c", "learner")); D$ex = D$sb - D$sb_or
M = D[D$grid == "main", ]
rates = c(0.1, 0.2, 0.3, 0.4, 0.5)
okabe = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7")
theme_set(theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
  strip.background = element_rect(fill = "grey95"), legend.key.width = unit(1.6, "lines")))


## Figure 2 (fig_strategies.pdf): missing-data strategies against the missing rate (no outcome heterogeneity)
S = c(oracle = "Oracle (full data)", cc = "Complete-case", native = "Native handling", mi_default = "MI")
N = M[M$row_rate > 0 & M$heterogeneity == "none" & M$strategy %in% names(S), ]
N$group = factor(ifelse(N$outdep, "Missingness depends on the outcome\n(outcome-dependent MAR and NMAR)",
                        "Missingness unrelated to the outcome\n(MCAR, MAR, NMAR)"),
                 levels = c("Missingness unrelated to the outcome\n(MCAR, MAR, NMAR)",
                            "Missingness depends on the outcome\n(outcome-dependent MAR and NMAR)"))
agg = function(y, v) { a = aggregate(y[[v]], list(group = y$group, strategy = y$strategy, rate = y$row_rate), mean)
  a$strategy = factor(S[a$strategy], levels = S); a }
f1 = list(ex = agg(N[N$c == 32, ], "ex"), rel = agg(N, "rel"), t1 = agg(N[N$beta1 == 0, ], "RejectRate"))
panel = function(a, ylab, ref = NULL, reflab = NULL, strip = TRUE, ylim = NULL) {
  g = ggplot(a, aes(100 * rate, x, colour = strategy, shape = strategy, linetype = strategy)) +
    facet_wrap(~ group, nrow = 1) + geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
    scale_colour_manual(values = okabe[c(1, 5, 2, 6)], name = NULL) + scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
    scale_linetype_manual(values = c("dashed", "solid", "solid", "solid"), name = NULL) +
    scale_x_continuous(breaks = 100 * rates) + labs(x = "Row-level missing rate (%)", y = ylab)
  if (!is.null(ref)) g = g + geom_hline(yintercept = ref, colour = "grey45", linetype = "dotted") +
    annotate("text", x = 50, y = ref, label = reflab, hjust = 1, vjust = 1.5, size = 2.6, colour = "grey35")
  if (!is.null(ylim)) g = g + expand_limits(y = ylim)
  if (!strip) g = g + theme(strip.text = element_blank())
  g }
g1 = (panel(f1$ex, "Excess source bias") + theme(axis.title.x = element_blank())) /
     (panel(f1$rel, "Relative MSE", 1, "trial alone", strip = FALSE) + theme(axis.title.x = element_blank())) /
     panel(f1$t1, "Type I error", 0.025, "nominal 0.025", strip = FALSE, ylim = 0.019) +
     plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(out_dir, "fig_strategies.pdf"), g1, width = 6.5, height = 6.2)
cat("wrote fig_strategies.pdf\n")


## Figure 1 (fig_relweight.pdf): relative weight of incomplete subjects against the missing rate, by missingness
## mechanism (c = 32, all scenarios with missingness)
P = c("MCAR", "MAR", "NMAR", "MAR+outdep", "NMAR+outdep")
Plab = c(MCAR = "MCAR", MAR = "MAR", NMAR = "NMAR", "MAR+outdep" = "MAR,\noutcome-dependent", "NMAR+outdep" = "NMAR,\noutcome-dependent")
L = c("logit native" = "Logistic regression with missingness indicators", "bart native" = "BARTm",
      "rf native" = "Random forest, native handling", "xgb native" = "XGBoost, native handling",
      "all mi_default" = "MI (all learners)", "all oracle" = "Oracle (all learners)")
W = M[M$c == 32 & M$row_rate > 0, ]
W$key = ifelse(W$strategy == "native", paste(W$learner, "native"), paste("all", W$strategy))
W = W[W$key %in% names(L), ]
a2 = aggregate(W$relw_inc, list(pattern = W$pattern, key = W$key, rate = W$row_rate), mean)
a2$method = factor(L[a2$key], levels = L); a2$pattern = factor(Plab[a2$pattern], levels = Plab)
g2 = ggplot(a2, aes(100 * rate, x, colour = method, shape = method, linetype = method)) +
  facet_wrap(~ pattern, nrow = 1) + geom_hline(yintercept = 0, colour = "grey45", linetype = "dotted") +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  scale_colour_manual(values = okabe[c(5, 2, 4, 3, 6, 1)], name = NULL) + scale_shape_manual(values = c(17, 15, 3, 4, 18, 16), name = NULL) +
  scale_linetype_manual(values = c(rep("solid", 4), "solid", "dashed"), name = NULL) +
  scale_x_continuous(breaks = 100 * rates) + scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Row-level missing rate (%)", y = "Relative weight of incomplete subjects") +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
ggsave(file.path(out_dir, "fig_relweight.pdf"), g2, width = 6.5, height = 3.2)
cat("wrote fig_relweight.pdf\n")


## Figure S1 (fig_size.pdf, Appendix F): sample-size experiment (results/summary_size.csv)
## MCAR 20% relative weight of incomplete subjects (A) and outcome-dependent MAR 20% excess source bias (B) against the
## number of subjects per propensity model (150 f trial + 200 f source); MI and the oracle / complete-case averaged over learners
z = read.csv("results/summary_size.csv")
z$n = 350 * z$f
Ls = c("logit native" = "Logistic regression with missingness indicators", "bart native" = "BARTm",
       "rf native" = "Random forest, native handling", "xgb native" = "XGBoost, native handling",
       "all mi_default" = "MI (all learners)", "all oracle" = "Oracle (all learners)", "all cc" = "Complete-case (all learners)")
z$key = ifelse(z$strategy == "native", paste(z$learner, "native"), paste("all", z$strategy))
zs = function(set, v, keys) { y = z[z$setting == set & z$key %in% keys, ]
  a = aggregate(y[[v]], list(key = y$key, n = y$n), mean); a$method = factor(Ls[a$key], levels = Ls); a }
cols_s = setNames(okabe[c(5, 2, 4, 3, 6, 1, 1)], Ls); shp_s = setNames(c(17, 15, 3, 4, 18, 16, 1), Ls)
lty_s = setNames(c(rep("solid", 5), "dashed", "dotted"), Ls)
ps = function(a, ylab, title) ggplot(a, aes(n, x, colour = method, shape = method, linetype = method)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  scale_colour_manual(values = cols_s, name = NULL, drop = FALSE) + scale_shape_manual(values = shp_s, name = NULL, drop = FALSE) +
  scale_linetype_manual(values = lty_s, name = NULL, drop = FALSE) +
  scale_x_continuous(trans = "log2", breaks = c(350, 700, 1400, 2800)) +
  labs(x = "Subjects per propensity score model", y = ylab, title = title) + theme(plot.title = element_text(size = 9))
aS = zs("mcar", "relw", names(Ls)); bS = zs("outdep", "ex", names(Ls))   #all 7 lines in both panels, so the legends merge
gS = (ps(aS, "Relative weight of incomplete subjects", "A. MCAR: relative weight") + geom_hline(yintercept = 0, colour = "grey45", linetype = "dotted")) +
     (ps(bS, "Excess source bias", "B. Outcome-dependent MAR: excess bias") + geom_hline(yintercept = 0, colour = "grey45", linetype = "dotted")) +
     plot_layout(guides = "collect") & theme(legend.position = "bottom") & guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2))
ggsave(file.path(out_dir, "fig_size.pdf"), gS, width = 6.5, height = 4)
cat("wrote fig_size.pdf\n")
cat("\n--- sample-size run ---\n"); print(xtabs(round(x, 3) ~ method + n, aS)); print(xtabs(round(x, 3) ~ method + n, bS))


## numbers for the text
cat("\n--- numbers for the text ---\n")
w = function(p, k) round(a2$x[a2$pattern == Plab[p] & a2$key == k & a2$rate %in% c(0.1, 0.5)], 2)
for (k in names(L)) cat(sprintf("%-16s MCAR 10%%/50%%: %s\n", k, paste(w("MCAR", k), collapse = " -> ")))
z = a2[a2$rate == 0.5 & a2$key %in% c("bart native", "rf native", "xgb native"), ]
cat("native ensembles at 50% by pattern:", paste(levels(z$pattern), round(tapply(z$x, z$pattern, mean), 2)), "\n")
for (v in names(f1)) { cat("\n", v, "\n"); print(round(xtabs(x ~ strategy + rate + group, f1[[v]]), 3)) }

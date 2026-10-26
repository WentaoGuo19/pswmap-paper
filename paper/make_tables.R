# LaTeX tables 3-5 of the Results section of Guo and Zhu, computed from the evaluated simulation results.
# Run from the repository root:  Rscript paper/make_tables.R
# Reads results/summary_main.csv and results/summary_pooled.csv (written by evaluation/evaluation.R) and writes one
# .tex file per table to paper/tables/. It also prints the numbers quoted in the text.
#
# Definitions (main grid unless stated):
#   relative MSE       = MSE of a method / MSE of the RCT-only analysis in the same scenario and c
#   excess source bias = bias of a source's weighted mean minus that of the oracle (same learner, scenario), averaged
#                        over the three sources: the bias added by the missing-data strategy
#   MI                 = MI with the default engine (the CART and random-forest engines give nearly the same results)
# Averages run over learners, shifts, selection models, missing rates and c, as stated in each caption.

out_dir = "paper/tables"
dir.create(out_dir, showWarnings = FALSE)
f = function(x, d = 3) formatC(x, d, format = "f")
fs = function(x, d = 3) ifelse(abs(x) < 0.5 * 10^-d, f(0, d), ifelse(x < 0, paste0("$-$", f(abs(x), d)), f(x, d)))   #signed: no "-0.000", true minus sign
write_tab = function(lines, file) { writeLines(lines, file.path(out_dir, file)); cat("wrote", file, "\n") }

d0 = read.csv("results/summary_main.csv")
rct = d0[d0$method == "rct_only", c("sid", "c", "MSE")]; names(rct)[3] = "MSE_rct"
D = merge(d0[d0$method != "rct_only", ], rct, by = c("sid", "c")); D$rel = D$MSE / D$MSE_rct
D$learner = sub("_.*", "", D$method); D$strategy = sub("^[a-z]+_", "", D$method)
D$sb = rowMeans(D[, paste0("bias_src_", 1:3)])
o = D[D$strategy == "oracle", c("sid", "c", "learner", "sb")]; names(o)[4] = "sb_or"
D = merge(D, o, by = c("sid", "c", "learner")); D$ex = D$sb - D$sb_or
M = D[D$grid == "main", ]
P = c("MCAR", "MAR", "NMAR", "MAR+outdep", "NMAR+outdep")
Plab = c(MCAR = "MCAR", MAR = "MAR", NMAR = "NMAR", "MAR+outdep" = "MAR, outcome-dependent", "NMAR+outdep" = "NMAR, outcome-dependent")
S = c(oracle = "Oracle (full data)", cc = "Complete-case", native = "Native handling", mi_default = "MI")


## (the relative weight of incomplete subjects is shown in Figure 1, see make_figures.R)


## Table 3 (tab_main.tex): missing-data strategies by missingness mechanism (no outcome heterogeneity)
N = M[M$row_rate > 0 & M$heterogeneity == "none", ]
cell = function(p, s) { y = N[N$pattern == p & N$strategy == s, ]; yd = y[y$c == 32, ]
  c(ex = mean(yd$ex), rel = mean(y$rel), t1 = mean(y$RejectRate[y$beta1 == 0]), ess = mean(y$avg_ESS)) }
body = unlist(lapply(P, function(p) c(
  sprintf("\\multicolumn{5}{l}{\\textit{%s}} \\\\", Plab[p]),
  sapply(names(S), function(s) { v = cell(p, s)
    sprintf("\\quad %s & %s & %s & %s & %s \\\\", S[s], if (s == "oracle") "(reference)" else fs(v["ex"]), f(v["rel"]), f(v["t1"]), f(v["ess"], 1)) }))))
r0 = d0[d0$method == "rct_only" & d0$grid == "main" & d0$row_rate > 0 & d0$heterogeneity == "none" & d0$beta1 == 0, ]
write_tab(c("\\begin{table}[ht]", "\\centering",
  sprintf("\\caption{Missing-data strategies by missingness mechanism, without outcome heterogeneity. Excess source bias is the bias of a source's weighted mean beyond that of the oracle with the same learner, as a signed mean over the three sources. Prior ESS is the prior effective sample size of the robust PSW-MAP prior, and relative MSE is the MSE of the treatment effect divided by that of the RCT-only analysis (type I error of RCT-only: %s). Averages over learners, covariate shift, selection model, missing rates 10--50\\%%, $c$ and 1{,}000 replicates per scenario. MI uses the default imputation engine.}", f(mean(r0$RejectRate))),
  "\\label{tab:main}", "\\begin{tabular}{lcccc}", "\\toprule",
  "Strategy & Excess source bias & Relative MSE & Type I error & Prior ESS \\\\", "\\midrule", body,
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), "tab_main.tex")


## Table 4 (tab_dose.tex): strength of outcome dependence (sensitivity grid; the missing-rate dose-response is Figure 2)
Sg = D[D$outdep & D$shift == "low" & D$selection == "linear" & D$c == 32 & D$row_rate %in% c(0.2, 0.4) & D$heterogeneity == "none", ]
Sg$or = ifelse(is.na(Sg$gamma_y), 2, round(exp(Sg$gamma_y), 2)); ors = c(1.25, 1.5, 2, 3)
pb = function(s, v) sapply(ors, function(g) { y = Sg[Sg$strategy == s & Sg$or == g, ]
  switch(v, ex = mean(y$ex), rel = mean(y$rel), t1 = mean(y$RejectRate[y$beta1 == 0])) })
S3 = S[c("cc", "native", "mi_default")]; S4 = S[c("oracle", "cc", "native", "mi_default")]
blk = function(fun, v, d, SS = S3) sapply(names(SS), function(s) sprintf("\\quad %s & %s \\\\", SS[s], paste(fs(fun(s, v), d), collapse = " & ")))
write_tab(c("\\begin{table}[ht]", "\\centering",
  "\\caption{Strength of outcome dependence: excess source bias, relative MSE and type I error by the odds ratio of missingness per standard deviation of the outcome, for outcome-dependent MAR and NMAR with low covariate shift, linear selection, missing rates 20\\% and 40\\% and $c = 32$. No outcome heterogeneity. Averages over learners and 1{,}000 replicates per scenario. The dependence on the missing rate is shown in Figure~\\ref{fig:strategies}.}",
  "\\label{tab:dose}", "\\begin{tabular}{lcccc}", "\\toprule",
  " & \\multicolumn{4}{c}{Odds ratio of missingness} \\\\", "\\cmidrule(lr){2-5}",
  " & 1.25 & 1.5 & 2 & 3 \\\\", "\\midrule",
  "\\multicolumn{5}{l}{\\textit{Excess source bias}} \\\\", blk(pb, "ex", 3),
  "\\multicolumn{5}{l}{\\textit{Relative MSE}} \\\\", blk(pb, "rel", 3, S4),
  "\\multicolumn{5}{l}{\\textit{Type I error}} \\\\", blk(pb, "t1", 3, S4),
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), "tab_dose.tex")


## Table 5 (tab_pooled.tex): source-specific versus pooled propensity score (16 MAR 20% scenarios)
p0 = read.csv("results/summary_pooled.csv"); sids = unique(p0$sid)
mk = function(d, scope) { d = d[d$method != "rct_only", ]; d$scope = scope; d$learner = sub("_.*", "", d$method)
  d$strategy = sub("^[a-z]+_(pooled_)?", "", d$method); d[d$strategy %in% c("oracle", "cc", "mi_default"), ] }
Q = rbind(mk(d0[d0$sid %in% sids, ], "source"), mk(p0, "pooled")[, names(mk(d0[d0$sid %in% sids, ], "source"))])
Q = merge(Q, rct, by = c("sid", "c")); Q$rel = Q$MSE / Q$MSE_rct
qrow = function(sh, sc) { y = Q[Q$shift == sh & Q$scope == sc, ]; yd = y[y$c == 32, ]; yn = y[y$heterogeneity == "none", ]
  c(colMeans(yd[, paste0("bias_src_", 1:3)]), mean(yd$mae_src), mean(yd$kish), mean(yn$avg_ESS), mean(yn$rel)) }
lab = c(source = "Source-specific", pooled = "Pooled")
tb = unlist(lapply(c("low", "high"), function(sh) c(sprintf("\\multicolumn{8}{l}{\\textit{%s covariate shift}} \\\\", if (sh == "low") "Low" else "High"),
  sapply(c("source", "pooled"), function(sc) { v = qrow(sh, sc)
    sprintf("\\quad %s & %s & %s & %s & %s & %s & %s & %s \\\\", lab[sc], fs(v[1]), fs(v[2]), fs(v[3]), f(v[4]), f(v[5], 1), f(v[6], 1), f(v[7])) }))))
write_tab(c("\\begin{table}[ht]", "\\centering",
  "\\caption{Source-specific versus pooled propensity scores, applied to the same simulated datasets (16 scenarios with MAR missingness at a 20\\% row-level rate). Source bias is the bias of each source's weighted mean, and Kish ESS is the effective sample size of a source's weights. Prior ESS and relative MSE (MSE of the treatment effect divided by that of the RCT-only analysis) are without outcome heterogeneity. Averages over learners, the oracle, complete-case and MI strategies, $c$ and 1{,}000 replicates per scenario.}",
  "\\label{tab:pooled}", "\\begin{tabular}{lccccccc}", "\\toprule",
  " & \\multicolumn{3}{c}{Source bias} & Mean & Kish & Prior & Relative \\\\", "\\cmidrule(lr){2-4}",
  "Propensity score & $H_1$ & $H_2$ & $H_3$ & $|$bias$|$ & ESS & ESS & MSE \\\\", "\\midrule", tb,
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), "tab_pooled.tex")


## numbers quoted in the text (printed, not written to a file)
cat("\n--- numbers for the text ---\n")
y = M[M$row_rate > 0 & M$strategy == "native" & M$learner != "logit" & M$c == 32, ]
cat("native relw (BART/RF/XGB) at 50% by pattern:", paste(P, f(tapply(y$relw_inc[y$row_rate == 0.5], y$pattern[y$row_rate == 0.5], mean)[P], 2)), "\n")
y = M[M$row_rate > 0 & !M$outdep & M$heterogeneity == "none", ]
cat("M indep Y|X, relative MSE MI vs native vs cc (all):", f(mean(y$rel[y$strategy == "mi_default"])), f(mean(y$rel[y$strategy == "native"])), f(mean(y$rel[y$strategy == "cc"])), "\n")
cat("   prior ESS MI vs native vs cc vs oracle:", f(mean(y$avg_ESS[y$strategy == "mi_default"]), 1), f(mean(y$avg_ESS[y$strategy == "native"]), 1), f(mean(y$avg_ESS[y$strategy == "cc"]), 1), f(mean(y$avg_ESS[y$strategy == "oracle"]), 1), "\n")
y = M[M$row_rate > 0 & M$outdep & M$c == 32, ]
cat("outdep: share of cells with all 3 sources biased down (rate >= 30%):", paste(names(S), f(tapply(rowSums(y[y$row_rate >= 0.3, paste0("bias_src_", 1:3)] < 0) == 3, y$strategy[y$row_rate >= 0.3], mean)[names(S)], 2)), "\n")
y = M[M$row_rate > 0, ]
spread = function(z, by) mean(tapply(z$rel, by, function(v) diff(range(v))))
for (pp in list(c("MCAR", "MAR"), "NMAR", c("MAR+outdep", "NMAR+outdep"))) { z = y[y$pattern %in% pp, ]
  cat("spread of relative MSE,", paste(pp, collapse = "/"), ": across strategies", f(spread(z, paste(z$sid, z$c, z$learner))), " across learners", f(spread(z, paste(z$sid, z$c, z$strategy))), "\n") }
B = M[M$strategy == "mi_default" | (M$row_rate == 0 & M$strategy == "oracle"), ]
cat("MI vs RCT-only: relative MSE by shift x heterogeneity:\n"); print(round(tapply(B$rel, list(B$shift, B$heterogeneity), mean), 3))
cat("share of cells MI beats RCT-only, by heterogeneity x shift x c:\n"); print(round(tapply(B$rel < 1, list(paste(B$heterogeneity, B$shift), B$c), mean), 3))
cat("MI type I by heterogeneity x shift:\n"); print(round(tapply(B$RejectRate[B$beta1 == 0], list(B$heterogeneity[B$beta1 == 0], B$shift[B$beta1 == 0]), mean), 3))
cat("MI power by heterogeneity x shift, and RCT-only power:\n"); print(round(tapply(B$RejectRate[B$beta1 == 0.7], list(B$heterogeneity[B$beta1 == 0.7], B$shift[B$beta1 == 0.7]), mean), 3))
r1 = d0[d0$method == "rct_only" & d0$grid == "main" & d0$beta1 == 0.7, ]; cat("RCT-only power:", f(mean(r1$RejectRate)), " type I:", f(mean(d0$RejectRate[d0$method == "rct_only" & d0$grid == "main" & d0$beta1 == 0])), "\n")
L = M[M$row_rate == 0 & M$strategy == "oracle", ]; L$setting = paste(L$shift, L$selection)
cat("learner relative MSE at 0% missing:\n"); print(round(tapply(L$rel, list(L$learner, L$setting), mean), 3))
cat("pooled vs source by learner, mean |source bias| (c32):\n"); qq = Q[Q$c == 32, ]; print(round(tapply(qq$mae_src, list(qq$learner, qq$scope), mean), 3))
qa = merge(Q[Q$scope == "source", ], Q[Q$scope == "pooled", ], by = c("sid", "c", "learner", "strategy"), suffixes = c(".s", ".p"))
qa$worse = qa$MSE.p - qa$MSE.s > 2 * sqrt(qa$MCSE_MSE.s^2 + qa$MCSE_MSE.p^2); qa$better = qa$MSE.s - qa$MSE.p > 2 * sqrt(qa$MCSE_MSE.s^2 + qa$MCSE_MSE.p^2)
cat("pooled significantly worse / better share by shift x heterogeneity:\n"); print(round(cbind(worse = tapply(qa$worse, paste(qa$shift.s, qa$heterogeneity.s), mean), better = tapply(qa$better, paste(qa$shift.s, qa$heterogeneity.s), mean)), 3))
cat("pooled prior ESS by shift x heterogeneity (source / pooled):\n"); print(round(tapply(Q$avg_ESS, list(paste(Q$shift, Q$heterogeneity), Q$scope), mean), 1))
cat("MCSE: median per-cell MCSE of type I", f(median(d0$MCSE_RejectRate[d0$beta1 == 0])), "; max", f(max(d0$MCSE_RejectRate[d0$beta1 == 0])), "\n")

# Runs a few replicates of one scenario of the simulation study on a laptop, without a cluster.
# Run from the repository root:  Rscript examples/run_one_scenario.R
#
# The example uses scenario 18 of the main grid (low covariate shift, linear selection, no outcome heterogeneity,
# outcome-dependent MAR at a 20% row-level rate, treatment effect 0) and, to keep it fast, three learners without
# BART. One replicate takes about a minute.

library(pswmap)
for (f in c("simulate_data.R", "scenarios.R", "diagnostics.R", "map_guard.R", "run_scenario.R")) source(file.path("R", f))

s = make_scenarios()[18, ]
print(s[, c("sid", "shift", "selection", "heterogeneity", "pattern", "row_rate", "beta1")])

res = run_scenario(n_treat = 100, n_control = 50, n_hist = rep(200, 3), mu_current = 1,
                   mu_hist = c(s$mu_h1, s$mu_h2, s$mu_h3), delta = c(s$d1, s$d2, s$d3), Vtimer = c(s$v1, s$v2, s$v3),
                   c = c(16, 4), trial_bounds = parse_bounds(s$trial_bounds), beta1 = s$beta1,
                   beta0_com = 0.1, beta0_mis = 0.3, mechanism = s$mechanism, row_rate = s$row_rate,
                   outdep = s$outdep, gamma_y = s$gamma_y, learners = c("logit", "rf", "xgb"),
                   mi_engines = "default", PS_cap = 0.95, reps = 1:2, sid = s$sid)

# Posterior mean of the treatment effect (true value 0), prior effective sample size, and the average relative
# weight of incomplete external subjects, for each method at c = 16, averaged over the replicates
tab = t(sapply(res$c16, function(m) c(effect = mean(m[, "mean"]), prior_ESS = mean(m[, "ESS"]),
                                      relweight_incomplete = mean(m[, paste0("relw_inc_", 1:3)]))))
print(round(tab, 3))
# Complete-case weighting and native handling give incomplete subjects little or no weight, and here (missingness
# depends on the outcome) they shift the borrowed information. MI keeps the incomplete subjects, like the oracle.
# With 1,000 replicates per scenario, evaluation/evaluation.R turns such results into the performance measures of
# the paper.

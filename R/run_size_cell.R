#' @title Run one cell of the sample-size experiment
#'
#' @description The sample-size experiment (Appendix F of Guo and Zhu) fixes the missing rate at 20% and scales the
#'    trial (100f treated, 50f control) and every external source (200f) by the factor f, so each propensity score
#'    model has 350f subjects. It uses low covariate shift, linear selection and no outcome heterogeneity. For each
#'    replicate, the propensity score weights of the four learners under the oracle, complete-case, native and MI
#'    (default engine) strategies are computed with pswmap::psw_weights(), and their weighting diagnostics are
#'    recorded. The MAP analysis is not run, because the weights and weighted source means are the quantities of
#'    interest.
#'
#'    Seeds: replicate r of cell q uses set.seed(10000000 + 10000 * q + r), which never coincides with a seed of the
#'    main simulation (10000 * sid + r with sid <= 464).
#'
#' @param f size factor (1, 2, 4 or 8 in the study)
#' @param mechanism,outdep missingness mechanism: "MCAR" with outdep = FALSE, or "MAR" with outdep = TRUE
#' @param cell cell number q, used in the seed
#' @param reps replicate indices
#' @param bart_args arguments for bartMachine (default: its own defaults)
#'
#' @return a data frame with one row per replicate and method: cell, f, setting ("mcar" or "outdep"), rep, method and the
#'    diagnostics of weight_diagnostics()
run_size_cell = function(f, mechanism, outdep, cell, reps, bart_args = list())
{
  xv = paste0("X", 1:10)
  learners = c("logit", "bart", "rf", "xgb")
  out = vector("list", length(reps))
  for (i in seq_along(reps)) {
    r = reps[i]
    set.seed(10000000 + 10000 * cell + r)
    dat_com = simdat_complete(n_treat = 100 * f, n_control = 50 * f, n_hist = rep(200 * f, 3), mu_current = 1,
                              mu_hist = c(1, 0.5, 1.2), delta = rep(0, 3), Vtimer = rep(1, 3), beta1 = 0,
                              beta0_com = 0.1, beta0_mis = 0.3)
    dat_mis = simdat_missing(dat_com, mechanism = mechanism, row_rate = 0.2, outdep = outdep)
    trial = data.frame(Y = dat_mis$RCT$Y, trt = dat_mis$RCT$treatment, dat_mis$RCT[, xv])
    ext_com = lapply(dat_com$RWD, function(d) data.frame(Y = d$Y, d[, xv]))
    ext_mis = lapply(dat_mis$RWD, function(d) data.frame(Y = d$Y, d[, xv]))
    incomplete = lapply(ext_mis, function(d) !complete.cases(d[, xv]))

    w = list()
    for (l in learners) {
      la = if (l == "bart") bart_args else list()
      w[[paste0(l, "_oracle")]] = pswmap::psw_weights(trial, ext_com, learner = l, learner_args = la)
      w[[paste0(l, "_cc")]] = pswmap::psw_weights(trial, ext_mis, learner = l, learner_args = la, missing = "cc")
      w[[paste0(l, "_native")]] = pswmap::psw_weights(trial, ext_mis, learner = l, learner_args = la, missing = "native")
      w[[paste0(l, "_mi_default")]] = pswmap::psw_weights(trial, ext_mis, learner = l, learner_args = la, missing = "mi")
      if (l == "bart") rJava::.jgc(R.gc = TRUE)
    }
    d = t(sapply(w, weight_diagnostics, incomplete = incomplete, rwd_full = dat_com$RWD, rct_full = dat_com$RCT))
    out[[i]] = data.frame(cell = cell, f = f, setting = if (outdep) "outdep" else "mcar", rep = r, method = names(w), d,
                          row.names = NULL)
    if (i %% 10 == 0) cat(i, "of", length(reps), "replicates done, current index", r, "\n")
  }
  do.call(rbind, out)
}

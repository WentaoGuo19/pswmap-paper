#' @title Run one simulation scenario
#'
#' @description Runs the replicates `reps` of one scenario of the simulation study (Section 3 of Guo and Zhu). In
#'    each replicate it:
#'    1. generates a trial and K external sources (simdat_complete()), and makes covariates missing in the external
#'       sources (simdat_missing());
#'    2. computes the source-specific propensity score weights of every method with pswmap::psw_weights(). The
#'       methods cross the learners with the strategies oracle (complete data before missingness), complete-case,
#'       native handling and multiple imputation with each engine, 24 methods by default;
#'    3. for each method and each value of the heterogeneity constant c, builds the robust PSW-MAP prior
#'       (pswmap::psw_map_prior(), with a time limit, see map_prior_guarded()) and the posterior of the treatment
#'       effect (pswmap::psw_posterior());
#'    4. analyses the trial alone (pswmap::pswmap(borrow = FALSE)).
#'
#'    The propensity scores and weights do not depend on c, so they are computed once and analysed under every c
#'    (paired comparison). With no missing covariates every strategy of a learner equals its oracle, so those
#'    results are copied. Performance measures are computed later by evaluation/evaluation.R, after the chunks of a
#'    scenario are combined.
#'
#'    Seeds: replicate r of scenario sid uses set.seed(10000 * sid + r) right before data generation. Every later
#'    random step (imputation, learners, MCMC, posterior draws) continues the same stream, except bartMachine, which
#'    uses its own Java random number generator.
#'
#' @param n_treat,n_control,n_hist,mu_current,mu_hist,delta,Vtimer,trial_bounds,beta1,beta0_com,beta0_mis data
#'    generation, see simdat_complete()
#' @param mechanism,row_rate,outdep,gamma_y missingness, see simdat_missing()
#' @param c vector of heterogeneity constants (e.g. c(32, 16, 8, 4)); the prior on tau is half-normal with scale
#'    sigma_ref / c
#' @param learners propensity score learners, subset of c("logit", "bart", "rf", "xgb")
#' @param mi_engines imputation engines, subset of c("default", "cart", "rf")
#' @param mi_m number of imputed datasets
#' @param native include the native missing-data handling (FALSE in the pooled-score comparison)
#' @param ps_scope "source" (one propensity score per source, the PSW-MAP method) or "pooled" (one score for the
#'    trial against all sources; methods are then named <learner>_pooled_<strategy>)
#' @param PS_cap cap on the propensity score
#' @param bart_args arguments for bartMachine (default: its own defaults)
#' @param map MCMC and mixture settings, see pswmap::map_control()
#' @param timeout,max_tries time limit and attempts of each MAP fit, see map_prior_guarded()
#' @param reps replicate indices to run (whole numbers below 10,000)
#' @param sid scenario id (row of make_scenarios()), used in the seed; NULL runs unseeded replicates
#'
#' @return a list with one element per value of c (named c32, c16, ...). Each is a list with one matrix per method
#'    (<learner>_<strategy>, plus rct_only), with one row per replicate and columns: mean, sd, 2.5%, 50%, 97.5% (the
#'    posterior of the treatment effect), ESS (prior effective sample size; 0 for rct_only) and the per-source
#'    diagnostics of weight_diagnostics() (NA for rct_only). If a MAP fit fails in every attempt, its estimates are
#'    NA and its diagnostics are kept.
run_scenario = function(n_treat, n_control, n_hist, mu_current, mu_hist, delta, Vtimer, c, trial_bounds = NULL,
                        beta1, beta0_com, beta0_mis, mechanism, row_rate, outdep = FALSE, gamma_y = log(2),
                        learners = c("logit", "bart", "rf", "xgb"), mi_engines = c("default", "cart", "rf"),
                        mi_m = 5, native = TRUE, ps_scope = "source", PS_cap = 0.95, bart_args = list(),
                        map = pswmap::map_control(cores = 2), timeout = 60, max_tries = 3, reps, sid = NULL)
{
  if (!is.numeric(reps) || length(reps) == 0 || anyNA(reps) || any(reps < 1) || any(reps != round(reps)) ||
      anyDuplicated(reps)) stop("reps must be distinct positive whole numbers")
  if (!is.null(sid) && max(reps) >= 10000) stop("reps must be below 10000, so that seeds stay unique across scenarios")
  K = length(n_hist)
  for (nm in c("mu_hist", "delta", "Vtimer"))
    if (length(get(nm)) != K) stop(nm, " must have length K = length(n_hist) = ", K)
  if (!all(learners %in% c("logit", "bart", "rf", "xgb"))) stop("unknown learner")
  if (!all(mi_engines %in% c("default", "cart", "rf"))) stop("unknown imputation engine")

  xv = paste0("X", 1:10)
  strategies = c("oracle", "cc", if (native) "native", if (length(mi_engines)) paste0("mi_", mi_engines))
  tag = if (ps_scope == "pooled") "_pooled_" else "_"
  methods = c(as.vector(t(outer(learners, strategies, paste, sep = tag))), "rct_only")
  diag_cols = paste0(rep(c("ybar_", "se_", "kish_", "relw_inc_", "smd_mean_", "smd_max_"), each = K), seq_len(K))
  cols = c("mean", "sd", "2.5%", "50%", "97.5%", "ESS", diag_cols)
  c_names = paste0("c", c)
  results = setNames(lapply(c_names, function(cn)
    setNames(lapply(methods, function(m) matrix(NA_real_, length(reps), length(cols), dimnames = list(reps, cols))),
             methods)), c_names)

  #data in the format of the pswmap package: a trial data frame (Y, trt, X1-X10) and a list of external sources
  as_trial = function(d) data.frame(Y = d$Y, trt = d$treatment, d[, xv])
  as_external = function(rwd) lapply(rwd, function(d) data.frame(Y = d$Y, d[, xv]))
  weights = function(ext, ...) {
    args = list(...)
    if (identical(args$learner, "bart")) args$learner_args = bart_args
    w = do.call(pswmap::psw_weights, c(list(trial = trial, external = ext, ps_cap = PS_cap, ps_scope = ps_scope), args))
    if (identical(args$learner, "bart")) rJava::.jgc(R.gc = TRUE)   #free Java memory of finished BART fits
    w
  }

  for (i in seq_along(reps)) {
    r = reps[i]
    #1. data
    if (!is.null(sid)) set.seed(sid * 10000 + r)
    dat_com = simdat_complete(n_treat = n_treat, n_control = n_control, n_hist = n_hist, mu_current = mu_current,
                              mu_hist = mu_hist, delta = delta, Vtimer = Vtimer, beta1 = beta1, beta0_com = beta0_com,
                              beta0_mis = beta0_mis, trial_bounds = trial_bounds)
    dat_mis = simdat_missing(dat_com, mechanism = mechanism, row_rate = row_rate, outdep = outdep, gamma_y = gamma_y)
    trial = as_trial(dat_mis$RCT)
    ext_com = as_external(dat_com$RWD)
    ext_mis = as_external(dat_mis$RWD)
    incomplete = lapply(ext_mis, function(d) !complete.cases(d[, xv]))
    no_missing = !any(unlist(incomplete))

    #2. propensity score weights of every method
    #(warnings of the pooled score, which is deliberately part of the comparison, are not shown)
    w = list()
    suppressWarnings({
      for (l in learners) {
        w[[paste0(l, tag, "oracle")]] = weights(ext_com, learner = l)
        if (!no_missing) {
          w[[paste0(l, tag, "cc")]] = weights(ext_mis, learner = l, missing = "cc")
          if (native) w[[paste0(l, tag, "native")]] = weights(ext_mis, learner = l, missing = "native")
        }
      }
      if (!no_missing) for (e in mi_engines) for (l in learners)
        w[[paste0(l, tag, "mi_", e)]] = weights(ext_mis, learner = l, missing = "mi", mi_engine = e, mi_m = mi_m)
    })

    #3. PSW-MAP prior and posterior, for every method and value of c
    sigma_ref = pswmap::psw_reference_scale(trial, ext_mis)
    y0 = trial$Y[trial$trt == 0]; y1 = trial$Y[trial$trt == 1]
    for (m in names(w)) {
      diag_m = weight_diagnostics(w[[m]], incomplete, dat_com$RWD, dat_com$RCT)
      summaries_m = pswmap::psw_summaries(w[[m]])
      for (k in seq_along(c)) {
        prior = map_prior_guarded(summaries_m, sigma_ref, c[k], control = map, timeout = timeout,
                                  max_tries = max_tries)
        if (is.null(prior)) {
          message("run_scenario: no MAP prior for method ", m, ", c = ", c[k], ", replicate ", r, " -> NA")
          results[[k]][[m]][i, names(diag_m)] = diag_m
          next
        }
        post = pswmap::psw_posterior(prior$robust, y0, y1, sigma_ref)
        e = post$effect
        results[[k]][[m]][i, ] = c(e[["mean"]], e[["sd"]], e[["lower"]], e[["median"]], e[["upper"]], prior$ess, diag_m)
      }
    }
    #no missing covariates: every strategy of a learner equals its oracle
    if (no_missing) for (cn in c_names) for (l in learners) for (st in strategies[-1])
      results[[cn]][[paste0(l, tag, st)]][i, ] = results[[cn]][[paste0(l, tag, "oracle")]][i, ]

    #4. the trial alone (does not depend on c)
    e = pswmap::pswmap(trial, ext_mis, borrow = FALSE)$effect
    for (cn in c_names) results[[cn]]$rct_only[i, c("mean", "sd", "2.5%", "50%", "97.5%", "ESS")] =
      c(e[["mean"]], e[["sd"]], e[["lower"]], e[["median"]], e[["upper"]], 0)

    if (i %% 10 == 0) cat(i, "of", length(reps), "replicates done, current index", r, "\n")
  }
  attr(results, "sid") = sid
  attr(results, "reps") = reps
  attr(results, "c") = c
  results
}

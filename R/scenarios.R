#' @title Simulation scenario table
#'
#' @description Builds the scenario grid of the simulation study (Section 3.4 and Appendix F of Guo and Zhu). The
#'    row number is the scenario id (sid), which also determines the seeds. Used by the cluster scripts (to run a
#'    scenario) and by evaluation/evaluation.R (to label the results).
#'
#'    Main grid, 416 data-generating scenarios:
#'      covariate shift (low, high) x selection model (linear, nonlinear) x outcome heterogeneity (none, medium)
#'      x missingness condition (26) x treatment effect (0, 0.7).
#'    The 26 missingness conditions are no missingness plus 5 mechanisms (MCAR, MAR, NMAR, outcome-dependent MAR and
#'    NMAR, the latter with gamma_y = log 2) at row-level rates of 10% to 50%.
#'    Selection model: "linear" has no trial eligibility criteria, so the true log-odds of trial membership is close
#'    to linear in the covariates. "nonlinear" applies the trial eligibility criteria X5 >= 0, X6 <= 2 and X8 >= 0.5
#'    to the trial only, so the true propensity score changes abruptly at the bounds. The cut-offs on the prognostic
#'    covariates X6 and X8, which may be missing in the external data, make this misspecification of a main-effects
#'    logistic propensity score matter for the outcome.
#'    The heterogeneity constant c is not a scenario factor: it only enters the MAP step, so every simulated dataset
#'    is analysed under each value in c_values (main grid: 32, 16, 8, 4).
#'
#'    Sensitivity grid, 48 scenarios appended after the main grid: the outcome-dependent mechanisms with
#'    gamma_y in {log 1.25, log 1.5, log 3}, row-level rates 20% and 40%, both treatment effects and both outcome
#'    heterogeneity levels, with low covariate shift, linear selection and c = 32.
#'
#' @return data.frame with one row per scenario:
#'    sid, grid ("main"/"sensitivity"), labels (shift, selection, heterogeneity, pattern),
#'    mechanism, outdep, gamma_y (NA if not outcome-dependent), row_rate, beta1,
#'    c_values (tau prior constants analysed in this scenario, comma-separated, e.g. "32,16,8,4"),
#'    mu_h1..mu_h3 (covariate means), v1..v3 (Vtimer), d1..d3 (delta),
#'    trial_bounds (eligibility criteria as text, "" for none; parse with parse_bounds()),
#'    target_mean (E_trial[beta0'X], the covariate-adjusted target of each source's weighted mean when delta_k = 0)
make_scenarios = function()
{
  #factor levels (mu_current = 1 for the trial)
  shift = list(low = c(1, 0.5, 1.2), high = c(0, 0.25, 1.5))             #mu_hist
  selection = list(linear = c(1, 1, 1), nonlinear = c(1, 1, 1))         #Vtimer (same covariance in all studies)
  bounds = c(linear = "", nonlinear = "X5>=0;X6<=2;X8>=0.5")             #trial eligibility criteria
  #E_trial[beta0'X] (mu_current = 1, beta0 = 0.1/0.3): 0.8 * pnorm(1) + 1.2 without bounds. With X5 >= 0, X6 <= 2 and
  #X8 >= 0.5, 2.0175 by Monte Carlo (2e7 draws, Monte Carlo SE 0.0002; 48.6% of trial-like subjects are eligible)
  target = c(linear = 0.8 * pnorm(1) + 1.2, nonlinear = 2.0175)
  heterogeneity = list(none = c(0, 0, 0), medium = c(0, -0.3, -0.5))     #delta
  patterns = data.frame(pattern = c("MCAR", "MAR", "NMAR", "MAR+outdep", "NMAR+outdep"),
                        mechanism = c("MCAR", "MAR", "NMAR", "MAR", "NMAR"),
                        outdep = c(FALSE, FALSE, FALSE, TRUE, TRUE))

  #helper: expand factor labels into the numeric columns run_scenario() needs
  add_values = function(g) {
    g$mechanism = patterns$mechanism[match(g$pattern, patterns$pattern)]
    g$mechanism[g$pattern == "none"] = "MCAR"   #row_rate = 0, mechanism unused
    g$outdep = g$pattern %in% c("MAR+outdep", "NMAR+outdep")
    mu = t(sapply(g$shift, function(x) shift[[x]]))
    v = t(sapply(g$selection, function(x) selection[[x]]))
    d = t(sapply(g$heterogeneity, function(x) heterogeneity[[x]]))
    cbind(g, mu_h1 = mu[,1], mu_h2 = mu[,2], mu_h3 = mu[,3],
          v1 = v[,1], v2 = v[,2], v3 = v[,3], d1 = d[,1], d2 = d[,2], d3 = d[,3],
          trial_bounds = unname(bounds[g$selection]), target_mean = unname(target[g$selection]))
  }

  #1.main grid
  miss = rbind(data.frame(pattern = "none", row_rate = 0),
               expand.grid(row_rate = c(0.1, 0.2, 0.3, 0.4, 0.5), pattern = patterns$pattern,
                           stringsAsFactors = FALSE)[, c("pattern", "row_rate")])
  main = expand.grid(miss_id = seq_len(nrow(miss)), heterogeneity = names(heterogeneity),
                     selection = names(selection), shift = names(shift),
                     beta1 = c(0, 0.7), stringsAsFactors = FALSE)
  main = cbind(main[, c("shift", "selection", "heterogeneity")], miss[main$miss_id, ], beta1 = main$beta1)
  main$c_values = "32,16,8,4"
  main = add_values(main)
  main$gamma_y = ifelse(main$outdep, log(2), NA)
  main$grid = "main"

  #2.sensitivity grid for gamma_y
  sens = expand.grid(gamma_y = log(c(1.25, 1.5, 3)), row_rate = c(0.2, 0.4),
                     pattern = c("MAR+outdep", "NMAR+outdep"), heterogeneity = names(heterogeneity),
                     beta1 = c(0, 0.7), stringsAsFactors = FALSE)
  sens$shift = "low"; sens$selection = "linear"; sens$c_values = "32"
  sens = add_values(sens)
  sens$grid = "sensitivity"

  #combine
  cols = c("grid", "shift", "selection", "heterogeneity", "pattern", "mechanism", "outdep", "gamma_y",
           "row_rate", "beta1", "c_values", "mu_h1", "mu_h2", "mu_h3", "v1", "v2", "v3", "d1", "d2", "d3",
           "trial_bounds", "target_mean")
  out = rbind(main[, cols], sens[, cols])
  out = cbind(sid = seq_len(nrow(out)), out)
  rownames(out) = NULL
  return(out)
}



#' @title Parse the trial eligibility criteria of a scenario
#'
#' @param txt text such as "X5>=0;X6<=2" (criteria separated by ";", each "Xj>=a" or "Xj<=b"), or "" for none
#'
#' @return NULL if txt is empty, otherwise a named list of c(lower, upper) bounds for simdat_complete(trial_bounds = )
parse_bounds = function(txt)
{
  if (is.na(txt) || txt == "") return(NULL)
  out = list()
  for (cr in strsplit(txt, ";")[[1]]) {
    v = sub("^(X[0-9]+).*", "\\1", cr); val = as.numeric(sub("^X[0-9]+(>=|<=)", "", cr))
    if (is.null(out[[v]])) out[[v]] = c(-Inf, Inf)
    if (grepl(">=", cr, fixed = TRUE)) out[[v]][1] = val else if (grepl("<=", cr, fixed = TRUE)) out[[v]][2] = val else stop("bad criterion: ", cr)
  }
  out
}

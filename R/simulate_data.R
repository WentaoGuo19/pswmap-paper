library(MASS)



#' @title Simulate a complete trial and K external control sources
#'
#' @description Generates one dataset of the simulation study (Section 3.1 and Appendix C of Guo and Zhu). Ten
#'    covariates are drawn from a multivariate normal distribution with unit variances and pairwise correlation 0.1,
#'    and X1-X4 are dichotomized at 0. The outcome is Y = beta0'X + beta1 T + e in the trial and
#'    Y = delta_k + beta0'X + e in external source k, with e ~ N(0, 1).
#'
#' @param n_treat,n_control numbers of treated and control subjects in the trial
#' @param n_hist vector of length K: number of subjects in each external source
#' @param mu_current covariate mean in the trial (all covariates share it)
#' @param mu_hist vector of length K: covariate mean in each external source (sets the covariate shift)
#' @param delta vector of length K: study-level outcome shift of each source (0 = exchangeable given the covariates;
#'    nonzero values create prior-data conflict)
#' @param Vtimer vector of length K: the covariance matrix of source k is that of the trial times Vtimer[k]
#'    (1 in every scenario of the study)
#' @param beta1 treatment effect (0 for type I error, 0.7 for power)
#' @param beta0_com coefficient of the always-observed covariates X1, X3, X5, X7, X9 (0.1 in the study)
#' @param beta0_mis coefficient of the covariates that may be missing, X2, X4, X6, X8, X10 (0.3 in the study)
#' @param trial_bounds NULL (no eligibility criteria: linear selection) or a named list of c(lower, upper) bounds on
#'    continuous covariates (X5-X10), e.g. list(X5 = c(0, Inf), X6 = c(-Inf, 2)). The trial then enrols only subjects
#'    within all bounds (rejection sampling), while the external sources are unrestricted, so the true propensity
#'    score changes abruptly at the bounds (nonlinear selection)
#'
#' The function has no seed argument: call set.seed() before it for reproducible data. simdat_missing() then
#' continues the same random number stream.
#'
#' @return a list with RCT (data frame of the trial, both arms, columns X1-X10, treatment, Y, study), RWD (list of K
#'    data frames, one per external source) and params (the generating parameters)
simdat_complete = function(n_treat, n_control, n_hist,
                           mu_current, mu_hist,
                           delta, Vtimer,
                           beta1, beta0_com, beta0_mis, trial_bounds = NULL)
{
  #fixed parameters
  p = 10  #covariate dimension
  rho = 0.1; sigma2_beta0 = 1  #pairwise covariate correlation / covariate marginal variance
  beta0 = rep(c(beta0_com, beta0_mis), p/2)  #covariates effects. (X2,X4,X6,X8,X10) will have missingness
  sigma_eps = 1  #outcome residual SD
  K = length(n_hist)  #number of external sources (lengths of mu_hist, delta, Vtimer are checked in run_scenario())

  #helper: binarize first 4 covariates at cutpoint 0
  binarize_covariates = function(X) {
    X[,1:4] = ifelse(X[,1:4]>0, 1, 0)
    return(X)
  }
  #helper: generate covariates for a single study
  generate_covariates = function(n, mu_val, V) {
    mu_vec = rep(mu_val, p)
    X = mvrnorm(n = n, mu = mu_vec, Sigma = V)
    X = binarize_covariates(X)
    colnames(X) = paste0("X", 1:p)
    return(X)
  }
  
  #RCT covariates covariance matrix
  V_beta0 = matrix(rho * sigma2_beta0, nrow = p, ncol = p)
  diag(V_beta0) = sigma2_beta0
  
  #1.generate RCT data
  n_current = n_treat + n_control
  if (is.null(trial_bounds)) {
    X_current = generate_covariates(n_current, mu_current, V_beta0)
  } else {
    #eligibility criteria: draw batches and keep rows within all bounds until n_current rows are accepted
    if (!all(names(trial_bounds) %in% paste0("X", 5:p))) stop("trial_bounds must be on continuous covariates X5-X10")
    X_current = NULL
    while (NROW(X_current) < n_current) {
      Xb = generate_covariates(n_current, mu_current, V_beta0)
      keep = rep(TRUE, nrow(Xb))
      for (v in names(trial_bounds)) keep = keep & Xb[, v] >= trial_bounds[[v]][1] & Xb[, v] <= trial_bounds[[v]][2]
      X_current = rbind(X_current, Xb[keep, , drop = FALSE])
    }
    X_current = X_current[seq_len(n_current), , drop = FALSE]
  }
  T_current = c(rep(1, n_treat), rep(0, n_control))  #1=treatment; 0=control
  eps_current = rnorm(n_current, mean = 0, sd = sigma_eps)
  Y_current   = as.numeric(X_current %*% beta0) + beta1 * T_current + eps_current
  current_df = data.frame(X_current, treatment = T_current, Y = Y_current, study = "current")
  
  #2.generate the external sources (controls only)
  historical_list = vector("list", K)
  for (k in 1:K) {
    n_k = n_hist[k]
    mu_k = mu_hist[k]
    delta_k = delta[k]
    Vtimer_k = Vtimer[k]
    
    X_k = generate_covariates(n_k, mu_k, V_beta0*Vtimer_k)
    eps_k = rnorm(n_k, mean = 0, sd = sigma_eps)
    Y_k   = delta_k + as.numeric(X_k %*% beta0) + eps_k
    historical_list[[k]] = data.frame(X_k, treatment = 0, Y = Y_k, study = paste0("historical_", k))
  }
  names(historical_list) = paste0("H", 1:K)
  
  #return
  return(list(RCT = current_df, RWD = historical_list,
    params = list(n_treat = n_treat, n_control = n_control, n_hist = n_hist,
                  mu_current = mu_current, mu_hist = mu_hist,
                  delta = delta, Vtimer = Vtimer,
                  beta1 = beta1, beta0 = beta0,
                  p = p, K = K,
                  rho = rho, sigma2_beta0 = sigma2_beta0,
                  sigma_eps = sigma_eps, trial_bounds = trial_bounds)))
}



#' @description helper function used in simdat_missing().
#'    given target row-level missing rate, calculate per-variable marginal missing probability pi_var,
#'    such that if p_M variables are each missing independently with probability pi_var,
#'    the resulting row-level missing rate equals row_rate.
#'
#' @param row_rate target row-level missing rate (between [0,1])
#' @param p_M number of variables subject to missingness
#' 
#' @return per-variable marginal missing probability pi_var
mcar_pi_from_row_rate = function(row_rate, p_M) {
  if (row_rate == 0) return(0)
  return(1 - (1-row_rate) ^ (1/p_M))
}


#' @description helper function used in simdat_missing() MAR/NMAR scenario. binary search for the value g at which
#'    the average missing probability over the subjects of one external source, mean(miss_prob_fn(g)), equals target.
#'    used twice per dataset:
#'      1. per-variable intercept gamma0_j: miss_prob_fn(g) = P(X_j missing) per subject, target = per-variable rate pi_var
#'      2. common shift of all intercepts: miss_prob_fn(g) = P(row has >= 1 missing) per subject, target = row_rate
#'
#' @param miss_prob_fn function of g returning a vector of per-subject missing probabilities, increasing in g
#' @param target target average missing probability
#' @param lo lower bound for binary search on g
#' @param hi upper bound for binary search on g
#' @param tol convergence tolerance on the average missing probability
#' @param max_iter maximum number of binary search iterations
#' 
#' @return g such that mean(miss_prob_fn(g)) is approximately target
find_gamma0 = function(miss_prob_fn, target,
                       lo = -20, hi = 20, tol = 1e-4, max_iter = 100) {
  if (target == 0) return(-Inf)   # effectively zero missingness
  
  for (i in seq_len(max_iter)) {
    mid = (lo + hi) / 2
    achieved = mean(miss_prob_fn(mid))
    if (abs(achieved - target) < tol) return(mid)
    if (achieved < target) lo = mid else hi = mid
  }
  return((lo + hi) / 2)
}



#' @title Make covariates missing in the external sources
#'
#' @description Makes X2, X4, X6, X8 and X10 partly missing in every external source. The trial stays complete
#'    (Section 3.2 and Appendix D of Guo and Zhu). Mechanisms:
#'      MCAR: each entry is missing with the same probability
#'      MAR:  logit P(X_j missing) = gamma0_j + s + 1.0 * driver_j, with drivers X1->X2, X3->X4, X5->X6, X7->X8, X9->X10
#'      NMAR: logit P(X_j missing) = gamma0_j + s + 1.5 * X_j
#'    outdep = TRUE adds gamma_y * Ytilde (Y standardized within the source) to the MAR or NMAR predictor. The
#'    intercepts gamma0_j are calibrated so that each covariate's average missing probability equals the per-variable
#'    rate implied by row_rate, and the common shift s so that the expected share of subjects with at least one missing
#'    covariate equals row_rate.
#'    Call it right after simdat_complete() without re-seeding, so that the random number stream continues.
#'
#' @param complete_data output of simdat_complete()
#' @param mechanism "MCAR", "MAR" or "NMAR"
#' @param row_rate target share of external subjects with at least one missing covariate (0 = no missingness)
#' @param outdep TRUE adds outcome dependence (MAR or NMAR only). Outcome-dependent MAR remains missing at random,
#'    because the outcome is always observed
#' @param gamma_y log odds ratio of missingness per standard deviation of the outcome (default log(2): the odds of a
#'    covariate being missing double for each one-standard-deviation increase in the outcome)
#'
#' @return complete_data with missing values in RWD, plus miss_params (mechanism, outcome dependence, target and
#'    achieved row-level rates)
simdat_missing = function(complete_data, mechanism, row_rate, outdep = FALSE, gamma_y = log(2)) 
{
  #sanity check
  mechanism = toupper(mechanism)
  if (!mechanism %in% c("MCAR", "MAR", "NMAR"))
    stop("mechanism must be 'MCAR', 'MAR', or 'NMAR'")
  if (row_rate < 0 | row_rate > 1)
    stop("row_rate must be between 0 and 1")
  if (!is.logical(outdep) || length(outdep) != 1 || is.na(outdep))
    stop("outdep must be TRUE or FALSE")
  if (outdep && mechanism == "MCAR")
    stop("outdep = TRUE requires mechanism 'MAR' or 'NMAR' (MCAR cannot depend on the outcome)")
  
  
  #make covariates missing in each external source independently
  miss_vars = c("X2", "X4", "X6", "X8", "X10")
  p_M = length(miss_vars)
  K = length(complete_data$RWD)
  hist_missing = complete_data$RWD
  achieved_row_rate = rep(NA, K)
  
  for (k in seq_len(K)) {
    df = hist_missing[[k]]
    n_k = nrow(df)
    
    #extract the columns that can go missing
    X_miss = as.matrix(df[, miss_vars])   # n_k x p_M matrix
    
    #initialize missingness indicator matrix (0 = observed, 1 = missing)
    M = matrix(0, nrow = n_k, ncol = p_M, dimnames = list(NULL, miss_vars))
    
    #no missingness: skip all
    if (row_rate == 0) {
      achieved_row_rate[k] = 0
      next
    }

    #per-variable marginal missing prob to achieve overall target row rate
    pi_var = mcar_pi_from_row_rate(row_rate, p_M)


    #scenario 1: MCAR
    # each entry missing independently with prob pi_var
    if (mechanism == "MCAR") {
      M = matrix(rbinom(n_k * p_M, size = 1, prob = pi_var),
                 nrow = n_k, ncol = p_M,
                 dimnames = list(NULL, miss_vars))
    }
    
    
    #scenario 2/3: MAR / NMAR (optionally outcome-dependent)
    # MAR : missingness of each variable driven by a distinct always-observed covariate:
    #         X1->X2, X3->X4, X5->X6, X7->X8, X9->X10
    #       logit P(M_j = 1) = gamma0_j + gamma1 * driver_j, gamma1 fixed at 1.0 (moderate dependence)
    # NMAR: missingness of X_j depends on X_j's own (unobserved) value
    #       logit P(M_j = 1) = gamma0_j + gamma1 * X_j, gamma1 fixed at 1.5 (stronger, to make NMAR visibly distinct)
    # outdep = TRUE adds gamma_y * Ytilde to every variable's linear predictor (Ytilde = standardized Y within source;
    #   gamma_y = log(2) by default: odds of missingness double per SD of Y)
    # gamma0_j tuned per variable so each variable's marginal miss rate equals the per-variable rate implied by row_rate,
    # then a common shift of all intercepts so that the expected row-level miss rate equals row_rate (see below)
    if (mechanism %in% c("MAR", "NMAR")) {
      #linear predictor without intercept, n_k x p_M
      if (mechanism == "MAR") {
        gamma1 = 1.0
        eta = gamma1 * as.matrix(df[, c("X1", "X3", "X5", "X7", "X9")])
      } else {
        gamma1 = 1.5
        eta = gamma1 * X_miss
      }
      if (outdep) eta = eta + gamma_y * as.numeric(scale(df[["Y"]]))
      
      #per-variable intercepts
      gamma0 = sapply(seq_len(p_M), function(j) find_gamma0(function(g0) plogis(g0 + eta[,j]), pi_var))
      lp = sweep(eta, 2, gamma0, "+")
      
      #the per-variable missing probabilities are positively related across variables (strongly through the shared
      #outcome term when outdep = TRUE, weakly through the 0.1 covariate correlation otherwise), so missing entries
      #cluster in the same rows and the per-variable calibration undershoots row_rate. shift all intercepts by a
      #common amount so that the expected row-level missing rate equals row_rate
      shift = find_gamma0(function(g) 1 - apply(1 - plogis(lp + g), 1, prod), row_rate)
      lp = lp + shift
      
      for (j in seq_len(p_M)) M[,j] = rbinom(n_k, size = 1, prob = plogis(lp[,j]))
    }
    
    
    #apply missingness: replace observed values with NA where M = 1
    for (j in seq_len(p_M)) {
      var_j = miss_vars[j]
      df[M[,j] == 1, var_j] = NA
    }
    
    #store achieved row-level missing rate for diagnostics
    achieved_row_rate[k] = mean(apply(M, 1, function(r) any(r == 1)))
    attr(df, "achieved_row_rate") = achieved_row_rate[k]
    attr(df, "mechanism") = mechanism
    attr(df, "target_row_rate") = row_rate
    hist_missing[[k]] = df
  }
  
  
  #same structure as complete_data, with missing values in the external sources
  return(list(RCT = complete_data$RCT,
              RWD = hist_missing,
              params = complete_data$params,
              miss_params = list(mechanism = mechanism, outdep = outdep, gamma_y = if (outdep) gamma_y else NA,
                                 target_row_rate = rep(row_rate, K), achieved_row_rate = achieved_row_rate,
                                 miss_vars = miss_vars, p_M = p_M)))
}

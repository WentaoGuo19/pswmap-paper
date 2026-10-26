#' @title Weighting diagnostics of one propensity score method
#'
#' @description Per-source diagnostics recorded in the simulation study (Appendix F of Guo and Zhu), computed from
#'    the weights of a pswmap::psw_weights() object. Unlike pswmap::psw_summaries(), they use information that is
#'    available only in a simulation:
#'    - the missingness pattern of the observed data, also for the oracle, whose weights are fitted to the complete
#'      data. Its relative weight of incomplete subjects is the benchmark for the other methods;
#'    - the true covariate values before missingness, for covariate balance.
#'
#'    For each source k:
#'      ybar_k, se_k : weighted mean and its standard error, as entered into the MAP prior (Rubin's rules under MI)
#'      kish_k       : Kish effective sample size of the weights, (sum w)^2 / sum w^2
#'      relw_inc_k   : share of the weight carried by incomplete subjects / their share of subjects
#'                     (0 = complete-case weighting, NA if the source is complete)
#'      smd_mean_k, smd_max_k : mean and maximum over X1-X10 of |weighted mean in source k - trial mean| / trial SD,
#'                     on the true covariate values
#'    With multiple imputation, kish, relw and smd are averaged over the imputed datasets.
#'
#' @param w a pswmap::psw_weights() object
#' @param incomplete list of K logical vectors: which external subjects have a missing covariate in the observed data
#' @param rwd_full list of K data frames with the true (pre-missingness) covariates of the external sources
#' @param rct_full data frame of the trial (both arms) with its covariates
#'
#' @return named numeric vector: ybar_1..K, se_1..K, kish_1..K, relw_inc_1..K, smd_mean_1..K, smd_max_1..K
weight_diagnostics = function(w, incomplete, rwd_full, rct_full)
{
  xv = paste0("X", 1:10)
  K = length(w$sources)
  s = pswmap::psw_summaries(w)
  m1 = colMeans(rct_full[, xv])
  s1 = apply(rct_full[, xv], 2, sd)
  per_source = sapply(seq_len(K), function(k) {
    W = w$sources[[k]]$w
    inc = incomplete[[k]]
    Xk = as.matrix(rwd_full[[k]][, xv])
    per_imp = apply(W, 2, function(wj) {
      smd = abs(colSums(wj * Xk) / sum(wj) - m1) / s1
      c(kish = sum(wj)^2 / sum(wj^2),
        relw = if (any(inc)) (sum(wj[inc]) / sum(wj)) / mean(inc) else NA_real_,
        smd_mean = mean(smd), smd_max = max(smd))
    })
    rowMeans(matrix(per_imp, nrow = 4, dimnames = list(c("kish", "relw", "smd_mean", "smd_max"), NULL)))
  })
  out = c(s$mean, s$se, per_source["kish", ], per_source["relw", ], per_source["smd_mean", ], per_source["smd_max", ])
  names(out) = paste0(rep(c("ybar_", "se_", "kish_", "relw_inc_", "smd_mean_", "smd_max_"), each = K), seq_len(K))
  out
}

#' @title Simulation evaluation: core performance metrics
#'
#' @description computes the core performance metrics from Zhu et al. (2024) across simulation replications,
#'    plus the rejection rate and the average prior effective sample size (ESS) of the robust PSW-MAP prior.
#'    Metrics: ATE, Bias, %Bias, SE, MSE, EmpVar, Coverage, RejectRate, avg_ESS, computed over the
#'    replications (rows of results) with true treatment effect beta1.
#'    call on the combined results of all chunks of a scenario (rbind), not on single chunks.
#'
#' @param results matrix of per-replicate results (one row per replicate), one method of run_scenario().
#'    Must have columns "mean", "sd", "2.5%", "97.5%" (the ATE posterior mean, sd, and 95% credible interval)
#'    and "ESS"
#' @param beta1 true treatment effect (the target estimand) for this scenario
#'
#' @return metrics: a vector with ATE, Bias, PctBias (NA when beta1 = 0), SE, MSE, EmpVar, Coverage,
#'    RejectRate, avg_ESS, and R (number of replications used, without those whose MAP fit failed).
#'    RejectRate = proportion of replications with the lower 95% credible limit above 0, i.e. posterior
#'    probability of a positive treatment effect above 0.975 (one-sided). it is the type I error rate when
#'    beta1 = 0 and the power when beta1 > 0.
#'    followed by the Monte Carlo standard error (MCSE) of each metric, i.e. its simulation uncertainty from using
#'    R replications (Morris, White & Crowther, 2019, Stat Med, Table 6): MCSE_Bias (also the MCSE of ATE),
#'    MCSE_PctBias, MCSE_SE, MCSE_MSE, MCSE_EmpVar (assumes approximately normal posterior means), MCSE_Coverage,
#'    MCSE_RejectRate and MCSE_avg_ESS
post_eva = function(results, beta1)
{
  #replicates with an estimate (a MAP fit that failed in every attempt of map_prior_guarded() is NA and dropped)
  results = results[!is.na(results[, 'mean']), , drop = FALSE]
  R = nrow(results)

  ATE_r = results[,'mean']
  sd_r = results[,'sd']
  CI_lo_r = results[,'2.5%']
  CI_hi_r = results[,'97.5%']
  ESS_r = results[,'ESS']

  #core metrics, following Zhu et al. (2024)
  ATE = mean(ATE_r)
  Bias = ATE - beta1
  PctBias = if (beta1 == 0) NA_real_ else Bias / beta1 * 100  #undefined when the true effect is 0
  SE = mean(sd_r)
  MSE = mean((ATE_r - beta1)^2)
  EmpVar = sum((ATE_r - ATE)^2) / (R - 1)
  Coverage = mean(beta1 >= CI_lo_r & beta1 <= CI_hi_r)

  #rejection rate: type I error (beta1 = 0) or power (beta1 > 0)
  RejectRate = mean(CI_lo_r > 0)

  avg_ESS = mean(ESS_r)
  metrics = c(ATE = ATE, Bias = Bias, PctBias = PctBias, SE = SE, MSE = MSE, EmpVar = EmpVar,
              Coverage = Coverage, RejectRate = RejectRate, avg_ESS = avg_ESS, R = R)

  #Monte Carlo standard errors (Morris, White & Crowther, 2019)
  MCSE_Bias = sqrt(EmpVar / R)
  mcse = c(MCSE_Bias = MCSE_Bias,
           MCSE_PctBias = if (beta1 == 0) NA_real_ else MCSE_Bias / abs(beta1) * 100,
           MCSE_SE = sd(sd_r) / sqrt(R),
           MCSE_MSE = sqrt(sum(((ATE_r - beta1)^2 - MSE)^2) / (R * (R - 1))),
           MCSE_EmpVar = EmpVar * sqrt(2 / (R - 1)),
           MCSE_Coverage = sqrt(Coverage * (1 - Coverage) / R),
           MCSE_RejectRate = sqrt(RejectRate * (1 - RejectRate) / R),
           MCSE_avg_ESS = sd(ESS_r) / sqrt(R))
  metrics = c(metrics, mcse)

  return(metrics)
}

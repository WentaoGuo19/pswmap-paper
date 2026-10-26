# Combines the chunk files of each scenario and computes the performance measures (Appendix F of Guo and Zhu).
# Run with Rscript from inside the folder that holds the files sid_<sid>_chunk_<chunk>.rds written by
# cluster/task_main.R or cluster/task_pooled.R, e.g. from the repository root:
#   cd sim_main && N_REPS=1000 Rscript ../evaluation/evaluation.R
# post_eva.R (next to this script) and R/scenarios.R are found from this script's own location.
# The script evaluates whatever values of c and methods the chunk files contain. Only scenarios with all N_REPS
# replicates present exactly once are evaluated. The others are reported in status.csv.
#
# Output:
#   summary.csv: one row per scenario x c x method, with the scenario settings, the performance measures of
#     post_eva() and the averaged weighting diagnostics (diag_summary() below; NA for rct_only)
#   status.csv: per scenario, the number of chunk files and replicates found, whether it is complete, and the run
#     time per replicate (sec_per_rep)


#folder of this script (from Rscript's --file argument)
script_file = sub('^--file=', '', grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(script_file) != 1) stop('run this script with Rscript (e.g. Rscript ../evaluation/evaluation.R), not source()')
script_file = gsub('~+~', ' ', script_file, fixed = TRUE)   #Rscript passes spaces in the script path as ~+~
script_dir = dirname(normalizePath(script_file))

#scenario table (416 main + 48 sensitivity scenarios) and the performance measures
source(file.path(script_dir, '..', 'R', 'scenarios.R'))
source(file.path(script_dir, 'post_eva.R'))
scenarios = make_scenarios()
n_reps = as.integer(Sys.getenv('N_REPS', '1000'))   #replicates per scenario


#covariate-adjusted target of each source's weighted mean: delta_k + E_trial[beta0'X], where E_trial[beta0'X] is the
#scenario's target_mean column (1.873 without trial eligibility criteria, 2.0175 with them; see scenarios.R)


#' averaged per-source weighting diagnostics for one method (matrix of combined replicates, see run_scenario())
#'    bias_src_k : mean of ybar_k minus its target delta_k + target_mean (residual covariate imbalance / missing-data bias)
#'    mae_src    : mean over sources of |bias_src_k|
#'    se_src, kish : mean over sources and replications of the SE of ybar_k and of the Kish ESS
#'    relw_inc   : mean relative weight of incomplete subjects (1 = like complete subjects, 0 = complete-case)
#'    smd_mean, smd_max : mean over sources and replications of the mean / max absolute SMD against the trial
diag_summary = function(res, delta, target_mean)
{
  K = length(delta)
  cols = function(prefix) paste0(prefix, seq_len(K))
  if (!all(cols("ybar_") %in% colnames(res)) || all(is.na(res[, cols("ybar_")]))) {
    out = rep(NA_real_, K + 6)
  } else {
    #na.rm: a replicate whose propensity score step failed entirely is an all-NA row and is skipped. Replicates
    #whose MAP fit failed keep their diagnostics, so they are still included here
    bias = colMeans(res[, cols("ybar_"), drop = FALSE], na.rm = TRUE) - (delta + target_mean)
    out = c(bias, mean(abs(bias)),
            mean(res[, cols("se_")], na.rm = TRUE), mean(res[, cols("kish_")], na.rm = TRUE),
            mean(res[, cols("relw_inc_")], na.rm = TRUE),   #NaN when no subject is incomplete (0% missing)
            mean(res[, cols("smd_mean_")], na.rm = TRUE), mean(res[, cols("smd_max_")], na.rm = TRUE))
  }
  names(out) = c(cols("bias_src_"), "mae_src", "se_src", "kish", "relw_inc", "smd_mean", "smd_max")
  out[is.nan(out)] = NA
  return(out)
}


#chunk files and their scenario ids
files = list.files(pattern = '^sid_[0-9]+_chunk_[0-9]+\\.rds$')
file_sid = as.integer(sub('^sid_([0-9]+)_chunk_[0-9]+\\.rds$', '\\1', files))


#combine chunks and evaluate, scenario by scenario
status = data.frame(sid = scenarios$sid, n_files = 0L, n_reps = 0L, complete = FALSE, sec_per_rep = NA_real_)
summary_list = list()
for (sid in scenarios$sid) {
  f = files[file_sid == sid]
  if (length(f) == 0) next

  #read chunks; all must contain the same c values (c32, c16, ...) and, within each, the same methods
  chunks = lapply(f, readRDS)
  c_names = names(chunks[[1]])
  methods = names(chunks[[1]][[1]])
  if (!all(sapply(chunks, function(x) identical(names(x), c_names) && all(sapply(x, function(y) identical(names(y), methods))))))
    stop('scenario ', sid, ': chunk files contain different c values or methods')

  #stack the replications of all chunks, per c value and method
  combined = lapply(setNames(c_names, c_names), function(cn)
    lapply(setNames(methods, methods), function(m) do.call(rbind, lapply(chunks, function(x) x[[cn]][[m]]))))
  reps = as.integer(rownames(combined[[1]][[1]]))
  if (anyDuplicated(reps)) stop('scenario ', sid, ': the same replication appears in more than one chunk file')
  status$n_files[sid] = length(f)
  status$n_reps[sid] = length(reps)
  sec = sapply(chunks, function(x) if (is.null(attr(x, "elapsed_sec"))) NA_real_ else attr(x, "elapsed_sec"))
  status$sec_per_rep[sid] = sum(sec) / length(reps)   #run time per replication (NA for files saved without timing)
  status$complete[sid] = setequal(reps, seq_len(n_reps))
  if (!status$complete[sid]) next

  #metrics and averaged diagnostics per c value and method
  delta = unlist(scenarios[sid, c("d1", "d2", "d3")])
  for (cn in c_names) {
    met = t(sapply(combined[[cn]], function(res) c(post_eva(res, beta1 = scenarios$beta1[sid]), diag_summary(res, delta, scenarios$target_mean[sid]))))
    summary_list[[length(summary_list) + 1]] =
      cbind(scenarios[rep(sid, nrow(met)), ], c = as.numeric(sub("^c", "", cn)), method = rownames(met), met, row.names = NULL)
  }
}


#save
write.csv(status, 'status.csv', row.names = FALSE)
cat(sum(status$complete), 'of', nrow(scenarios), 'scenarios complete;',
    sum(status$n_files > 0 & !status$complete), 'incomplete;', sum(status$n_files == 0), 'with no files\n')
if (length(summary_list) > 0) {
  summary = do.call(rbind, summary_list)
  write.csv(summary, 'summary.csv', row.names = FALSE)
  cat('summary.csv:', nrow(summary), 'rows (scenario x c x method)\n')
}

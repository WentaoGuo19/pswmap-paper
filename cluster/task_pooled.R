# Pooled propensity score comparison (Slurm array task; submitted with cluster/run_pooled.sh from the repository root).
# The 16 main-grid scenarios with MAR missingness at a 20% row-level rate are analysed with ONE propensity score
# fitted to the trial against all external sources pooled (methods <learner>_pooled_<strategy>). Scenario ids, and
# therefore seeds, are those of the main simulation, so the pooled and source-specific analyses use identical
# datasets (paired comparison). Strategies: oracle, complete-case and MI with the default engine. Native handling is
# left out because it converges to complete-case weighting, and the other imputation engines because they behave
# like the default engine in the main simulation (Appendix F of Guo and Zhu).
# Environment variables as in task_main.R; RESULTS_DIR defaults to sim_pooled.

options(java.parameters = "-Xmx8g")
library(bartMachine)
set_bart_machine_num_cores(4)
library(pswmap)

for (f in c("simulate_data.R", "scenarios.R", "diagnostics.R", "map_guard.R", "run_scenario.R")) source(file.path("R", f))
scenarios = make_scenarios()
pooled_sids = scenarios$sid[scenarios$grid == "main" & scenarios$pattern == "MAR" & scenarios$row_rate == 0.2]
stopifnot(length(pooled_sids) == 16)

n_chunks = as.integer(Sys.getenv("N_CHUNKS", "10"))
n_reps = as.integer(Sys.getenv("N_REPS", "1000"))
results_dir = Sys.getenv("RESULTS_DIR", "sim_pooled")
if (is.na(n_reps) || is.na(n_chunks) || n_chunks < 1 || n_chunks > n_reps) stop("need 1 <= N_CHUNKS <= N_REPS")
task = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) + as.integer(Sys.getenv("TASK_OFFSET", "0"))
if (is.na(task) || task < 1 || task > length(pooled_sids) * n_chunks)
  stop("task must be between 1 and ", length(pooled_sids) * n_chunks, " (got ", task, ")")
sid = pooled_sids[(task - 1) %/% n_chunks + 1]
chunk = (task - 1) %% n_chunks + 1
reps = split(seq_len(n_reps), sort(rep_len(seq_len(n_chunks), n_reps)))[[chunk]]
s = scenarios[sid, ]
cat("pooled score | task", task, "| scenario", sid, "| chunk", chunk, "of", n_chunks, "| replicates", min(reps), "-",
    max(reps), "\n")

dir.create(results_dir, showWarnings = FALSE)
out_file = file.path(results_dir, paste0("sid_", sid, "_chunk_", chunk, ".rds"))
if (file.exists(out_file)) { cat("already done:", out_file, "\n"); quit(save = "no") }

t0 = proc.time()[["elapsed"]]
result = run_scenario(n_treat = 100, n_control = 50, n_hist = rep(200, 3), mu_current = 1,
                      mu_hist = c(s$mu_h1, s$mu_h2, s$mu_h3), delta = c(s$d1, s$d2, s$d3), Vtimer = c(s$v1, s$v2, s$v3),
                      c = as.numeric(strsplit(s$c_values, ",")[[1]]), trial_bounds = parse_bounds(s$trial_bounds),
                      beta1 = s$beta1, beta0_com = 0.1, beta0_mis = 0.3, mechanism = s$mechanism,
                      row_rate = s$row_rate, outdep = s$outdep, gamma_y = if (s$outdep) s$gamma_y else log(2),
                      mi_engines = "default", native = FALSE, ps_scope = "pooled", PS_cap = 0.95, reps = reps, sid = sid)
attr(result, "elapsed_sec") = proc.time()[["elapsed"]] - t0
cat("finished in", round(attr(result, "elapsed_sec")), "s,", round(attr(result, "elapsed_sec") / length(reps), 1),
    "s per replicate\n")
saveRDS(result, file = out_file)

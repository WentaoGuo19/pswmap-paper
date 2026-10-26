# Main simulation (Slurm array task; submitted with cluster/run_main.sh from the repository root).
# 464 scenarios (416 main grid + 48 sensitivity scenarios, see R/scenarios.R), 1,000 replicates each, split into chunks.
# Task t runs chunk ((t - 1) %% n_chunks) + 1 of scenario ((t - 1) %/% n_chunks) + 1 and saves
# <RESULTS_DIR>/sid_<sid>_chunk_<chunk>.rds. Evaluate the finished folder with evaluation/evaluation.R.
# Environment variables: N_REPS (default 1000), N_CHUNKS (default 10), RESULTS_DIR (default sim_main), TASK_OFFSET
# (added to the task id, for submissions split across several arrays).

#Java heap for bartMachine: set it before the package is loaded, and keep it below the job's memory (run_main.sh: 14G).
#On some systems Java also needs "--add-modules=jdk.incubator.vector", e.g.
#options(java.parameters = c("-Xmx8g", "--add-modules=jdk.incubator.vector"))
options(java.parameters = "-Xmx8g")
library(bartMachine)
set_bart_machine_num_cores(4)   #one chain per core
library(pswmap)

for (f in c("simulate_data.R", "scenarios.R", "diagnostics.R", "map_guard.R", "run_scenario.R")) source(file.path("R", f))
scenarios = make_scenarios()

n_chunks = as.integer(Sys.getenv("N_CHUNKS", "10"))
n_reps = as.integer(Sys.getenv("N_REPS", "1000"))
results_dir = Sys.getenv("RESULTS_DIR", "sim_main")
if (is.na(n_reps) || is.na(n_chunks) || n_chunks < 1 || n_chunks > n_reps) stop("need 1 <= N_CHUNKS <= N_REPS")
task = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) + as.integer(Sys.getenv("TASK_OFFSET", "0"))
if (is.na(task) || task < 1 || task > nrow(scenarios) * n_chunks)
  stop("task must be between 1 and ", nrow(scenarios) * n_chunks, " (got ", task, ")")
sid = (task - 1) %/% n_chunks + 1
chunk = (task - 1) %% n_chunks + 1
reps = split(seq_len(n_reps), sort(rep_len(seq_len(n_chunks), n_reps)))[[chunk]]   #contiguous, near-equal chunks
s = scenarios[sid, ]
cat("task", task, "| scenario", sid, "| chunk", chunk, "of", n_chunks, "| replicates", min(reps), "-", max(reps), "\n")
cat(s$grid, "|", s$shift, s$selection, s$heterogeneity, "|", s$pattern, s$row_rate, "| beta1", s$beta1,
    "| c", s$c_values, "| gamma_y", s$gamma_y, "| trial bounds", s$trial_bounds, "\n")

#a finished chunk is skipped, so failed or timed-out tasks can simply be resubmitted
dir.create(results_dir, showWarnings = FALSE)
out_file = file.path(results_dir, paste0("sid_", sid, "_chunk_", chunk, ".rds"))
if (file.exists(out_file)) { cat("already done:", out_file, "\n"); quit(save = "no") }

t0 = proc.time()[["elapsed"]]
result = run_scenario(n_treat = 100, n_control = 50, n_hist = rep(200, 3), mu_current = 1,
                      mu_hist = c(s$mu_h1, s$mu_h2, s$mu_h3), delta = c(s$d1, s$d2, s$d3), Vtimer = c(s$v1, s$v2, s$v3),
                      c = as.numeric(strsplit(s$c_values, ",")[[1]]), trial_bounds = parse_bounds(s$trial_bounds),
                      beta1 = s$beta1, beta0_com = 0.1, beta0_mis = 0.3, mechanism = s$mechanism,
                      row_rate = s$row_rate, outdep = s$outdep, gamma_y = if (s$outdep) s$gamma_y else log(2),
                      PS_cap = 0.95, reps = reps, sid = sid)
attr(result, "elapsed_sec") = proc.time()[["elapsed"]] - t0
cat("finished in", round(attr(result, "elapsed_sec")), "s,", round(attr(result, "elapsed_sec") / length(reps), 1),
    "s per replicate\n")
saveRDS(result, file = out_file)

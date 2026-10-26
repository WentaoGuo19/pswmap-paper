# Sample-size experiment (Slurm array task; submitted with cluster/run_size.sh from the repository root).
# 8 cells = size factor f in {1, 2, 4, 8} x {MCAR, outcome-dependent MAR}, both at a 20% row-level missing rate,
# 1,000 replicates each, split into chunks. Task t runs chunk ((t - 1) %% n_chunks) + 1 of cell ((t - 1) %/% n_chunks) + 1;
# cells 1-4 are MCAR with f = 1, 2, 4, 8, cells 5-8 outcome-dependent MAR with f = 1, 2, 4, 8. Output:
# <RESULTS_DIR>/cell_<q>_chunk_<k>.rds; evaluate with evaluation/evaluate_size.R.
# Environment variables as in task_main.R; RESULTS_DIR defaults to sim_size.

options(java.parameters = "-Xmx8g")
library(bartMachine)
set_bart_machine_num_cores(4)
library(pswmap)

for (f in c("simulate_data.R", "diagnostics.R", "run_size_cell.R")) source(file.path("R", f))
cells = expand.grid(f = c(1, 2, 4, 8), mechanism = c("MCAR", "MAR"), stringsAsFactors = FALSE)
cells$outdep = cells$mechanism == "MAR"
cells$cell = seq_len(nrow(cells))

n_chunks = as.integer(Sys.getenv("N_CHUNKS", "10"))
n_reps = as.integer(Sys.getenv("N_REPS", "1000"))
results_dir = Sys.getenv("RESULTS_DIR", "sim_size")
if (is.na(n_reps) || is.na(n_chunks) || n_chunks < 1 || n_chunks > n_reps) stop("need 1 <= N_CHUNKS <= N_REPS")
task = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) + as.integer(Sys.getenv("TASK_OFFSET", "0"))
if (is.na(task) || task < 1 || task > nrow(cells) * n_chunks)
  stop("task must be between 1 and ", nrow(cells) * n_chunks, " (got ", task, ")")
cl = cells[(task - 1) %/% n_chunks + 1, ]
chunk = (task - 1) %% n_chunks + 1
reps = split(seq_len(n_reps), sort(rep_len(seq_len(n_chunks), n_reps)))[[chunk]]
cat("sample size | task", task, "| cell", cl$cell, "(f =", cl$f, ",", cl$mechanism,
    if (cl$outdep) "outcome-dependent", ") | chunk", chunk, "of", n_chunks, "| replicates", min(reps), "-", max(reps), "\n")

dir.create(results_dir, showWarnings = FALSE)
out_file = file.path(results_dir, paste0("cell_", cl$cell, "_chunk_", chunk, ".rds"))
if (file.exists(out_file)) { cat("already done:", out_file, "\n"); quit(save = "no") }

t0 = proc.time()[["elapsed"]]
result = run_size_cell(f = cl$f, mechanism = cl$mechanism, outdep = cl$outdep, cell = cl$cell, reps = reps)
attr(result, "elapsed_sec") = proc.time()[["elapsed"]] - t0
cat("finished in", round(attr(result, "elapsed_sec")), "s,", round(attr(result, "elapsed_sec") / length(reps), 1),
    "s per replicate\n")
saveRDS(result, file = out_file)

# Evaluates the sample-size experiment (cluster/task_size.R).
# Run with Rscript from inside the folder that holds the files cell_<q>_chunk_<k>.rds, e.g. from the repository root:
#   cd sim_size && N_REPS=1000 Rscript ../evaluation/evaluate_size.R
# It combines the chunk files, checks that every cell has replicates 1..N_REPS, and writes summary_size.csv:
# one row per cell x method with the mean and Monte Carlo SE over replicates of
#   relw      relative weight of incomplete subjects (mean over the 3 sources)
#   ex        excess source bias: weighted mean minus that of the oracle with the same learner, mean over the 3 sources
#   kish      Kish effective sample size (mean over the 3 sources)

n_reps = as.integer(Sys.getenv('N_REPS', '1000'))
files = list.files(".", pattern = "^cell_[0-9]+_chunk_[0-9]+\\.rds$")
if (length(files) == 0) stop("no cell_*_chunk_*.rds files here; run from inside the results folder")
d = do.call(rbind, lapply(files, readRDS))

#completeness: every cell needs replicates 1..n_reps, each with all methods
reps_by_cell = tapply(d$rep, d$cell, function(r) sort(unique(r)))
for (cl in names(reps_by_cell)) if (!identical(reps_by_cell[[cl]], seq_len(n_reps)))
  warning("cell ", cl, ": ", length(reps_by_cell[[cl]]), " of ", n_reps, " replicates present")
if (anyDuplicated(d[c("cell", "rep", "method")])) stop("duplicated cell x replicate x method rows")

K = 3
d$learner = sub("_.*", "", d$method)
d$strategy = sub("^[a-z]+_", "", d$method)
d$relw = rowMeans(d[paste0("relw_inc_", 1:K)])
d$kish = rowMeans(d[paste0("kish_", 1:K)])
d$ybar = rowMeans(d[paste0("ybar_", 1:K)])
o = d[d$strategy == "oracle", c("cell", "rep", "learner", "ybar")]; names(o)[4] = "ybar_or"
d = merge(d, o, by = c("cell", "rep", "learner")); d$ex = d$ybar - d$ybar_or

agg = function(v) aggregate(d[[v]], d[c("cell", "f", "setting", "learner", "strategy")],
                            function(x) c(mean = mean(x), mcse = sd(x) / sqrt(length(x)), R = length(x)))
s = Reduce(function(a, b) merge(a, b, by = c("cell", "f", "setting", "learner", "strategy")),
           lapply(c("relw", "ex", "kish"), function(v) { a = agg(v); x = a$x; a$x = NULL
             a[[v]] = x[, "mean"]; a[[paste0("MCSE_", v)]] = x[, "mcse"]; if (v == "relw") a$R = x[, "R"]; a }))
s = s[order(s$setting, s$learner, s$strategy, s$f), ]
write.csv(s, "summary_size.csv", row.names = FALSE)
cat("wrote summary_size.csv:", nrow(s), "rows,", length(unique(d$cell)), "cells,", length(files), "chunk files\n")

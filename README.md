# pswmap-paper

Code and results of the simulation study in Guo and Zhu, *Propensity Score Integrated MAP Prior with Historical
Missingness*.

The method itself is implemented in the R package [pswmap](https://github.com/WentaoGuo19/pswmap). This repository
contains everything around it:

- generating the simulated trials and external sources;
- running the 464 scenarios on a computing cluster;
- evaluating the replicates;
- rebuilding every results table and figure of the paper from the evaluated results.

## Contents

| Folder | Contents |
|---|---|
| `R/` | `simulate_data.R` (trial, external sources, missing covariates), `scenarios.R` (the scenario grid), `run_scenario.R` (all methods for the replicates of one scenario), `run_size_cell.R` (sample-size experiment), `diagnostics.R` (weighting diagnostics), `map_guard.R` (time limit for the MCMC fits) |
| `cluster/` | For the main simulation, the pooled propensity score comparison and the sample-size experiment: Slurm job scripts (`run_main.sh`, `run_pooled.sh`, `run_size.sh`) and the R script that each array task runs (`task_main.R`, `task_pooled.R`, `task_size.R`), plus a throttled submitter (`submit_main.sh`) |
| `evaluation/` | Combine replicate files into performance measures (`evaluation.R`, `post_eva.R`, `evaluate_size.R`) |
| `results/` | The evaluated results, which form the full-results supplement of the paper (see below) |
| `paper/` | `make_tables.R` and `make_figures.R`, which rebuild Tables 3-5 and Figures 1, 2 and S1 into `paper/tables/` and `paper/figures/` |
| `examples/` | `run_one_scenario.R`: a few replicates of one scenario on a laptop |

## Requirements

- R (>= 4.1) and the package **pswmap** 1.0.0: `remotes::install_github("WentaoGuo19/pswmap@v1.0.0")`. It installs RBesT
  and mice.
- **bartMachine** (which needs Java), **ranger** and **xgboost** for the propensity score learners, and **MASS** for
  the data generation.
- **ggplot2** and **patchwork** for the figures.

The study was run with R 4.4.0, bartMachine 1.3.4.1, ranger 0.18.0, xgboost 3.2.1.1, mice 3.19.0 and RBesT 1.10.0.

## Rebuild the tables and figures

The evaluated results are included, so the tables and figures of the paper can be rebuilt in seconds. From the
repository root:

```bash
Rscript paper/make_tables.R     # Tables 3-5 -> paper/tables/*.tex
Rscript paper/make_figures.R    # Figures 1, 2, S1 -> paper/figures/*.pdf
```

Both scripts also print the numbers quoted in the text.

## Rerun the simulation

All scripts run from the repository root. Each task writes one file per chunk of replicates, and a finished chunk is
skipped, so failed or timed-out tasks can simply be resubmitted.

**Main simulation.** 464 scenarios × 1,000 replicates, split into 10 chunks of 100 replicates: 4,640 array tasks.

```bash
mkdir -p logs_main
bash cluster/submit_main.sh              # or: sbatch --export=NONE cluster/run_main.sh
cd sim_main && N_REPS=1000 Rscript ../evaluation/evaluation.R    # -> summary.csv, status.csv
```

**Pooled propensity score comparison.** 16 scenarios, 160 tasks.

```bash
mkdir -p logs_pooled && sbatch --export=NONE cluster/run_pooled.sh
cd sim_pooled && N_REPS=1000 Rscript ../evaluation/evaluation.R
```

**Sample-size experiment.** 8 cells, 80 tasks. See `cluster/run_size.sh` for submitting each size with its own time
limit.

```bash
mkdir -p logs_size && sbatch --export=NONE cluster/run_size.sh
cd sim_size && N_REPS=1000 Rscript ../evaluation/evaluate_size.R   # -> summary_size.csv
```

**Cluster settings.** The job scripts contain the account, partition, module and conda lines of the Alpine cluster
(University of Colorado Boulder / Colorado State University). Replace them with those of your cluster. On Alpine:

- a main-simulation chunk took about 10 hours (at most 14);
- the main simulation used about 40,000 task-hours with 4 cores each.

**Without a cluster.** `Rscript examples/run_one_scenario.R` runs two replicates of one scenario in a few minutes.

## Seeds and reproducibility

Replicate *r* of scenario *sid* generates its data after `set.seed(10000 * sid + r)`. The sample-size experiment uses
`set.seed(10000000 + 10000 * cell + r)`. All later random steps continue the same stream: the imputations, random
forests, XGBoost folds, MCMC sampling and posterior draws. On the same platform, with the same versions of R and the
packages, a rerun gives the same numbers for every learner except BART.

A rerun reproduces the files in `results/` up to Monte Carlo error.

## Results files

These form the full-results supplement of the paper.

- **`summary_main.csv`.** 42,800 rows: 464 scenarios × the values of the heterogeneity constant `c` × 25 methods.
  - Methods are named `<learner>_<strategy>`:
    - learner: `logit`, `bart`, `rf`, `xgb`;
    - strategy: `oracle` (complete data), `cc` (complete-case), `native`, or `mi_default`, `mi_cart`, `mi_rf` (multiple imputation by engine);
    - plus `rct_only` (the trial alone).
  - Columns:
    - the scenario settings (see `R/scenarios.R`);
    - the performance measures of the treatment effect: `ATE`, `Bias`, `PctBias`, `SE`, `MSE`, `EmpVar`, `Coverage`, `RejectRate` (type I error when `beta1 = 0`, power otherwise), `avg_ESS` (prior effective sample size) and `R` (replicates used), each with its Monte Carlo standard error (`MCSE_*`);
    - the averaged weighting diagnostics: `bias_src_1`-`bias_src_3` (bias of each source's weighted mean), `mae_src`, `se_src`, `kish` (Kish effective sample size), `relw_inc` (relative weight of incomplete subjects), `smd_mean` and `smd_max` (covariate balance).
- **`summary_pooled.csv`.** The same columns for the pooled propensity score comparison: methods `<learner>_pooled_<strategy>` and `rct_only`. The scenarios and seeds are those of the main simulation, so the two files are paired.
- **`summary_size.csv`.** The sample-size experiment: one row per cell × learner × strategy.
  - `f`: size factor;
  - `setting`: `mcar`, or `outdep` for outcome-dependent MAR;
  - `relw`: relative weight of incomplete subjects;
  - `ex`: excess source bias beyond the oracle;
  - `kish`: Kish effective sample size;
  - each with its Monte Carlo standard error, and `R`.

## Citation

Guo W, Zhu AY. Propensity score integrated MAP prior with historical missingness. Manuscript submitted.

## License

GPL (>= 3). See [LICENSE.md](LICENSE.md).

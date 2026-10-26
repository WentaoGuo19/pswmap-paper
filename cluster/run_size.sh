#!/bin/bash -l
# Sample-size experiment: one array task per cell x chunk (cluster/task_size.R).
# Submit from the repository root, e.g.  mkdir -p logs_size && sbatch --export=NONE cluster/run_size.sh
# Pass settings as --export=VAR=value (e.g. --export=RESULTS_DIR=sim_size,N_REPS=1000), never --export=ALL: the job must
# start from a clean login environment so that the module system works on the compute node.
#
# SITE-SPECIFIC: the account, partition, QOS, module and conda lines below are those of the Alpine cluster
# (University of Colorado Boulder / Colorado State University). Replace them with your cluster's equivalents.
#SBATCH --account=csu93_alpine1
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=pswmap_size
#SBATCH --output=logs_size/task_%A_%a.out
#SBATCH --array=1-80   # 8 cells x 10 chunks; cells 1-4 MCAR, 5-8 outcome-dependent MAR, each with f = 1, 2, 4, 8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4      # BART runs one chain per core, so more cores do not speed up a task
#SBATCH --mem=14G               # peak use about 4 GB; must exceed the Java heap set in the R script (8 GB)
#SBATCH --time=07:00:00   # enough for f = 8; smaller sizes need less, see below

# Time per 100-replicate chunk grows with the size factor (on Alpine roughly 0.6, 1.1, 2.0 and 3.9 h for f = 1, 2, 4,
# 8). To request only what each size needs, submit the sizes separately:
#   sbatch --export=NONE --array=1-10,41-50  --time=01:30:00 cluster/run_size.sh   # f = 1
#   sbatch --export=NONE --array=11-20,51-60 --time=02:30:00 cluster/run_size.sh   # f = 2
#   sbatch --export=NONE --array=21-30,61-70 --time=04:00:00 cluster/run_size.sh   # f = 4
#   sbatch --export=NONE --array=31-40,71-80 --time=07:00:00 cluster/run_size.sh   # f = 8

module purge
module load gcc/11.2.0          # C++ runtime needed by xgboost
module load anaconda
conda activate /curc/sw/anaconda3/2019.03/envs/rstudio
module load R/4.4.0
echo "node: $(hostname)"; module list 2>&1
Rscript cluster/task_size.R

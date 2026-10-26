#!/bin/bash -l
# Pooled propensity score comparison: one array task per scenario x chunk (cluster/task_pooled.R).
# Submit from the repository root, e.g.  mkdir -p logs_pooled && sbatch --export=NONE cluster/run_pooled.sh
# Pass settings as --export=VAR=value (e.g. --export=RESULTS_DIR=sim_pooled,N_REPS=1000), never --export=ALL: the job must
# start from a clean login environment so that the module system works on the compute node.
#
# SITE-SPECIFIC: the account, partition, QOS, module and conda lines below are those of the Alpine cluster
# (University of Colorado Boulder / Colorado State University). Replace them with your cluster's equivalents.
#SBATCH --account=csu93_alpine1
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=pswmap_pooled
#SBATCH --output=logs_pooled/task_%A_%a.out
#SBATCH --array=1-160   # 16 scenarios x 10 chunks
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4      # BART runs one chain per core, so more cores do not speed up a task
#SBATCH --mem=14G               # peak use about 4 GB; must exceed the Java heap set in the R script (8 GB)
#SBATCH --time=05:30:00   # about 3.7 h per 100-replicate chunk on Alpine

module purge
module load gcc/11.2.0          # C++ runtime needed by xgboost
module load anaconda
conda activate /curc/sw/anaconda3/2019.03/envs/rstudio
module load R/4.4.0
echo "node: $(hostname)"; module list 2>&1
Rscript cluster/task_pooled.R

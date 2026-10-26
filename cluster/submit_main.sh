#!/bin/bash
# Submits the main simulation in groups, without exceeding a per-user limit on queued array tasks (each array task
# counts as one job on many clusters). Run from the repository root.
#
# Each of the 10 chunks is submitted as three groups with their own time limit, because scenarios without missing
# covariates and the sensitivity scenarios run much faster:
#   miss: 400 main-grid scenarios with missing covariates   --time 14:30:00
#   sens:  48 sensitivity scenarios                        --time 06:30:00
#   none:  16 main-grid scenarios without missing covariates --time 02:30:00
# A group is submitted only when (tasks in the queue) + (its tasks) <= LIMIT. Submitted groups are recorded in the
# ledger submitted_main.txt ("chunk group jobid date"), so no group is submitted twice.
#
# Usage:
#   bash cluster/submit_main.sh            one pass: submit every pending group that fits now
#   nohup bash cluster/submit_main.sh --loop > submit_main.log 2>&1 &
#                                          repeat every 30 minutes until all 30 groups are submitted
# LIMIT (default 1000), WAIT (seconds between passes, default 1800), RESULTS_DIR (default sim_main) and CHUNKS
# (default "1 2 3 4 5 6 7 8 9 10") can be set, e.g.  LIMIT=900 bash cluster/submit_main.sh
LIMIT=${LIMIT:-1000}
WAIT=${WAIT:-1800}
RESULTS_DIR=${RESULTS_DIR:-sim_main}
CHUNKS=${CHUNKS:-"1 2 3 4 5 6 7 8 9 10"}
LEDGER=submitted_main.txt
N_GROUPS=$(( $(echo $CHUNKS | wc -w) * 3 ))
touch "$LEDGER"
mkdir -p logs_main

one_pass() {
  for ch in $CHUNKS; do
    for grp in miss sens none; do
      grep -q "^$ch $grp " "$LEDGER" && continue
      # task = 10 * (sid - 1) + chunk; sids 1-416 form the main grid (every 26th, starting at sid 1, has no missing
      # covariates) and sids 417-464 the sensitivity grid
      case $grp in
        miss) list=$(seq $ch 10 4160 | awk -v c=$ch '($0-c)%260!=0' | paste -sd, -); t=14:30:00 ;;
        sens) list=$(seq $((4160+ch)) 10 4640 | paste -sd, -); t=06:30:00 ;;
        none) list=$(seq $ch 260 4160 | paste -sd, -); t=02:30:00 ;;
      esac
      n=$(echo "$list" | tr ',' '\n' | wc -l)
      queued=$(squeue --me -r -h | wc -l)
      if [ $((queued + n)) -le "$LIMIT" ]; then
        if out=$(sbatch --array="$list" --time=$t --export=RESULTS_DIR=$RESULTS_DIR cluster/run_main.sh); then
          echo "$ch $grp ${out##* } $(date '+%F %T')" >> "$LEDGER"
          echo "submitted chunk $ch $grp ($n tasks, --time $t): $out"
        else
          echo "sbatch failed for chunk $ch $grp: $out"
        fi
      fi
    done
  done
  echo "$(date '+%F %T'): $(wc -l < "$LEDGER") of $N_GROUPS groups submitted; $(squeue --me -r -h | wc -l) tasks in the queue"
}

one_pass
if [ "$1" = "--loop" ]; then
  while [ "$(wc -l < "$LEDGER")" -lt "$N_GROUPS" ]; do sleep "$WAIT"; one_pass; done
  echo "all $N_GROUPS groups submitted"
fi

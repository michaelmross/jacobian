#!/bin/bash
# resume.sh CELL... — re-run named cells with a bigger budget, WITH LOGGING.
#   MIRACLE_TIMEOUT=600 ./resume.sh "2 4 2" "3 3 3"
#   tail -F resume.log
# Output goes to both the screen and resume.log (appended, timestamped), so
# a long run is not lost to the scrollback buffer.
set -u
export MIRACLE_TIMEOUT=${MIRACLE_TIMEOUT:-300}
export MIRACLE_SLOW_WARN=${MIRACLE_SLOW_WARN:-30}
LOG=${RESUME_LOG:-resume.log}
{
  echo "=============================================================="
  echo " resume run  $(date)   MIRACLE_TIMEOUT=${MIRACLE_TIMEOUT}s"
  echo " cells: $*"
  echo "=============================================================="
} | tee -a "$LOG"
for cell in "$@"; do
  echo "=== cell ($cell)  started $(date +%H:%M:%S)" | tee -a "$LOG"
  start=$(date +%s)
  python3 -u miracle.py $cell 2>&1 | tee -a "$LOG"
  echo "    cell ($cell) took $(( $(date +%s) - start ))s" | tee -a "$LOG"
done
echo " resume finished $(date)" | tee -a "$LOG"

#!/bin/bash
# sweep.sh — run miracle.py over every reasonable (j,k,ncuts) and flag any
# MIRACLE, i.e. a candidate new Jacobian counterexample of covering degree
# C(j+k, j).
#
#   chmod +x sweep.sh && ./sweep.sh                     # j+k <= 6
#   MIRACLE_TIMEOUT=120 nohup ./sweep.sh > sweep.log 2>&1 &
#
# Environment knobs (passed through to miracle.py):
#   MIRACLE_TIMEOUT    seconds per COMBINATION (default 120)
#   MIRACLE_MODULUS    Groebner modulus, 32003 by default
#   MIRACLE_SLOW_WARN  print a line for combinations slower than this
#   CELL_CAP           seconds per (j,k,ncuts) cell (default 3600)
#
# Ranges. The map is Sym^j x Sym^k -> Sym^{j+k} and multiplication is
# commutative, so only j <= k is scanned. {Res = 1} has dimension j+k+1 and
# each cut drops one, so ncuts runs 0 .. j+k-2, stopping at C^3 (the smallest
# dimension a counterexample can live in). ncuts = 0 is a real case: it would
# give a counterexample in C^{j+k+1} with no cut at all, and it is where the
# (1,1) SL_2(C) obstruction appears.
#
# WHY THIS DRIVES miracle.py AND NOT symmult.py: symmult.py's scan accumulated
# conflicting patches, including an argument-order bug that passed the TIMEOUT
# as the Groebner MODULUS -- so every call silently ran mod 60 and threw
# SymmetricModularIntegerMod60 errors that masqueraded as real failures. Its
# setup/partition helpers are sound and miracle.py imports them; its scan is
# not. Do not resurrect it.
#
# The per-COMBINATION timeout matters. Earlier sweeps used a per-RUN timeout,
# so one pathological pair (e.g. (3,2)+(2,1,1,1)) consumed the whole budget
# and the entire cell was reported TIMEOUT even though most of it had finished
# cleanly. miracle.py caps each combination and reports UNRESOLVED only for
# those that genuinely did not finish.

set -u
echo "sweep.sh starting (pid $$) at $(date)"      # first line, always
if ! command -v python3 >/dev/null; then echo "FATAL: no python3"; exit 1; fi
if [ ! -f miracle.py ]; then echo "FATAL: miracle.py not in $(pwd)"; exit 1; fi
if ! command -v Singular >/dev/null; then
  echo "WARNING: Singular not found -- falling back to sympy, much slower"
fi
DMAX=${1:-6}
export MIRACLE_TIMEOUT=${MIRACLE_TIMEOUT:-120}
export MIRACLE_SLOW_WARN=${MIRACLE_SLOW_WARN:-5}
CELL_CAP=${CELL_CAP:-3600}
OUT=sweep_results.txt
HITS=sweep_miracles.txt
: > "$OUT"; : > "$HITS"
log() { echo "$@" | tee -a "$OUT"; }

log "=============================================================="
log " miracle sweep   j+k <= $DMAX"
log " per-combination timeout ${MIRACLE_TIMEOUT}s, per-cell cap ${CELL_CAP}s"
log " modulus ${MIRACLE_MODULUS:-32003}  (mod-p NEGATIVES are trustworthy;"
log "   a HIT must be re-run with --exact before it means anything)"
log " started $(date)"
log "=============================================================="

# ---- POSITIVE CONTROL ---------------------------------------------------
# (1,2) ncuts=1 must find exactly one miracle, at partition (2,1). If it does
# not, the test is broken and every negative below is worthless.
log ""
log "POSITIVE CONTROL: (1,2) ncuts=1 must find exactly 1 miracle"
ctl=$(timeout "$CELL_CAP" python3 miracle.py 1 2 1 2>&1)
echo "$ctl" | grep -E 'MIRACLE|miracle\(s\)' | tee -a "$OUT"
ctl_n=$(echo "$ctl" | grep -oE '[0-9]+ miracle\(s\)' | grep -oE '^[0-9]+' | head -1)
if [ "${ctl_n:-0}" != "1" ]; then
  log ""
  log "!!! CONTROL FAILED (expected 1, got ${ctl_n:-none}). Aborting:"
  log "!!! negatives from a broken test are worse than no data."
  exit 1
fi
log "control PASSED"

# ---- SWEEP --------------------------------------------------------------
total=0; new=0; unres=0; capped=0
for (( j=1; j<=DMAX/2; j++ )); do
  for (( k=j; j+k<=DMAX; k++ )); do
    for (( n=0; n<=j+k-2; n++ )); do
      total=$((total+1))
      label="(j,k)=($j,$k) ncuts=$n"
      start=$(date +%s)
      # Stream the cell's output instead of capturing it with $(...).
      # Command substitution buffers everything until the cell FINISHES, so a
      # long cell looks hung for up to CELL_CAP seconds with nothing in the
      # log. -u makes Python unbuffered; tee gives both live output and a
      # copy to parse.
      tmpf=$(mktemp)
      log "  ....      $label  running (live progress below)"
      # tee -a FILE1 FILE2 writes to BOTH files and still passes the stream
      # through to stdout. Writing `tee -a "$OUT" > "$tmpf"` instead sends
      # tee's stdout into the temp file, so nothing ever reaches the log --
      # which is exactly how a working sweep comes to look like a dead one.
      timeout "$CELL_CAP" python3 -u miracle.py "$j" "$k" "$n" 2>&1 \
        | tee -a "$OUT" "$tmpf"
      rc=${PIPESTATUS[0]}
      res=$(cat "$tmpf"); rm -f "$tmpf"
      dur=$(( $(date +%s) - start ))
      if [ $rc -eq 124 ]; then
        capped=$((capped+1))
        partial=$(grep -c 'UNRESOLVED' "$tmpf" 2>/dev/null || echo 0)
        unres=$((unres+partial))
        rm -f "$tmpf"
        log "  CELL-CAP  $label after ${dur}s -- INCOMPLETE, not a negative."
        log "            raise CELL_CAP and re-run: ./resume.sh \"$j $k $n\""
        continue
      fi
      cov=$(echo "$res" | grep -oE 'covering degree [0-9]+' | head -1)
      summ=$(echo "$res" | grep -oE '[0-9]+ miracle\(s\), [0-9]+ unresolved.*' | head -1)
      m=$(echo "$summ" | grep -oE '^[0-9]+')
      u=$(echo "$summ" | grep -oE '[0-9]+ unresolved' | grep -oE '^[0-9]+')
      m=${m:-0}; u=${u:-0}; unres=$((unres+u))
      if [ "$m" -gt 0 ] && [ "$j" = "1" ] && [ "$k" = "2" ] && [ "$n" = "1" ]; then
        log "  KNOWN     $label  $cov  [${dur}s]  <- Tao's counterexample"
      elif [ "$m" -gt 0 ]; then
        new=$((new+m))
        log ""
        log "  *** NEW MIRACLE ***  $label  $cov  [${dur}s]"
        echo "$res" | grep -E 'MIRACLE' | tee -a "$OUT" | tee -a "$HITS"
        { echo "=== $label  $cov"; echo "$res"; echo; } >> "$HITS"
        log ""
      else
        extra=""
        [ "$u" -gt 0 ] && extra="  ($u unresolved)"
        log "  none      $label  $cov  [${dur}s]$extra"
      fi
    done
  done
done

log ""
log "=============================================================="
log " cells swept        : $total"
log " NEW miracles       : $new   (the known (1,2) hit is excluded)"
log " unresolved combos  : $unres   <- NOT negatives; re-run with a larger"
log "                       MIRACLE_TIMEOUT before drawing conclusions"
log " cells hitting cap  : $capped"
log " finished $(date)"
log "=============================================================="
if [ "$new" -gt 0 ]; then
  log ""
  log " A MIRACLE IS ONLY A CANDIDATE. Before believing it:"
  log "  1. re-run that cell with --exact (mod-p positives are NOT proof;"
  log "     mod-p negatives are);"
  log "  2. this test examines ONE degenerate locus (leading coeff of L) --"
  log "     necessary, not sufficient;"
  log "  3. build the explicit map, verify det J is a nonzero constant;"
  log "  4. compute a generic fiber and confirm the covering degree;"
  log "  5. show the WHOLE variety is affine space, not just that fiber."
  log " Details in $HITS"
fi

#!/bin/bash
# checkall.sh — clean-room verification runbook. Run from a FRESH directory
# containing only the final files (see FILES below), on the machine that will
# be cited in the paper.
#
#   mkdir ~/cleanroom && cd ~/cleanroom
#   (copy the 8 files in)
#   chmod +x checkall.sh sweep.sh resume.sh
#   ./checkall.sh            # tiers 0-3, ~15 min
#   ./checkall.sh --sweep    # also tier 4, the full overnight sweep
#
# TIERS, in order. STOP AT THE FIRST FAILURE -- later tiers are meaningless
# after an earlier one fails.
#   0  environment (python3, sympy, Singular)
#   1  positive controls: the machinery must FIND the known counterexample
#   2  theorem verifications: moduli.py spot checks + multicut.py (internal
#      92-config criterion-vs-Groebner check, unit demo)
#   3  fast negative cells: every covering degree <= 10, minutes total
#   4  (--sweep) the full sweep, j+k <= 6, hours; gaps expected ONLY at
#      covering degrees 15 (2,4) and 20 -- and (3,3) resolved fully on the
#      reference machine, so unresolved lines outside (2,4) are a regression
#
# FILES required in this directory:
#   symmult.py miracle.py reconstruct.py moduli.py multicut.py
#   sweep.sh resume.sh checkall.sh
# NOTE: symmult.py is used ONLY for its helpers (form, resultant, apply_D,
# PARTITIONS). Do not run its own scan; it retains known-broken patches.

set -u
LOG=checkall.log
: > "$LOG"
say()  { echo "$@" | tee -a "$LOG"; }
fail() { say ""; say "*** FAILED: $* -- do not proceed to later tiers."; exit 1; }

say "clean-room verification  $(date)"
say "host: $(uname -a | cut -c1-70)"
say ""

# ---- TIER 0: environment ------------------------------------------------
say "TIER 0: environment"
command -v python3 >/dev/null || fail "python3 missing"
python3 -c 'import sympy' 2>/dev/null || fail "sympy missing (pip install sympy)"
command -v Singular >/dev/null || fail "Singular missing (apt install singular)"
# NOTE: never call `Singular --version` -- it enters the interactive REPL
# and hangs the script. Query the version through a piped quit instead.
SINGVER=$(echo 'exit;' | timeout 10 Singular -q 2>&1 | head -1 | cut -c1-40)
say "  python3 $(python3 -V 2>&1 | cut -d' ' -f2), sympy $(python3 -c 'import sympy;print(sympy.__version__)'), Singular ${SINGVER:-present}"
for f in symmult.py miracle.py reconstruct.py moduli.py multicut.py; do
  [ -f "$f" ] || fail "$f missing"
done
python3 - << 'EOF' || exit 1
import symmult
assert sorted(symmult.PARTITIONS) == [3,4,5,6], 'PARTITIONS wrong'
print('  symmult helpers import, PARTITIONS degrees 3-6')
EOF
say "  TIER 0 PASSED"
say ""

# ---- TIER 1: positive controls -----------------------------------------
say "TIER 1: positive controls (the machinery must find Alpoge's counterexample)"
out=$(timeout 600 python3 miracle.py 1 2 1 2>&1)
echo "$out" >> "$LOG"
echo "$out" | grep -q '1 miracle(s), 0 unresolved' || fail "miracle.py control: expected exactly 1 miracle"
echo "$out" | grep -q 'MIRACLE' || fail "miracle.py control: no MIRACLE line"
echo "$out" | grep -q '>>> EXPLICIT COUNTEREXAMPLE' || fail "reconstruction did not produce an explicit counterexample"
echo "$out" | grep -qE 'det J = -?1/[0-9]+ +constant and nonzero: True' || fail "det J not a nonzero constant"
echo "$out" | grep -q 'a tested fiber has 3 distinct points' \
  || fail "tested fiber does not have 3 points"
say "  miracle.py 1 2 1: exactly 1 MIRACLE at (2,1), auto-reconstructed:"
say "    det J constant, a tested fiber has 3 distinct points (noninjective)"
say "    >>> EXPLICIT COUNTEREXAMPLE"
say "  TIER 1 PASSED"
say ""

# ---- TIER 2: theorem verification (Euler characteristic) ---------------
# REVISED after review. The paper's j=1 theorem rests on the GLOBAL Euler-
# characteristic argument, so Tier 2 now runs verify_chi_tables.py, which
# reproduces the 20/108/110 tables with assertions. multicut.py has been
# demoted to a diagnostic and is run separately, labeled as such: it tests
# the DEGENERATE FIBER only and cannot certify global non-affineness.
say "TIER 2: theorem verification -- Euler characteristic tables"
# --full so the archived log contains all 238 configurations, making the
# aggregate assertions auditable line by line rather than on trust.
chi=$(timeout 3600 python3 verify_chi_tables.py --full 2>&1)
echo "$chi" >> "$LOG"
echo "$chi" | grep -q 'ALL ASSERTIONS PASSED' || fail "verify_chi_tables.py"
for line in "count == 20" "count == 108" "count == 110"; do
  echo "$chi" | grep -q "$line" || fail "missing configuration count: $line"
done
echo "$chi" | grep -c 'chi(V) = 1 occurs nowhere    : OK' | grep -q '^3$' \
  || fail "chi(V)=1 exclusion not asserted for all three k"
say "  verify_chi_tables.py: 20/108/110 independent configurations,"
say "    all resolved, chi values {0,k-1}, chi(V)=1 nowhere. ASSERTIONS PASSED"
m1=$(timeout 300 python3 moduli.py 2 "2,1" 2>&1 | tail -1)
say "  moduli control (1,2)(2,1):   $m1"
echo "$m1" | grep -q 'MIRACLE for all moduli' || fail "moduli positive control"
m2=$(timeout 300 python3 moduli.py 2 "1,1,1" 2>&1 | tail -1)
say "  moduli control (1,2)(1,1,1): $m2"
echo "$m2" | grep -q 'NEVER a miracle' || fail "moduli negative control"
eu=$(timeout 300 python3 euler.py 2>&1)
echo "$eu" >> "$LOG"
echo "$eu" | grep -q 'V is NOT C^3' || fail "euler.py (1,1,1) result"
say "  euler.py: (1,1,1) chi(V)=0 != 1; (2,1) control chi(V)=1"
au=$(timeout 300 python3 audit.py 2>&1)
echo "$au" >> "$LOG"
echo "$au" | grep -q 'UNCONDITIONAL' || fail "audit.py units certificate"
say "  audit.py: global-unit certificates for the full-derivative classes"
dg=$(timeout 900 python3 multicut.py 2>&1)
echo "$dg" >> "$LOG"
echo "$dg" | grep -q '92 configurations: ALL AGREE' || fail "multicut diagnostic"
say "  multicut.py [DIAGNOSTIC ONLY -- degenerate fiber, not global]:"
say "    92/92 left-kernel vs Groebner agreement"
say "  TIER 2 PASSED"
say ""

# ---- TIER 3: fast negative cells (covering degrees <= 10) ---------------
say "TIER 3: EXPLORATORY SCREENING, covering degrees 2-10 (expect 0 candidates)"
say "  NOTE: the triangular-basis criterion is a SUFFICIENT certificate for"
say "  affineness, not an iff test, so a null result here is screening --"
say "  NOT rigorous exclusion. The j=1 theorem does not rely on it."
for cell in "1 1 0" "1 3 0" "1 3 1" "1 3 2" "1 4 0" "1 4 1" "1 4 2" "1 4 3" \
            "1 5 0" "1 5 1" "1 5 2" "1 5 3" "1 5 4" \
            "2 2 0" "2 2 1" "2 2 2" "2 3 0" "2 3 1" "2 3 2" "2 3 3"; do
  out=$(MIRACLE_SLOW_WARN=1000 timeout 3600 python3 -u miracle.py $cell 2>&1)
  echo "$out" >> "$LOG"
  summ=$(echo "$out" | grep -oE '[0-9]+ miracle\(s\), [0-9]+ unresolved.*' | head -1)
  say "  ($cell): ${summ:-NO SUMMARY}"
  m=$(echo "$summ" | grep -oE '^[0-9]+'); u=$(echo "$summ" | grep -oE '[0-9]+ unresolved' | grep -oE '^[0-9]+')
  [ "${m:-1}" = "0" ] || fail "cell ($cell) reports a miracle -- investigate before anything else"
  [ "${u:-1}" = "0" ] || fail "cell ($cell) has unresolved combos -- raise MIRACLE_TIMEOUT and re-run"
done
say "  TIER 3 PASSED: no candidates surfaced at covering degrees <= 10"
say ""

# ---- TIER 4: the full sweep (optional, overnight) -----------------------
if [ "${1:-}" = "--sweep" ]; then
  say "TIER 4: full sweep j+k <= 6 (hours; log streams to sweep.log)"
  MIRACLE_TIMEOUT=${MIRACLE_TIMEOUT:-600} MIRACLE_MEM_MB=${MIRACLE_MEM_MB:-4000} ./sweep.sh 6
  say "  compare sweep_results.txt to the reference: 0 NEW miracles;"
  say "  unresolved lines acceptable ONLY in (2,4) cells (covering degree 15)."
else
  say "TIER 4 skipped (run ./checkall.sh --sweep for the full overnight sweep)"
fi
say ""
say "clean-room verification complete $(date) -- log in $LOG"

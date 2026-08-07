"""
verify_chi_tables.py — the clean-room verification for the j = 1 theorem.

THIS is the script that reproduces the computational claims of the paper.
It replaces the earlier reliance on multicut.py and the Groebner sweep, both
of which test the DEGENERATE FIBER only and therefore cannot certify global
non-affineness (see the note in multicut.py).

What it does, deterministically and with assertions:
  1. builds the canonical cut classes for each k (partitions of k+1);
  2. keeps only LINEARLY INDEPENDENT cut collections -- required by Lemma 2
     of the paper, which is false for dependent cuts;
  3. runs both value patterns (1,..,1) and (1,0,..,0);
  4. computes chi(V) with chi_plane.py (exact over Q, no modulus, no
     hypothesis (H), no fiber bridge);
  5. asserts the configuration counts, that every case RESOLVED, that the
     value set is exactly {0, k-1}, and that chi(V) = 1 occurs NOWHERE.

ENUMERATION SCOPE (stated exactly, because the totals are otherwise
mysterious). Dimension n = k + 2 - m >= 3 forces m <= k - 1.
    k = 3: m = 2 only.  C(5,2) = 10 pairs x 2 value patterns   =  20
    k = 4: m = 2 and 3. (C(7,2) + C(7,3)) x 2 = 112, minus 4 cases from
           two linearly DEPENDENT cut triples                  = 108
    k = 5: m = 2 only (pairs), C(11,2) = 55 x 2                = 110
The k = 5 row is pairs only; m = 3, 4 at k = 5 are NOT enumerated here and
the paper must not claim them. Everything reported is exactly what is run.

Usage:  python3 verify_chi_tables.py            (summary + assertions)
        python3 verify_chi_tables.py --full     (also prints all 238
                                                 configurations, for the
                                                 archived provenance log)
Exit code 0 iff every assertion passes.
"""
import sys
from itertools import combinations

import sympy as sp

from symmult import apply_D, PARTITIONS
from chi_plane import chi_plane, controls

z, w = sp.symbols('z w')

# (k, list of m values enumerated, expected configuration count)
SCOPE = [(3, [2], 20), (4, [2, 3], 108), (5, [2], 110)]


def functional_vector(gammas, d):
    """Coefficient vector of the cut functional on Sym^d, for independence
    testing. Two cuts are dependent iff these vectors are."""
    Cs = sp.symbols(f'C0:{d+1}')
    G = sum(Cs[i] * z**(d - i) * w**i for i in range(d + 1))
    DG = sp.expand(apply_D(gammas, G))
    return [sp.diff(DG, C) for C in Cs]


def independent(gam_list, d):
    M = sp.Matrix([functional_vector(g, d) for g in gam_list])
    return M.rank() == len(gam_list)


def run(full=False):
    print('verify_chi_tables.py — Euler-characteristic tables for j = 1')
    print('=' * 66)
    print('CONTROLS first (chi_plane against euler.py and the known miracle):')
    if not controls():
        print('CONTROLS FAILED — every number below would be worthless.')
        return 1
    print()

    failures = []
    for k, ms, expected in SCOPE:
        labs = list(PARTITIONS[k + 1])
        seen, dropped = [], 0
        for m in ms:
            assert m <= k - 1, 'scope error: n = k+2-m >= 3 requires m <= k-1'
            for combo in combinations(labs, m):
                gam = [PARTITIONS[k + 1][c] for c in combo]
                if not independent(gam, k + 1):
                    dropped += 1
                    continue
                for vals in [(1,) * m, (1,) + (0,) * (m - 1)]:
                    seen.append((combo, vals))
        results = {}
        unresolved = []
        rows_out = []
        for combo, vals in seen:
            v = chi_plane(k, list(combo), list(vals), verbose=False)
            if v is None:
                unresolved.append((combo, vals))
            else:
                results[(combo, vals)] = v
            rows_out.append((combo, vals, v))
        if full:
            # Print EVERY configuration so the aggregate assertions below can
            # be spot-checked rather than taken on trust. Without this a
            # reviewer must rebuild the enumeration harness to audit any
            # single entry -- which is precisely what happened.
            print(f'   --- all {len(rows_out)} configurations at k = {k} ---')
            for combo, vals, v in rows_out:
                tag = 'UNRESOLVED' if v is None else f'chi(V) = {v}'
                flag = '   <== chi = 1 !!' if v == 1 else ''
                print(f'   {"+".join(combo):<42} {str(tuple(vals)):<12} '
                      f'{tag}{flag}')
            print()
        vals_set = sorted(set(results.values()))
        n_cfg = len(seen)
        ok_count = (n_cfg == expected)
        ok_resolved = (not unresolved)
        ok_set = (vals_set == sorted({0, k - 1}))
        ok_no_one = (1 not in results.values())
        print(f'k = {k}:  m in {ms},  {n_cfg} independent configurations '
              f'({dropped} dependent cut sets dropped)')
        print(f'   count == {expected}          : '
              f'{"OK" if ok_count else f"FAIL ({n_cfg})"}')
        print(f'   all resolved                 : '
              f'{"OK" if ok_resolved else f"FAIL ({len(unresolved)})"}')
        print(f'   chi values == {{0, {k-1}}}        : '
              f'{"OK" if ok_set else f"FAIL {vals_set}"}')
        print(f'   chi(V) = 1 occurs nowhere    : '
              f'{"OK" if ok_no_one else "FAIL"}')
        dist = {v: list(results.values()).count(v) for v in vals_set}
        print(f'   distribution                 : {dist}')
        for name, ok in [('count', ok_count), ('resolved', ok_resolved),
                         ('value set', ok_set), ('no chi=1', ok_no_one)]:
            if not ok:
                failures.append(f'k={k} {name}')
        print()

    print('=' * 66)
    if failures:
        print('FAILURES:', ', '.join(failures))
        return 1
    print('ALL ASSERTIONS PASSED.')
    print('Every enumerated configuration with n >= 3 and independent cuts')
    print('has chi(V) in {0, k-1}; chi(V) = 1 occurs only at the control')
    print('(1,2) class (2,1) with nonzero cut value — Alpoge\'s example.')
    return 0


if __name__ == '__main__':
    # Accept -full, --full and full: a flag that silently does nothing when
    # the user is one dash off is a trap, and an unrecognized argument should
    # say so rather than being ignored.
    args = [a.lstrip('-').lower() for a in sys.argv[1:]]
    unknown = [a for a in args if a != 'full']
    if unknown:
        print(f'unrecognized argument(s): {unknown}. '
              f'Usage: verify_chi_tables.py [--full]')
        sys.exit(2)
    sys.exit(run(full=('full' in args)))

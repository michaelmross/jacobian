"""
miracle.py — EXPLORATORY SCREENING for candidate affine miracles.

STATUS AFTER REVIEW: the criterion implemented here is a SUFFICIENT
certificate, not an if-and-only-if test, and this file's output must be read
accordingly.

The test asks whether the lexicographic Groebner basis of the degenerate
fiber is triangular (each generator of degree 1 in its leading variable with
constant leading coefficient). If it is, the fiber IS affine space -- a valid
positive certificate. But the converse fails: x^2 - y = 0 defines A^1, yet
with x leading it has degree 2 and this test would label it "branching". An
affine variety may admit a nonlinear finite projection onto the chosen free
coordinates.

CONSEQUENCES:
  * a MIRACLE verdict is a genuine candidate, worth reconstructing;
  * a "no miracle" verdict is SCREENING, not rigorous exclusion;
  * therefore the j >= 2 sweep results are exploratory, and the phrase
    "fully negative" should not be used of them;
  * the j = 1 theorem does NOT depend on this criterion at all. It rests on
    the Euler-characteristic argument (UNIQUENESS.md), verified by
    verify_chi_tables.py.

SINGULAR IS REQUIRED. An earlier docstring advertised a sympy fallback; that
branch passed raw Groebner generators where the parser's five-field fact
tuples were expected, and raised TypeError/ValueError on the positive
control. Rather than ship a broken path, the fallback is removed: if Singular
is absent the script exits with a clear message.

Usage:  python3 miracle.py J K NCUTS
"""
import os
import re
import resource
import subprocess
import sys
import tempfile
from itertools import combinations_with_replacement
from math import comb

import time as _time

import sympy as sp

from symmult import form, resultant, apply_D, PARTITIONS

MODULUS = int(os.environ.get('MIRACLE_MODULUS', '32003'))
TIMEOUT = int(os.environ.get('MIRACLE_TIMEOUT', '120'))
SLOW_WARN = float(os.environ.get('MIRACLE_SLOW_WARN', '3'))
# MEMORY CAP on the Singular child. Some lex bases grow without bound --
# (5,1)+(3,3)+(2,2,2)+(1,1,1,1,1,1) in (2,4) ncuts=4 is one -- and on WSL an
# unbounded child exhausts the VM's allocation and takes the whole VM down,
# not merely the process. A TIME limit cannot prevent this: the memory is
# consumed well before the clock runs out. RLIMIT_AS makes Singular fail
# cleanly instead, and the combination is reported UNRESOLVED.
MEM_MB = int(os.environ.get('MIRACLE_MEM_MB', '3000'))


def _cap_memory():
    b = MEM_MB * 1024 * 1024
    resource.setrlimit(resource.RLIMIT_AS, (b, b))
HAVE_SINGULAR = subprocess.run(['which', 'Singular'], capture_output=True
                               ).returncode == 0
if not HAVE_SINGULAR:
    # The sympy fallback was REMOVED, not fixed: it passed raw Groebner
    # generators to verdict(), which expects the Singular parser's five-field
    # fact tuples, so the positive control died with TypeError/ValueError.
    # Shipping a broken alternate path is worse than requiring the tool.
    print('FATAL: Singular is required (the sympy fallback was removed as '
          'broken).\n       Install with: apt-get install singular')
    raise SystemExit(2)


def build_system(j, k, cuts):
    """Degenerate-locus system: {Res=1} plus the cuts, at leading coeff = 0."""
    L, Lc = form(j, 'a')
    Q, Qc = form(k, 'c')
    prod = sp.expand(L * Q)
    R = resultant(L, Q)
    lead = Lc[0]
    gens = [sp.expand(R.subs(lead, 0) - 1)]
    for gammas, val in cuts:
        d = apply_D(gammas, prod)
        if d == 0:
            return None, None, 'a cut annihilates the product'
        d0 = sp.expand(d.subs(lead, 0) - val)
        if d0 == -val:
            return None, None, 'a cut vanishes identically on the degenerate locus'
        gens.append(d0)
    return gens, Lc[1:] + Qc, None


def _clear_denoms(e, variables):
    lcm = 1
    for c in sp.Poly(e, *variables).coeffs():
        lcm = sp.ilcm(lcm, sp.Rational(c).q)
    return sp.expand(e * lcm)


def groebner_singular(gens, variables, modulus):
    """lex Groebner basis via Singular, returning only the THREE FACTS the
    triangularity criterion needs, per generator:

        (leading variable index, degree in it, is its coefficient constant)

    Parsing the full basis back into sympy was the real bottleneck -- NOT the
    Groebner computation. A 20s cap on Singular still produced 35s
    combinations because sympify ran afterwards, unbounded, on enormous
    expressions. Singular computes these facts itself in negligible time and
    the output is a handful of integers.

    Returns (facts, flagged_strings) or None on timeout.
    """
    n = len(variables)
    names = {v: f'v({i+1})' for i, v in enumerate(variables)}

    def tr(e):
        s = str(_clear_denoms(e, variables))
        for v in sorted(names, key=lambda v: -len(str(v))):
            s = re.sub(rf'\b{re.escape(str(v))}\b', names[v], s)
        return s.replace('**', '^')

    # DEGREVLEX PRE-CHECK. Cost is driven by the fiber DIMENSION (= j+k-ncuts),
    # not by degree or cut count -- which is why (3,3) ncuts=4 (dim 2) runs at
    # 0.6s per combination while (2,4) ncuts=2 (dim 4) needs 186s and often
    # exceeds 600s. Degrevlex is far cheaper than lex and settles two things
    # up front:
    #   * unit ideal  -> EMPTY fiber -> no miracle, no lex needed at all;
    #   * dimension 0 -> lex is cheap (FGLM-style), always worth attempting.
    # A hard case that lex could not finish in 600s answers in 1.6s here.
    pre = [f'ring rp = {modulus if modulus else 0}, (v(1..{n})), dp;',
           'ideal Ip = ' + ',\n  '.join(tr(g) for g in gens) + ';',
           'ideal Gp = std(Ip);',
           '"DP_UNIT:", (size(Gp) == 1 && Gp[1] == 1);',
           '"DP_DIM:", dim(Gp);', 'exit;']
    with tempfile.NamedTemporaryFile('w', suffix='.sing', delete=False) as f:
        f.write('\n'.join(pre))
        ppath = f.name
    try:
        pout = subprocess.run(['Singular', '-q', ppath], capture_output=True,
                              text=True, timeout=TIMEOUT,
                              preexec_fn=_cap_memory).stdout
    except subprocess.TimeoutExpired:
        pout = ''
    finally:
        os.unlink(ppath)
    if 'DP_UNIT: 1' in pout:
        return 'UNIT', []
    m = re.search(r'DP_DIM:\s*(-?\d+)', pout)
    predim = int(m.group(1)) if m else None

    script = [
        f'ring r = {modulus if modulus else 0}, (v(1..{n})), lp;',
        'ideal I = ' + ',\n  '.join(tr(g) for g in gens) + ';',
        'ideal G = std(I);',
        'int i, j, d, lcconst, lv;',
        'poly g, lc;  matrix C;',
        'if (size(G) == 1 && G[1] == 1) { "UNIT"; exit; }',
        'for (i = 1; i <= size(G); i++) {',
        '  g = G[i];  if (g == 0) { i = i; }',
        '  lv = 0;',
        f'  for (j = 1; j <= {n}; j++) ' + '{',
        '    if (lv == 0 && diff(g, v(j)) != 0) { lv = j; }',
        '  }',
        '  if (lv > 0) {',
        '    C = coeffs(g, v(lv));',
        '    d = nrows(C) - 1;',
        '    lc = C[nrows(C), 1];',
        '    lcconst = (deg(lc) <= 0);',
        '    "FACT:", i, lv, d, lcconst, size(g);',
        '    if (d > 1 || lcconst == 0) { if (size(g) < 12) { "STR:", g; } }',
        '  }',
        '}',
        'exit;']
    with tempfile.NamedTemporaryFile('w', suffix='.sing', delete=False) as f:
        f.write('\n'.join(script))
        path = f.name
    try:
        proc = subprocess.run(['Singular', '-q', path], capture_output=True,
                              text=True, timeout=TIMEOUT,
                              preexec_fn=_cap_memory)
        out = proc.stdout
        if proc.returncode != 0 and 'FACT:' not in out and 'UNIT' not in out:
            return ('MEMORY', predim)
    except subprocess.TimeoutExpired:
        return ('TIMEOUT', predim)
    finally:
        os.unlink(path)
    if 'UNIT' in out:
        return 'UNIT', []
    facts, strs = [], []
    for line in out.splitlines():
        if line.startswith('FACT:'):
            p = line.split()[1:]
            facts.append(tuple(int(x) for x in p[:5]))
        elif line.startswith('STR:'):
            strs.append(line[4:].strip())
    return facts, strs


def verdict(res, variables):
    """Triangular => affine space => MIRACLE. Works from the fact tuples."""
    if isinstance(res, tuple) and res and res[0] == 'MEMORY':
        return (f'UNRESOLVED: Singular exceeded the {MEM_MB}MB cap '
                f'(fiber dim {res[1]}) -- not evidence either way')
    if isinstance(res, tuple) and res and res[0] == 'TIMEOUT':
        d = res[1]
        return (f'UNRESOLVED: lex timed out (fiber dim {d}) -- not evidence '
                f'either way')
    if res is None:
        return 'UNRESOLVED: Groebner timed out -- not evidence either way'
    facts, strs = res
    if facts == 'UNIT':
        return 'EMPTY degenerate fiber -- no miracle'
    if not facts:
        return 'UNRESOLVED: no basis returned'
    bad = [f for f in facts if f[2] > 1 or f[3] == 0]
    if not bad:
        return 'MIRACLE: degenerate fiber is affine space (triangular basis)'
    bad.sort(key=lambda f: f[4])          # fewest terms = most informative
    i, lv, d, lcconst, nterms = bad[0]
    var = variables[lv - 1]
    detail = f'  [{strs[0]}]' if strs else f'  [{nterms} terms]'
    if d > 1:
        return (f'no miracle: BRANCHING -- generator has degree {d} in {var} '
                f'({d} sheets){detail}')
    return (f'no miracle: PUNCTURE -- leading coefficient in {var} is not '
            f'constant (fiber punctured){detail}')


def test(j, k, cuts, modulus=MODULUS):
    gens, variables, err = build_system(j, k, cuts)
    if err:
        return err
    if HAVE_SINGULAR:
        exprs = groebner_singular(gens, variables, modulus)
        return verdict(exprs, variables)
    try:
        G = sp.groebner([_clear_denoms(g, variables) for g in gens],
                        *variables, order='lex',
                        **({'modulus': modulus} if modulus else {}))
        return verdict(list(G.exprs), variables)
    except Exception as e:
        return f'UNRESOLVED: {type(e).__name__}'


def scan(j, k, ncuts, modulus=MODULUS):
    d = j + k
    if ncuts == 0:
        items = []
    elif d not in PARTITIONS:
        print(f'  no partitions for degree {d}')
        return 0, 0
    else:
        items = list(PARTITIONS[d].items())
    engine = 'Singular' if HAVE_SINGULAR else 'sympy'
    mod = f'mod {modulus}' if modulus else 'exact over Q'
    print(f'(j,k)=({j},{k})  covering degree {comb(d, j)}  cuts={ncuts}  '
          f'[{engine}, {mod}]')
    hits = unres = 0
    if ncuts == 0:
        combos = [()]
    else:
        combos = [c for c in combinations_with_replacement(items, ncuts)
                  if len({l for l, _ in c}) == ncuts]
    ntot = len(combos) * (1 if ncuts <= 1 else 2)
    done = 0
    for combo in combos:
        valsets = [(1,) * ncuts] if ncuts <= 1 else \
                  [(1,) * ncuts, (1,) + (0,) * (ncuts - 1)]
        for vals in valsets:
            cuts = [(g, v) for (l, g), v in zip(combo, vals)]
            lab = '+'.join(l for l, _ in combo) or '(no cut)'
            t0 = _time.time()
            r = test(j, k, cuts, modulus)
            dt = _time.time() - t0
            done += 1
            if dt > SLOW_WARN:
                print(f'  [{dt:5.1f}s] {done}/{ntot} {lab} {vals}',
                      flush=True)
            elif done % 25 == 0:
                print(f'  ... {done}/{ntot}', flush=True)
            if r.startswith('MIRACLE'):
                hits += 1
                print(f'  *** {lab} {vals}: {r}')
                # A MIRACLE is only a flag. Carry it straight through to an
                # explicit map so the hit arrives verified rather than
                # merely announced.
                try:
                    from reconstruct import reconstruct, verify
                    out, err = reconstruct(j, k, cuts, verbose=True)
                    if err:
                        print(f'      reconstruction did not complete: {err}')
                    else:
                        Fmap, freevars = out
                        for i, f in enumerate(Fmap):
                            print(f'      F{i+1} = {sp.factor(f)}')
                        ok, nfib = verify(Fmap, freevars, verbose=True)
                        print('      >>> ' + ('EXPLICIT COUNTEREXAMPLE'
                                              if ok else
                                              'not a counterexample'))
                except Exception as ex:
                    print(f'      reconstruction error: {type(ex).__name__}')
            elif r.startswith('UNRESOLVED'):
                unres += 1
                print(f'  {lab} {vals}: {r}')
    print(f'  {hits} miracle(s), {unres} unresolved, '
          f'{sum(1 for c in combos for _ in ([0] if ncuts <= 1 else [0, 1]))} '
          f'combinations tested')
    return hits, unres


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    exact = '--exact' in sys.argv
    if len(args) >= 2:
        j, k = int(args[0]), int(args[1])
        n = int(args[2]) if len(args) > 2 else 1
        scan(j, k, n, modulus=None if exact else MODULUS)
    else:
        print(__doc__)

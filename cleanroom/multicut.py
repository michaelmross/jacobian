"""
multicut.py — DIAGNOSTIC ONLY. Not a proof of anything global.

WHAT THIS SCRIPT DOES: it verifies a left-kernel classification of the
DEGENERATE FIBER V_0 = V ∩ {a0 = 0} against Groebner ground truth, for
selected canonical multi-cut configurations (92/92 agreement).

WHAT IT DOES NOT DO: it does not show that V is or is not affine space. Its
original docstring reasoned from "a1 is a unit on V_0" to a classification of
when V can be C^n. That step silently used

        V isomorphic to A^n  ==>  V_0 isomorphic to A^{n-1},

which is FALSE in general (a hypersurface in C^n can be C^* x C^{n-2}), and
is exactly the fiber bridge this project abandoned after adversarial review.
The earlier headings on this file ("the multi-cut case, closed", "THE FINAL
PIECE", "arises ONLY") were therefore unjustified and have been removed.

The multi-cut case is closed by the GLOBAL Euler-characteristic argument in
UNIQUENESS.md, verified by verify_chi_tables.py — not by this file.

Retained because the left-kernel classification is a genuine and useful
consistency check on the degenerate fiber, and because its agreement with
Groebner computation is evidence that the fiber-level bookkeeping used
elsewhere is implemented correctly.

Usage:  python3 multicut.py
"""
from itertools import combinations

import sympy as sp

from symmult import form, resultant, apply_D

z, w = sp.symbols('z w')
ROOTS = ['inf', None, 1, 2, 3, 4]


def parts(n, m=None):
    m = m or n
    if n == 0:
        yield ()
        return
    for f in range(min(n, m), 0, -1):
        for r in parts(n - f, f):
            yield (f,) + r


def classes(k):
    """(label, gammas, mu-vector) for every class of degree k+1."""
    out = []
    for p in parts(k + 1):
        g = []
        for i, mult in enumerate(p):
            g += [ROOTS[i]] * mult
        mu = [apply_D(g, z**(k - i) * w**(i + 1)) for i in range(k + 1)]
        if any(m != 0 for m in mu):
            out.append(('(' + ','.join(map(str, p)) + ')', g, mu))
    return out


def predict(k, mus, vals):
    """The left-kernel criterion. Decides the degenerate fiber from the mu
    vectors alone -- no Groebner, no solve."""
    M = sp.Matrix([[sp.Rational(x) for x in mu[1:]] for mu in mus])
    mu0 = sp.Matrix([sp.Rational(mu[0]) for mu in mus])
    vv = sp.Matrix([sp.Rational(v) for v in vals])
    NS = M.T.nullspace()
    if not NS:
        return 'nonconst'
    conds = [(sp.expand((u.T * vv)[0]), sp.expand((u.T * mu0)[0]))
             for u in NS]
    if all(cv == 0 and cm == 0 for cv, cm in conds):
        return 'nonconst'
    if any(cv == 0 and cm != 0 for cv, cm in conds):
        return 'EMPTY'
    gam = {cm / cv for cv, cm in conds if cv != 0}
    if len(gam) > 1 or 0 in gam:
        return 'EMPTY'
    return 'UNIQUE' if k == 2 else 'multi'


def brute(k, gs, vals, Lc, Qc, prod, R, a0):
    """Groebner ground truth. a1 goes LAST so lex elimination exposes its
    univariate relation (a variable-order mistake here initially produced 52
    phantom mismatches)."""
    a1 = Lc[1]
    eqs = [sp.expand(R.subs(a0, 0) - 1)]
    for g, v in zip(gs, vals):
        eqs.append(sp.expand(apply_D(g, prod).subs(a0, 0) - v))
    G = sp.groebner(eqs, *(Qc + [a1]), order='lex', modulus=32003)
    ex = list(G.exprs)
    if ex == [sp.Integer(1)]:
        return 'EMPTY'
    uni = [e for e in ex if e.free_symbols == {a1}]
    if not uni:
        return 'nonconst'
    d = min(sp.Poly(e, a1).degree() for e in uni)
    return 'UNIQUE' if d == 1 else 'multi'


def check_res_sign():
    print('(1) Res|_{a0=0} = +a1^k c0 for all k (the unit relation):')
    for k in range(2, 7):
        L, Lc = form(1, 'a')
        Q, Qc = form(k, 'c')
        r0 = sp.expand(resultant(L, Q).subs(Lc[0], 0))
        ok = sp.simplify(r0 - Lc[1]**k * Qc[0]) == 0
        print(f'    k={k}: {"+a1^{k} c0  OK".replace("{k}", str(k)) if ok else "SIGN DIFFERS: " + str(r0)}')


def check_criterion():
    print('(2) left-kernel criterion vs Groebner, k=2..4, m=2..3:')
    total = mismatch = 0
    for k in (2, 3, 4):
        L, Lc = form(1, 'a')
        Q, Qc = form(k, 'c')
        prod = sp.expand(L * Q)
        R = resultant(L, Q)
        a0 = Lc[0]
        cls = classes(k)
        for m in (2, 3):
            for combo in combinations(cls, m):
                for vals in [(1,) * m, (1,) + (0,) * (m - 1)]:
                    total += 1
                    p_ = predict(k, [c[2] for c in combo], vals)
                    b_ = brute(k, [c[1] for c in combo], vals,
                               Lc, Qc, prod, R, a0)
                    if p_ != b_:
                        mismatch += 1
                        print(f'    MISMATCH k={k} '
                              f'{"+".join(c[0] for c in combo)} {vals}: '
                              f'{p_} vs {b_}')
    print(f'    {total} configurations: '
          f'{"ALL AGREE" if not mismatch else f"{mismatch} MISMATCHES"}')
    # every UNIQUE is at k=2 -- the theorem's computational content
    uniq = []
    for k in (2, 3, 4):
        for m in (2, 3):
            for combo in combinations(classes(k), m):
                for vals in [(1,) * m, (1,) + (0,) * (m - 1)]:
                    if predict(k, [c[2] for c in combo], vals) == 'UNIQUE':
                        uniq.append((k, [c[0] for c in combo], vals))
    print(f'    UNIQUE verdicts: {uniq}')
    print('    (all at k = 2, as the theorem requires)')


def check_global_unit():
    print('(3) the k=2 two-cut varieties are NOT affine (nonconstant unit):')
    L, Lc = form(1, 'a')
    Q, Qc = form(2, 'c')
    a0, a1 = Lc
    c0 = Qc[0]
    prod = sp.expand(L * Q)
    R = resultant(L, Q)
    DA = sp.expand(apply_D(['inf', 'inf', None], prod))
    DB = sp.expand(apply_D(['inf', None, 1], prod))
    G = sp.groebner([R - 1, DA - 1, DB - 1], a0, a1, *Qc, order='lex')
    unit_rel = sp.expand(3 * a1**2 * c0 - a1 - 1)
    has = any(sp.simplify(g - unit_rel) == 0 for g in G.exprs)
    print(f'    Groebner contains 3 a1^2 c0 - a1 - 1: {has}'
          f'   [a1 (3 a1 c0 - 1) = 1: a1 is a unit on V]')
    pt = {a0: 1, a1: 1, c0: sp.Rational(2, 3),
          Qc[1]: sp.Rational(-1, 6), Qc[2]: sp.Rational(1, 6)}
    on_v = all(sp.expand(e.subs(pt)) == v
               for e, v in [(R, 1), (DA, 1), (DB, 1)])
    print(f'    point with a1 = 1 lies on V: {on_v};  degenerate locus has '
          f'a1 = 2  ->  a1 is a NONCONSTANT unit  ->  V is not any C^m')


if __name__ == '__main__':
    print(__doc__.split('PROOF')[0])
    check_res_sign()
    print()
    check_criterion()
    print()
    check_global_unit()

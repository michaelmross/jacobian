"""
moduli.py — close the cross-ratio gap SYMBOLICALLY, not by sampling.

The sweep tests one representative per hyperplane class. That is complete only
for classes with at most 3 distinct roots, since PGL_2 is 3-transitive and can
normalize three roots but no more. A class with m distinct roots is really an
(m-3)-dimensional family parametrized by cross-ratios, and the sweep sampled a
single point of it.

Sampling more points would never be conclusive. Instead, make the free roots
SYMBOLIC and analyze the eliminant as a polynomial in them. For j = 1 the
degenerate system at a0 = 0 is

    Res:  +-a1^k * c0 = 1        =>   c0 = +-1 / a1^k
    cut:  D(prod)|_{a0=0} = 1

so substituting c0 and clearing denominators gives a single ELIMINANT in a1
whose coefficients involve c1..ck and the cross-ratio parameters. Then

    MIRACLE  <=>  eliminant has degree 1 in a1  AND  its root is a constant
                  (degree > 1 => BRANCHING; root depending on c => PUNCTURE)

Because the coefficients are polynomials in the parameters, the degree is
determined once and for all, with special parameter values appearing as the
locus where a leading coefficient vanishes identically in c1..ck.

Result for the (1,3) class (1,1,1,1): the eliminant is

    4*a1^3*(c1*lam + c1 - 3*c2*lam/2) + a1^2 - 6

whose a1^2 coefficient is the CONSTANT 1, inherited from the resultant
normalization. It can never vanish, so the eliminant has degree >= 2 in a1 for
EVERY cross-ratio lam -- at least two branches, never a miracle. No sampling
required.

Usage:
  python3 moduli.py            # every class with moduli, j = 1, degrees 4-6
"""
import sys

import sympy as sp

from symmult import form, resultant, apply_D

z, w = sp.symbols('z w')


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for first in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


def gammas_symbolic(part):
    """Roots for a partition, with the free ones SYMBOLIC.

    PGL_2 fixes three roots; take them as infinity, 0, 1. Every further
    distinct root is a modulus and becomes a symbol lam0, lam1, ...
    """
    fixed = ['inf', None, 1]
    out, params = [], []
    for idx, mult in enumerate(part):
        if idx < 3:
            r = fixed[idx]
        else:
            r = sp.Symbol(f'lam{idx-3}')
            params.append(r)
        out += [r] * mult
    return out, params


def analyze(k, part, verbose=True):
    """Analyze the (1,k) single-cut case for a class, over symbolic moduli."""
    L, Lc = form(1, 'a')
    Q, Qc = form(k, 'c')
    prod = sp.expand(L * Q)
    R = resultant(L, Q)
    lead, b = Lc[0], Lc[1]
    gam, params = gammas_symbolic(part)
    D = apply_D(gam, prod)
    if D == 0:
        return 'D annihilates the product', params
    cut0 = sp.expand(D.subs(lead, 0) - 1)
    if cut0 == -1:
        return 'cut vanishes identically on the degenerate locus', params
    R0 = sp.expand(R.subs(lead, 0) - 1)
    sol = sp.solve(sp.Eq(R0, 0), Qc[0])
    if len(sol) != 1:
        return f'resultant not linear in {Qc[0]}', params
    elim = sp.together(cut0.subs(Qc[0], sol[0]))
    num, _ = sp.fraction(sp.cancel(elim))
    num = sp.expand(num)
    # strip any bare factor of b (b = 0 is off the variety, since c0 = 1/b^k)
    p = sp.Poly(num, b)
    while p.degree() > 0 and p.coeff_monomial(1) == 0:
        p = sp.Poly(sp.quo(p.as_expr(), b, b), b)
    deg = p.degree()
    coeffs = {d: sp.factor(p.coeff_monomial(b**d)) for d in range(deg + 1)}
    if verbose:
        print(f'  eliminant degree {deg} in {b}')
        for d in sorted(coeffs, reverse=True):
            if coeffs[d] != 0:
                print(f'    [{b}^{d}] {coeffs[d]}')
    if deg == 1:
        root = sp.simplify(sp.solve(p.as_expr(), b)[0])
        if root.is_number:
            return f'MIRACLE for all moduli: {b} = {root}', params
        return f'degree 1 but PUNCTURE ({b} = {root})', params
    # Is there any parameter value making the degree drop to 1? That needs
    # every coefficient above degree 1 to vanish IDENTICALLY in c1..ck.
    cvars = [c for c in Qc[1:]]
    conds = []
    for d in range(2, deg + 1):
        cf = sp.expand(coeffs[d])
        if cf == 0:
            continue
        if not cf.free_symbols & set(params):
            if cf.is_number and cf != 0:
                return (f'degree {deg}: the [{b}^{d}] coefficient is the '
                        f'CONSTANT {cf}, so no parameter value can reduce '
                        f'the degree -- NEVER a miracle'), params
            conds.append(cf)
        else:
            conds.append(cf)
    if not conds:
        return f'degree {deg}, no reducing conditions found', params
    # solve: all higher coefficients vanish identically in the c-variables
    eqs = []
    for cf in conds:
        pc = sp.Poly(cf, *cvars) if cvars else None
        eqs += [sp.expand(t) for t in (pc.coeffs() if pc else [cf])]
    eqs = [e for e in set(eqs) if e != 0]
    if not eqs:
        return f'degree {deg}: reduces for all moduli (unexpected)', params
    if not params:
        return f'degree {deg} always -- NEVER a miracle (no moduli)', params
    sols = sp.solve(eqs, params, dict=True)
    if not sols:
        return (f'degree {deg} for EVERY value of {params} -- no cross-ratio '
                f'gives a miracle'), params
    return f'degree drops for {sols} -- CHECK THESE VALUES', params


def main():
    print('symbolic cross-ratio analysis, single cut, j = 1')
    print('(classes with <= 3 distinct roots have no moduli and are already')
    print(' settled by the sweep; only classes with 4+ appear below)')
    for k in (3, 4, 5):
        d = 1 + k
        print(f'\n=== (1,{k}), target degree {d}, covering degree {d}')
        for part in partitions(d):
            if len(part) < 4:
                continue
            lab = '(' + ','.join(map(str, part)) + ')'
            verdict, params = analyze(k, part, verbose=False)
            print(f'  {lab:>14} moduli {len(params)}: {verdict}')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        k = int(sys.argv[1])
        part = tuple(int(x) for x in sys.argv[2].split(','))
        print(analyze(k, part)[0])
    else:
        main()

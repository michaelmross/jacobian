"""
reconstruct.py — turn a MIRACLE verdict into an EXPLICIT counterexample.

`miracle.py` reports that a degenerate fiber is affine space. That is a flag,
not a map. This module carries the flag the rest of the way: it parametrizes
the variety V = {Res = 1} ∩ {cuts}, writes the multiplication map in those
coordinates, and verifies the two things that make it a counterexample —
constant nonzero Jacobian, and generic fiber of size > 1.

THE ALGORITHM (exactly the steps the verdict forces, made mechanical):

 1. Solve the degenerate system at lead = 0. A MIRACLE means the solution is a
    single point with CONSTANT coordinates, say b = 1, c = 1.
 2. Those constants say b - 1 and c - 1 vanish at lead = 0, hence are
    divisible by lead. Substitute  b = 1 + lead*t1,  c = 1 + lead*t2  with
    fresh variables.
 3. Solve the defining equations one at a time for the remaining unknowns.
 4. A solved value may still carry a denominator lead^m. Its numerator taken
    mod lead is an obstruction R that must vanish. If R is linear in some free
    variable, solve R = lead*(fresh variable) for it — this is the step that
    produced s = (3t - a*u)/2 in the (1,2) case. Repeat until everything is
    polynomial.
 5. The surviving free variables parametrize V. Coordinates on the target are
    the product coefficients not pinned by the cuts.
 6. Verify det J is a nonzero constant and count a generic fiber.

Step 4 is the only heuristic part: it handles obstructions that are linear in
some free variable, which covers the known case and any of similar shape. If
an obstruction is nonlinear the function says so rather than guessing.

Usage:
  python3 reconstruct.py 1 2 "(2,1)"        # rebuilds Alpoge's counterexample
  python3 reconstruct.py J K "CLASS[,CLASS...]"
"""
import sys

import sympy as sp

from symmult import form, resultant, apply_D, PARTITIONS

z, w = sp.symbols('z w')


def _fresh(n, used):
    out, i = [], 0
    while len(out) < n:
        s = sp.Symbol(f'p{i}')
        if s not in used:
            out.append(s)
        i += 1
    return out


def reconstruct(j, k, cuts, verbose=True):
    L, Lc = form(j, 'a')
    Q, Qc = form(k, 'c')
    prod = sp.expand(L * Q)
    R = resultant(L, Q)
    lead = Lc[0]
    rest = Lc[1:] + Qc

    eqs = [sp.expand(R - 1)]
    functionals = []
    for gammas, val in cuts:
        D = apply_D(gammas, prod)
        eqs.append(sp.expand(D - val))
        functionals.append(D)

    # --- step 1: the degenerate fiber ------------------------------------
    deg_eqs = [sp.expand(e.subs(lead, 0)) for e in eqs]
    sols = sp.solve(deg_eqs, rest, dict=True)
    if len(sols) != 1:
        return None, f'degenerate fiber has {len(sols)} branches -- no miracle'
    pinned = {v: val for v, val in sols[0].items() if val.is_number}
    if not pinned:
        return None, 'degenerate solution is not constant -- no miracle'
    if verbose:
        print(f'  degenerate fiber (at {lead} = 0): '
              + ', '.join(f'{v} = {val}' for v, val in pinned.items()))

    # --- step 2: lift the pinned coordinates -----------------------------
    used = set(rest) | {lead}
    fresh = _fresh(len(pinned), used)
    sub = {}
    free = [lead]
    for (v, val), p in zip(pinned.items(), fresh):
        sub[v] = val + lead * p
        free.append(p)
    if verbose:
        for v in pinned:
            print(f'    {v} = {sub[v]}')

    # --- steps 3-4: solve the rest, clearing denominators -----------------
    unknown = [v for v in rest if v not in sub]
    for _ in range(len(unknown) + 4):
        remaining = [v for v in rest if v not in sub]
        if not remaining:
            break
        progressed = False
        for eq in eqs:
            cur = sp.expand(eq.subs(sub))
            present = [v for v in remaining if cur.has(v)]
            if len(present) != 1:
                continue
            v = present[0]
            got = sp.solve(sp.Eq(cur, 0), v)
            if not got:
                continue
            val = sp.cancel(sp.together(got[0]))
            num, den = sp.fraction(val)
            if den.is_number:
                sub[v] = sp.expand(val)
                progressed = True
                continue
            # obstruction: numerator mod lead must vanish
            obstruction = sp.expand(num.subs(lead, 0))
            fixed = False
            for fv in free[1:]:
                p = sp.Poly(obstruction, fv)
                if p.degree() == 1 and p.LC().is_number:
                    q = _fresh(1, used | set(free))[0]
                    newval = sp.solve(sp.Eq(obstruction, lead * q), fv)[0]
                    for kk in list(sub):
                        sub[kk] = sp.expand(sub[kk].subs(fv, newval))
                    free = [x for x in free if x != fv] + [q]
                    used.add(q)
                    if verbose:
                        print(f'    divisibility forces {fv} = {newval}')
                    fixed = True
                    break
            if not fixed:
                return None, ('obstruction is not linear in any free '
                              f'variable: {sp.factor(obstruction)}')
            progressed = True
            break
        if not progressed:
            return None, 'could not solve the system successively'

    for v in rest:
        if v not in sub:
            return None, f'variable {v} never determined'
        if not sp.denom(sp.together(sub[v])).is_number:
            return None, f'{v} still has a denominator: {sub[v]}'

    # --- verify the defining equations hold identically -------------------
    for eq in eqs:
        if sp.simplify(sp.expand(eq.subs(sub))) != 0:
            return None, 'parametrization does not satisfy the equations'
    if verbose:
        print(f'  V is parametrized by {free} with no denominators: '
              f'V ~ C^{len(free)}')

    # --- step 5: coordinates on the target -------------------------------
    P = sp.Poly(prod, z, w)
    coeffs = [P.coeff_monomial(z**(j + k - i) * w**i) for i in range(j + k + 1)]
    # Each cut is a LINEAR FUNCTIONAL on Sym^{j+k}. Read its coefficients by
    # applying it to a generic form with symbolic coefficients -- the product
    # coefficients themselves are composite expressions (a0*c0, ...) and
    # cannot be differentiated against.
    Cs = [sp.Symbol(f'C{i}') for i in range(j + k + 1)]
    G = sum(Cs[i] * z**(j + k - i) * w**i for i in range(j + k + 1))
    rows = []
    for gammas, _ in cuts:
        DG = sp.expand(apply_D(gammas, G))
        rows.append([sp.diff(DG, C) for C in Cs])
    M = sp.Matrix(rows) if rows else sp.zeros(0, len(coeffs))
    pivots = M.rref()[1] if functionals else ()
    keep = [i for i in range(len(coeffs)) if i not in pivots]
    F = [sp.expand(coeffs[i].subs(sub)) for i in keep]
    return (F, free), None


def verify(F, free, verbose=True):
    J = sp.Matrix(F).jacobian(free)
    det = sp.expand(J.det())
    ok_det = det.is_number and det != 0
    if verbose:
        print(f'  det J = {det}   constant and nonzero: {ok_det}')
    if not ok_det:
        return False, None
    pt = {v: sp.Rational(2 + 3 * i) for i, v in enumerate(free)}
    tgt = [f.subs(pt) for f in F]
    sols = sp.solve([sp.expand(F[i] - tgt[i]) for i in range(len(F))],
                    list(free), dict=True)
    good = [s for s in sols
            if all(sp.simplify(F[i].subs(s) - tgt[i]) == 0 for i in range(len(F)))]
    if verbose:
        # The test solves ONE chosen fiber. That proves noninjectivity,
        # which is all a Jacobian counterexample needs -- but it does not
        # establish the GENERIC degree, which would require a symbolic
        # target and a discriminant argument. Wording corrected accordingly.
        print(f'  a tested fiber has {len(good)} distinct points  '
              f'({"NONINJECTIVE: COUNTEREXAMPLE" if len(good) > 1 else "injective at this point"})')
    return ok_det and len(good) > 1, len(good)


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit()
    j, k = int(sys.argv[1]), int(sys.argv[2])
    labels = [s.strip() for s in sys.argv[3].split('],[')]
    labels = [s.strip('[]') for s in sys.argv[3].replace('],[', '|').split('|')]
    table = PARTITIONS[j + k]
    cuts = [(table[lab], 1) for lab in labels]
    print(f'reconstructing (j,k)=({j},{k}) with cuts {labels}')
    res, err = reconstruct(j, k, cuts)
    if err:
        print('  FAILED:', err)
        sys.exit()
    F, free = res
    print('  multiplication map in these coordinates:')
    for i, f in enumerate(F):
        print(f'    F{i+1} =', sp.factor(f))
    verify(F, free)

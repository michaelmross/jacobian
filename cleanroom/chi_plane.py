"""
chi_plane.py — Route A for Lemma S: compute chi(V) as a TWO-VARIABLE
constructible-set Euler characteristic, with no fiber-constancy hypothesis.

WHY THIS EXISTS. hcheck.py computes chi(V_0) and the generic-t fiber chi
separately, and then needs hypothesis (H) ("every t != 0 fiber has chi = 1")
to assemble them. Fibering instead over the whole (t, a1)-plane eliminates
(H): with a = t fixed AND a1 fixed, every remaining equation is linear in c,
so the fiber is an affine subspace (chi = 1) whenever the system is
consistent and empty otherwise. Hence

    chi(V) = chi(S),   S = { (t, a1) in C^2 : M(t,a1) c = v is consistent },

with NO hypothesis. chi(S) = chi(C^2) - chi(Inc) = 1 - chi(Inc).

THE STRUCTURE THAT MAKES Inc COMPUTABLE (apolarity, verified in THEOREM.md):
row 1 of M is evaluation at p = [a1 : -t], the root of L; each cut row is
t*rho_r + a1*sigma_r. For a single cut (m = 1) the matrix has 2 rows, so

    Inc  ⊆  { rank M <= 1 }  =  V( 2x2 minors ),

and rank M <= 1 away from the origin means the cut functional is
PROPORTIONAL to evaluation at p — which happens only on the suspect lines
    t = 0        (the root g = 0, i.e. the operator d_z:  <L, .> = t)
    a1 = 0       (the root g = infinity, i.e. d_w:        <L, .> = a1)
    t = g*a1     (the finite root g, d_z - g d_w:          <L, .> = t - g a1)
NOTE ON LABELS: in the convention d_g = d_z - g d_w, the operator d_z is
g = 0 and d_w is g = infinity. This file's internal tags 'inf' and None mean
d_z and d_w respectively -- i.e. they are the OPPOSITE of the g-values. An
earlier draft of the comments (and of the paper) had these labels reversed;
the equations were always correct.
At a proportional point row2 = lam * row1, and the system reads
row1.c = 1 together with lam*(row1.c) = v: CONSISTENT iff lam = v.
At the origin row1 = 0 while the Res equation demands 1: always
INCONSISTENT. That origin point is the (t,a1)-plane's version of the
"a1 = 0 is always inconsistent" fact.

WHAT THIS FILE DOES NOT DO. It handles components that are LINES, which is
what the dictionary predicts and what every computed case produces. If a
component of higher degree appears it reports NONLINEAR COMPONENT and
refuses to guess, rather than silently assuming chi = 1 for it.

CONTROLS (run first, always):
  (1,2) (2,1)    -> Inc = {t=0} minus {a1=2} ~ C*,  chi(Inc)=0, chi(V)=1
  (1,2) (1,1,1)  -> Inc = {(0,0)},                  chi(Inc)=1, chi(V)=0
  (1,2) (3)      -> chi(V) = 0
matching euler.py and hcheck.py, which were derived independently.

Usage:
  python3 chi_plane.py            # controls, then the single-cut tables
"""
from itertools import combinations

import sympy as sp

from symmult import apply_D, PARTITIONS

z, w = sp.symbols('z w')
t, s = sp.symbols('t a1')          # a = t, a1 = s


def rows(k, gammas_list, vals):
    """Res row = evaluation at [a1 : -t]; cut rows = t*rho + a1*sigma."""
    M = [[s**(k - i) * (-t)**i for i in range(k + 1)]]
    rhs = [sp.Integer(1)]
    for g, v in zip(gammas_list, vals):
        rho = [apply_D(g, z**(k - i + 1) * w**i) for i in range(k + 1)]
        sig = [apply_D(g, z**(k - i) * w**(i + 1)) for i in range(k + 1)]
        M.append([sp.expand(t * rho[i] + s * sig[i]) for i in range(k + 1)])
        rhs.append(sp.Integer(v))
    return sp.Matrix(M), sp.Matrix(rhs)


def _minors(Mat, size):
    out = []
    R, C = Mat.shape
    if size == 0 or size > min(R, C):
        return out
    for ri in combinations(range(R), size):
        for ci in combinations(range(C), size):
            m = sp.expand(Mat[ri, ci].det())
            if m != 0:
                out.append(m)
    return out


def _rank_at(Mat, sub):
    return Mat.subs(sub).applyfunc(sp.simplify).rank()


def _line_components(polys):
    """Codim-1 part of V(polys): factor the gcd. Returns (linear, nonlinear)."""
    if not polys:
        return [], []
    g = sp.Integer(0)
    for p in polys:
        g = sp.gcd(g, p)
    if g == 0 or g.is_number:
        return [], []
    lin, non = [], []
    for f, _ in sp.factor_list(g)[1]:
        fe = sp.expand(f)
        if fe.is_number:
            continue
        deg = sp.Poly(fe, t, s).total_degree()
        (lin if deg == 1 else non).append(sp.factor(fe))
    return lin, non


def _param_line(f):
    """Parametrize the line f = 0 by one variable. Returns (sub, param)."""
    if sp.Poly(f, s).degree() == 1:
        return {s: sp.solve(sp.Eq(f, 0), s)[0]}, t
    return {t: sp.solve(sp.Eq(f, 0), t)[0]}, s


def chi_inc(M, rhs, verbose=True):
    """chi of the inconsistency locus Inc ⊂ C^2, and a description."""
    Ma = M.row_join(rhs)
    r_gen, ra_gen = M.rank(), Ma.rank()
    notes = []
    if ra_gen > r_gen:
        return None, ['generically INCONSISTENT -- not handled by this route']
    # Inc lies where rank M drops below generic
    drop = _minors(M, r_gen)
    lin, non = _line_components(drop)
    if non:
        return None, [f'NONLINEAR COMPONENT {non} -- manual analysis required']
    chi_total = 0
    pieces = []
    handled_pts = set()
    for f in lin:
        sub, par = _param_line(f)
        Mv, Mav = M.subs(sub), Ma.subs(sub)
        # generic behavior along the line
        rg, rag = Mv.rank(), Mav.rank()
        generic_inc = (rag > rg)
        # special points ON the line where behavior flips
        cand = set()
        for Mat, r in ((Mv, rg), (Mav, rag)):
            mins = _minors(Mat, r)
            if mins:
                g = sp.Integer(0)
                for m in mins:
                    g = sp.gcd(g, sp.expand(m))
                if g != 0 and not g.is_number:
                    for ff, _ in sp.factor_list(sp.expand(g), par)[1]:
                        if sp.Poly(ff, par).degree() >= 1:
                            cand.add(sp.factor(ff))
        flips = 0
        for ff in sorted(cand, key=lambda e: (sp.Poly(e, par).degree(),
                                              sp.srepr(e))):
            pf = sp.Poly(ff, par)
            roots = ([sp.solve(sp.Eq(ff, 0), par)[0]] if pf.degree() == 1
                     else [sp.RootOf(ff, 0)])
            mult = pf.degree()
            for p in roots:
                rm = _rank_at(Mv, {par: p})
                rma = _rank_at(Mav, {par: p})
                inc_here = (rma > rm)
                if inc_here != generic_inc:
                    flips += mult
                    handled_pts.add((sp.factor(ff), str(sub)))
        # chi(line) = 1; remove/add the flipped points
        c = (1 - flips) if generic_inc else flips
        chi_total += c
        pieces.append(f'{f} = 0: {"INC" if generic_inc else "cons"} generically'
                      f', {flips} flip point(s), chi contribution {c}')
    # ISOLATED POINTS (codimension 2). The gcd only sees the curve part of
    # V(minors); a rank drop confined to a point contributes nothing to the
    # gcd and was silently dropped by the first version of this routine --
    # which is precisely why (1,1,1) came back chi(V) = 1 instead of 0. Its
    # entire inconsistency locus IS the origin.
    try:
        pts = sp.solve(drop, [t, s], dict=True)
    except Exception:
        pts = []
    for pt in pts:
        if not pt or not all(getattr(v, 'is_number', False)
                             for v in pt.values()):
            continue                      # positive-dimensional: a line
        if len(pt) < 2:
            continue                      # underdetermined: a line
        if any(sp.simplify(f.subs(pt)) == 0 for f in lin):
            continue                      # already counted on a line
        rm, rma = _rank_at(M, pt), _rank_at(Ma, pt)
        if rma > rm:
            chi_total += 1
            pieces.append(f'isolated point {pt}: INCONSISTENT, '
                          f'chi contribution 1')
    if verbose:
        for p in pieces:
            notes.append(p)
    return chi_total, notes


def chi_plane(k, labels, vals=None, verbose=True):
    vals = vals or [1] * len(labels)
    gs = [PARTITIONS[k + 1][lab] for lab in labels]
    M, rhs = rows(k, gs, vals)
    ci, notes = chi_inc(M, rhs)
    if ci is None:
        if verbose:
            print(f'(1,{k}) {"+".join(labels)} {tuple(vals)}: UNRESOLVED')
            for n in notes:
                print(f'      {n}')
        return None
    chiV = 1 - ci
    if verbose:
        print(f'(1,{k}) {"+".join(labels)} {tuple(vals)}:')
        for n in notes:
            print(f'      {n}')
        print(f'      chi(Inc) = {ci}   chi(V) = 1 - {ci} = {chiV}   '
              + ('<-- MIRACLE candidate' if chiV == 1 else
                 'NOT affine space (unconditional, no hypothesis (H))'))
    return chiV


def controls():
    print('CONTROLS -- must match euler.py and hcheck.py:')
    exp = {('(2,1)',): 1, ('(1,1,1)',): 0, ('(3)',): 0}
    ok = True
    for labs, want in exp.items():
        got = chi_plane(2, list(labs))
        if got != want:
            print(f'   !! expected chi(V) = {want}, got {got}')
            ok = False
    print('CONTROLS ' + ('PASSED' if ok else 'FAILED -- stop here'))
    print('=' * 68)
    return ok


if __name__ == '__main__':
    if not controls():
        raise SystemExit(1)
    print()
    for k in (3, 4, 5):
        if k + 1 not in PARTITIONS:
            continue
        print(f'--- single cut, k = {k}')
        for lab in PARTITIONS[k + 1]:
            chi_plane(k, [lab])
        print()

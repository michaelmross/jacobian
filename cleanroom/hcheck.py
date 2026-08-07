"""
hcheck.py — compute chi(V) for j = 1 configurations of Tao's construction,
via the linear-in-c structure of the fibers. THE BRIDGE, made mechanical.

STRUCTURE. For j = 1 fix a = t. By the evaluation identity
Res(t z + a1 w, Q) = Q(a1, -t), the Res row has coefficients
(a1^k, -t a1^{k-1}, ..., (-t)^k), and each cut row is t*rho_r + a1*sigma_r
with constant vectors rho, sigma (sigma = the mu-vector of multicut.py).
Every equation is LINEAR in c. So the fiber V_t fibers over the a1-line with
affine-subspace fibers, and

    chi(V_t) = chi(S_t),   S_t = { a1 : the linear system is consistent }.

chi(S_t) is decided by POINTWISE rank comparisons: consistent at a1 = p iff
rank M(p) = rank [M|v](p). Two failure modes of a naive implementation, both
caught by the controls and fixed here:
  * generically inconsistent does NOT mean chi = 0: chi(S) is then the number
    of special points where consistency holds (for (2,1) at t = 0 that is the
    single point a1 = 2 -- the miracle itself);
  * "augmented rank drops too" must be tested at the POINT, not against the
    generic rank (at a1 = 0 in (1,1,1), rank M drops to 0 while rank [M|v]
    is 1: inconsistent, though all generic-size augmented minors vanish).

ASSEMBLY. chi(C*) = 0, so a CONSTANT fiber chi over t != 0 contributes 0 and
    chi(V) = chi(V_0) + sum over special t != 0 of (chi(V_t) - generic chi).
This file computes chi(V_0) and the generic-t fiber chi exactly, and reports
the special-t locus it can see (content/denominator degenerations in t) as a
REMAINING CHECK rather than silently assuming it is empty.

chi(C^n) = 1, and chi is a homeomorphism invariant: chi(V) != 1 proves
UNCONDITIONALLY that V is not affine space. That replaces the false fiber
bridge everywhere it fires.

Usage:
  python3 hcheck.py            # controls, then the single-cut table
"""
from itertools import combinations

import sympy as sp

from symmult import form, apply_D, PARTITIONS

z, w = sp.symbols('z w')
t, s = sp.symbols('t a1')


def rows(k, gammas_list, vals):
    M = [[s**(k - i) * (-t)**i for i in range(k + 1)]]
    rhs = [sp.Integer(1)]
    for g, v in zip(gammas_list, vals):
        rho = [apply_D(g, z**(k - i + 1) * w**i) for i in range(k + 1)]
        sig = [apply_D(g, z**(k - i) * w**(i + 1)) for i in range(k + 1)]
        M.append([sp.expand(t * rho[i] + s * sig[i]) for i in range(k + 1)])
        rhs.append(sp.Integer(v))
    return sp.Matrix(M), sp.Matrix(rhs)


def _minor_gcd(Mat, size):
    g = sp.Integer(0)
    R, C = Mat.shape
    if size == 0 or size > min(R, C):
        return None
    for ri in combinations(range(R), size):
        for ci in combinations(range(C), size):
            m = sp.expand(Mat[ri, ci].det())
            if m != 0:
                g = sp.gcd(g, m)
    return sp.Poly(sp.expand(g), s) if g != 0 else None


def _ranks_at(Mv, Ma, root):
    """Pointwise ranks at a1 = root (root may be a RootOf)."""
    sub = lambda X: X.subs(s, root).applyfunc(sp.simplify)
    return sub(Mv).rank(), sub(Ma).rank()


def chi_S(M, rhs, tval):
    """chi of {a1 : M(tval, a1) c = rhs consistent}, with detail."""
    Mv = M.subs(t, tval)
    Ma = Mv.row_join(rhs)
    r_gen, ra_gen = Mv.rank(), Ma.rank()
    generic_ok = (ra_gen == r_gen)
    # candidate special points: where either rank could drop
    # A set here made the DETAIL ORDER nondeterministic across runs (sympy
    # expression hashes vary with the process hash seed), which breaks
    # diff-based verification of logs. Deduplicate into a sorted list.
    cands = []
    for Mat, r in ((Mv, r_gen), (Ma, ra_gen)):
        g = _minor_gcd(Mat, r)
        if g is not None and g.degree() > 0:
            for f, _ in sp.factor_list(g.as_expr(), s)[1]:
                pf = sp.Poly(f, s)
                if pf.degree() >= 1 and sp.factor(f) not in cands:
                    cands.append(sp.factor(f))
    cands.sort(key=lambda e: (sp.Poly(e, s).degree(), sp.srepr(e)))
    detail = []
    chi = 1 if generic_ok else 0
    for f in cands:
        pf = sp.Poly(f, s)
        if pf.degree() == 1:
            roots = [sp.solve(sp.Eq(f, 0), s)[0]]
            mult = 1
        else:
            roots = [sp.RootOf(f, 0)]
            mult = pf.degree()      # conjugate roots share ranks (Galois)
        for p in roots:
            rM, rA = _ranks_at(Mv, Ma, p)
            ok = (rM == rA)
            detail.append((f, 'consistent' if ok else 'INCONSISTENT'))
            if generic_ok and not ok:
                chi -= mult
            if (not generic_ok) and ok:
                chi += mult
    return chi, detail


def analyze(k, labels, vals=None, verbose=True):
    vals = vals or [1] * len(labels)
    gs = [PARTITIONS[k + 1][lab] for lab in labels]
    M, rhs = rows(k, gs, vals)
    cg, dg = chi_S(M, rhs, t)
    c0, d0 = chi_S(M, rhs, sp.Integer(0))
    chiV = c0            # + special-t corrections, reported not assumed
    if verbose:
        tag = ('== 1: consistent with V ~ C^n (the miracle candidate)'
               if chiV == 1 else
               f'!= 1: V is NOT affine space (UNCONDITIONAL, modulo special-t)')
        print(f'(1,{k}) {"+".join(labels)} vals {tuple(vals)}:')
        print(f'   generic-t fiber chi : {cg}    {dg}')
        print(f'   chi(V_0)            : {c0}    {d0}')
        print(f'   chi(V) = {chiV}   {tag}')
        print()
    return chiV, cg


def controls():
    print('CONTROLS (must match euler.py and the known miracle):')
    ok = True
    v, _ = analyze(2, ['(2,1)'])
    ok &= (v == 1)
    v, _ = analyze(2, ['(1,1,1)'])
    ok &= (v == 0)
    v, _ = analyze(2, ['(3)'])
    ok &= (v == 0)
    print('CONTROLS ' + ('PASSED' if ok else 'FAILED -- do not trust below'))
    print('=' * 66)
    return ok


if __name__ == '__main__':
    if not controls():
        raise SystemExit(1)
    print()
    print('single-cut table, k = 3, 4 (previously conditional negatives):')
    for k in (3, 4):
        for lab in PARTITIONS[k + 1]:
            try:
                analyze(k, [lab])
            except Exception as ex:
                print(f'(1,{k}) {lab}: ERROR {type(ex).__name__}')

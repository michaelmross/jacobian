"""
symmult.py — Tao's multiplication-map construction, generalized.

Tao's July 2026 "digestion" post explains the counterexample conceptually: it
is the polynomial MULTIPLICATION map

    F : Sym^1(C^2) x Sym^2(C^2) -> Sym^3(C^2),    F(L, Q) = L*Q

restricted to {Res(L,Q) = 1} and then to a hyperplane in the target. Its
three properties:

  (b) NOT globally injective, because a generic cubic C = L1 L2 L3 splits as
      (L1, L2 L3), (L2, L1 L3), (L3, L1 L2) — three preimages. The covering
      degree 3 is the number of splittings, not a coincidence.
  (a) locally injective, because Res != 0 separates the root of L from the
      roots of Q, so the splitting can be reconstructed from the product.
  (c) the resulting variety is isomorphic to C^3 — the "miracle", and it
      occurs precisely when the hyperplane's dual operator has the form
      D = D1^2 D2 (two coincident roots).

GENERALIZATION. Nothing forces (j, k) = (1, 2). For any j, k:

    F : Sym^j x Sym^k -> Sym^{j+k},   covering degree C(j+k, j)

since the domain is an ORDERED pair, so every choice of which j roots go to
the first factor is a distinct preimage. Dimensions: {Res = 1} has dimension
j+k+1, matching the target, so one may then cut by hyperplanes.

| (j,k) | {Res=1} dim | cover deg | status |
|-------|-------------|-----------|--------|
| (1,1) | 3 | 2 | FAILS: {Res=1} is SL_2(C), not affine space |
| (1,2) | 4 | 3 | the known counterexample |
| (1,3) | 5 | 4 | UNTESTED |
| (2,2) | 5 | 6 | UNTESTED |
| (1,4) | 6 | 5 | UNTESTED |
| (2,3) | 6 | 10 | UNTESTED |

**(1,1) explains our conic-cone result.** The census found the conic cone
EMPTY at every degree. Here is why: covering degree 2 needs (j,k) = (1,1),
where Res(L1, L2) = a1 b2 - a2 b1, so {Res = 1} is exactly SL_2(C). Its
dimension already equals the target's, so no hyperplane cut remains to fix
anything, and SL_2(C) retracts to SU(2) = S^3 — not contractible, hence not
isomorphic to C^3. Property (c) fails outright. Minimal covering degree 2 is
blocked for a topological reason, not an arithmetic one.

**What to run.** For each (j,k) beyond (1,2), and each number of cuts, decide
whether the sliced variety is isomorphic to affine space. One need NOT land
in C^3: a counterexample in C^n for any n >= 3 is still a counterexample, so
cutting fewer times is allowed and widens the search. If any (j,k) works, it
gives a counterexample of covering degree C(j+k, j) — an infinite FAMILY of
inequivalent counterexamples indexed by (j,k), which is what the deformation
campaign failed to find near the known specimen.

Necessary conditions that are cheap to check first (a variety isomorphic to
affine space must satisfy all of these):
  * smooth everywhere;
  * trivial Picard group and trivial units (only constants);
  * contractible — in particular the (1,1) obstruction above.
"""
import sys
import sympy as sp
from math import comb

z, w = sp.symbols('z w')


def form(deg, name):
    """Generic homogeneous form of given degree in (z,w), with coefficients."""
    cs = [sp.Symbol(f'{name}{i}') for i in range(deg + 1)]
    return sum(c * z**(deg - i) * w**i for i, c in enumerate(cs)), cs


def resultant(F1, F2):
    """Homogeneous resultant, via dehomogenisation in w."""
    f1 = sp.Poly(F1.subs(w, 1), z)
    f2 = sp.Poly(F2.subs(w, 1), z)
    return sp.expand(sp.resultant(f1, f2))


def setup(j, k, verbose=True):
    """Build the variety {Res = 1} and the multiplication map in coordinates."""
    L, Lc = form(j, 'a')
    Q, Qc = form(k, 'c')
    prod = sp.expand(L * Q)
    P = sp.Poly(prod, z, w)
    coeffs = [P.coeff_monomial(z**(j + k - i) * w**i) for i in range(j + k + 1)]
    R = resultant(L, Q)
    if verbose:
        print(f'(j,k) = ({j},{k})')
        print(f'  domain coords: {len(Lc)} + {len(Qc)} = {len(Lc)+len(Qc)}')
        print(f'  target dim   : {j+k+1}')
        print(f'  Res(L,Q)     : {sp.factor(R)}')
        print(f'  {{Res=1}} dim  : {len(Lc)+len(Qc)-1}')
        print(f'  covering deg : {comb(j+k, j)}')
    return dict(L=L, Q=Q, Lc=Lc, Qc=Qc, coeffs=coeffs, Res=R)


def check_11():
    """The (1,1) obstruction, made explicit."""
    d = setup(1, 1, verbose=False)
    print('(1,1): Res(L1,L2) =', sp.factor(d['Res']))
    print('  {Res = 1} is SL_2(C) in the four coordinates (a0,a1,c0,c1).')
    print('  dim 3 = target dim, so NO hyperplane cut is available.')
    print('  SL_2(C) retracts to SU(2) = S^3, not contractible,')
    print('  so it is not isomorphic to C^3. Covering degree 2 is BLOCKED.')



# ---------------------------------------------------------------------------
# MIRACLE TEST (the computation that decides each (j,k) case)
#
# Tao's property (c) turns on what happens at the degenerate locus a0 = 0,
# where L's root runs to infinity. For (1,2) the two equations there were
# c*b^2 = 1 and b*c = 1, with the UNIQUE affine solution b = c = 1 -- versus
# the six Bezout predicts, the other five sitting at infinity. That unique,
# CONSTANT solution is what makes the a0 = 0 fiber a clean copy of affine
# space and lets the global chart close up.
#
# So the test on any candidate cut is: at a0 = 0, does {Res = 1} + {cut = 1}
# have exactly one solution, with POLYNOMIAL (not Laurent) values? A Laurent
# solution such as a1 = 1/c1 means the fiber is punctured -- C^* x C^2, not
# C^3 -- and the miracle fails just as surely as multiple branches do.
#
# Hyperplanes are classified by their dual operator's root multiplicities.
# For cubics the partitions of 3 are (3), (2,1), (1,1,1), and Tao's miracle
# is at (2,1). For quartics the partitions of 4 are (4), (3,1), (2,2),
# (2,1,1), (1,1,1,1) -- and NOTE that a single-coefficient cut only reaches
# the two-root partitions (4), (3,1), (2,2). The mixed classes need genuine
# linear combinations, so scanning coefficients alone silently skips them.
#
# RESULT for (1,3), single cut (target C^4, covering degree 4):
#   (4)        cut vanishes identically at a0 = 0
#   (3,1)      2 branches (a1 = +-sqrt(6))
#   (2,2)      1 branch but LAURENT: a1 = 1/(4c1) -- punctured fiber
#   (2,1,1)    eliminant 4 a1^3 c1 + a1^2 - 6 = 0, cubic: 3 branches
#   (1,1,1,1)  eliminant 12 a1^3 (c1-c2) + a1^2 - 6 = 0, cubic: 3 branches
# No class produces the miracle, so (1,3) with one cut does not give a
# counterexample in C^4. Two-cut (1,3) and all of (2,2) remain UNTESTED.
# ---------------------------------------------------------------------------


def apply_D(gammas, expr):
    """Apply D = prod_i (d_z - gamma_i d_w); 'inf' means d_z, None means d_w."""
    e = expr
    for g in gammas:
        if g == 'inf':
            e = sp.diff(e, z)
        elif g is None:
            e = sp.diff(e, w)
        else:
            e = sp.diff(e, z) - g * sp.diff(e, w)
    return sp.expand(e)


def miracle_test(j, k, gammas, verbose=True):
    """Test one hyperplane class for the affine miracle at a0 = 0."""
    L, Lc = form(j, 'a')
    Q, Qc = form(k, 'c')
    prod = sp.expand(L * Q)
    R = resultant(L, Q)
    lead = Lc[0]
    cut = apply_D(gammas, prod)
    if cut == 0:
        return 'D annihilates the product'
    cut0 = sp.expand(cut.subs(lead, 0))
    if cut0 == 0:
        return 'cut vanishes identically on the degenerate locus'
    R0 = sp.expand(R.subs(lead, 0))
    # On the degenerate locus the resultant collapses to a1^k * c0 = 1, so
    # substitute c0 = 1/a1^k and clear denominators. This ELIMINATION is
    # robust and fast; sp.solve hangs on the mixed-root classes.
    b, q0 = Lc[1], Qc[0]
    c0_val = sp.solve(sp.Eq(R0, 1), q0)
    if len(c0_val) != 1:
        return f'degenerate resultant: {sp.factor(R0)} = 1'
    elim = sp.together(sp.expand((cut0 - 1).subs(q0, c0_val[0])))
    num, _ = sp.fraction(sp.cancel(elim))
    num = sp.expand(num)
    if not num.free_symbols or b not in num.free_symbols:
        return 'cut is independent of the degenerate coordinate'
    p = sp.Poly(num, b)
    # strip any factor of b (b = 0 is not on the variety, since c0 = 1/b^k)
    while p.degree() > 0 and p.coeff_monomial(1) == 0:
        p = sp.Poly(sp.quo(p.as_expr(), b, b), b)
    deg = p.degree()
    if deg == 0:
        return 'cut is independent of the degenerate coordinates'
    if deg > 1:
        return (f'{deg} branches -- no miracle  '
                f'(eliminant {sp.factor(p.as_expr())} = 0)')
    root = sp.solve(p.as_expr(), b)[0]
    if not sp.simplify(root).is_number:
        return (f'unique branch but LAURENT/parametric (punctured fiber) '
                f'-- no miracle   [b = {sp.simplify(root)}]')
    c0r = sp.simplify(c0_val[0].subs(b, root))
    return (f'MIRACLE: unique polynomial solution '
            f'{{{b}: {root}, {q0}: {c0r}}}')



# ---------------------------------------------------------------------------
# GENERAL miracle test: works for any (j,k) and ANY NUMBER of cuts.
#
# The elimination version above assumes Res is linear in c0, which holds for
# j = 1 but fails for (2,2), where Res|_{a0=0} is quadratic in c0. The general
# criterion avoids solving altogether:
#
#     the degenerate fiber is an affine SUBSPACE
#       <=>  its ideal has a Groebner basis of LINEAR forms.
#
# Check against Tao's case: at a0 = 0 the ideal is (c b^2 - 1, b c - 1), whose
# Groebner basis is (b - 1, c - 1) -- linear, so the fiber is a clean copy of
# affine space. Contrast (1,3) class (2,2), whose ideal (a1^3 c0 - 1,
# 4 a1 c1 - 1) keeps the nonlinear generator a1 c1 - 1/4: the fiber is
# punctured, not affine. Linearity of the Groebner basis captures BOTH the
# "multiple branches" and the "Laurent/punctured" failures in one test.
#
# HOW TO EXTEND:
#   * more partitions -> use partition_gammas(part) below, which builds the
#     gamma list from any partition of d. Distinct parts get distinct roots.
#   * more cuts -> pass a list of (gammas, value) pairs. To land in C^n from
#     {Res=1} of dimension d+1, use (d+1-n) cuts.
#   * j >= 2 -> nothing extra needed; this test does not assume j = 1.
# ---------------------------------------------------------------------------


def partition_gammas(part):
    """Build a gamma list from a partition, e.g. (2,1,1) -> distinct roots
    with multiplicities 2,1,1. Roots used: d_z, d_w, then 1, 2, 3, ..."""
    roots = ['inf', None, 1, 2, 3, 4, 5]
    out = []
    for idx, mult in enumerate(part):
        out += [roots[idx]] * mult
    return out


def miracle_test_general(j, k, cuts, verbose=False):
    """cuts: list of (gammas, value). Returns a verdict string.

    Criterion: at the degenerate locus (leading coefficient of L = 0), the
    ideal generated by {Res - 1} and {D_i(product) - value_i} must have a
    Groebner basis consisting of LINEAR forms.
    """
    L, Lc = form(j, 'a')
    Q, Qc = form(k, 'c')
    prod = sp.expand(L * Q)
    R = resultant(L, Q)
    lead = Lc[0]
    gens = [sp.expand(R.subs(lead, 0) - 1)]
    for gammas, val in cuts:
        d = apply_D(gammas, prod)
        if d == 0:
            return 'a cut annihilates the product'
        d0 = sp.expand(d.subs(lead, 0) - val)
        if d0 == -val:
            return 'a cut vanishes identically on the degenerate locus'
        gens.append(d0)
    variables = [v for v in Lc[1:] + Qc]
    try:
        G = sp.groebner(gens, *variables, order='lex')
    except Exception as e:
        return f'Groebner failed: {type(e).__name__}'
    exprs = list(G.exprs)
    if exprs == [sp.Integer(1)]:
        return 'EMPTY degenerate fiber (no solutions) -- no miracle'
    # CORRECT CRITERION: the fiber is isomorphic to affine space iff the lex
    # Groebner basis is TRIANGULAR -- every generator has degree 1 in its
    # leading variable AND constant leading coefficient, so each variable is
    # a polynomial function of the later ones with no branching and no
    # division. Requiring LINEAR generators is too strict: a graph such as
    # a1 - 13824*c0**2 is nonlinear yet describes affine space perfectly
    # well. The two failure modes it must still catch are
    #   * degree > 1 in the leading variable  -> BRANCHING (several points)
    #   * non-constant leading coefficient    -> PUNCTURE (Laurent solution)
    bad = []
    for e in exprs:
        p = sp.Poly(e, *variables)
        lead_var = next((v for v in variables if e.has(v)), None)
        if lead_var is None:
            continue
        pv = sp.Poly(e, lead_var)
        if pv.degree() > 1:
            bad.append(('BRANCHING', e, lead_var, pv.degree()))
        elif not pv.LC().is_number:
            bad.append(('PUNCTURE', e, lead_var, pv.LC()))
    if not bad:
        return f'MIRACLE: degenerate fiber is affine space  {exprs}'
    # report the most informative obstruction: fewest variables
    bad.sort(key=lambda b: len(b[1].free_symbols))
    kind, e, lv, extra = bad[0]
    if kind == 'BRANCHING':
        return (f'no miracle: BRANCHING -- {sp.factor(e)} has degree '
                f'{extra} in {lv} ({extra} sheets, fiber disconnected)')
    return (f'no miracle: PUNCTURE -- {sp.factor(e)} solves {lv} only where '
            f'{sp.factor(extra)} != 0 (fiber is punctured, not affine)')


def _partitions(n, maxpart=None):
    """All partitions of n as non-increasing tuples."""
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for first in range(min(n, maxpart), 0, -1):
        for rest in _partitions(n - first, first):
            yield (first,) + rest


def build_partitions(dmax=6):
    """PARTITIONS[d][label] = gamma list, for every partition of d."""
    out = {}
    for d in range(3, dmax + 1):
        out[d] = {}
        for part in _partitions(d):
            lab = '(' + ','.join(map(str, part)) + ')'
            out[d][lab] = partition_gammas(part)
    return out


PARTITIONS = build_partitions(6)


def scan(j, k):
    """Run the miracle test over every hyperplane class for this (j,k)."""
    d = j + k
    if d not in PARTITIONS:
        print(f'  no partition table for degree {d}; pass gammas manually '
              f'to miracle_test(j, k, gammas)')
        return
    print(f'(j,k) = ({j},{k})   target degree {d}, '
          f'covering degree {comb(d, j)}')
    for lab, gam in PARTITIONS[d].items():
        print(f'  {lab:>10}: {miracle_test(j, k, gam)}')


def scan_general(j, k, ncuts=1, deg=None):
    """Scan every hyperplane class (or pair, for ncuts=2) for the miracle."""
    from itertools import combinations_with_replacement
    d = deg or (j + k)
    if d not in PARTITIONS:
        print(f'  extend PARTITIONS for degree {d}, or use partition_gammas')
        return
    items = list(PARTITIONS[d].items())
    print(f'(j,k) = ({j},{k})  covering degree {comb(j+k, j)}  '
          f'cuts = {ncuts}')
    if ncuts == 0:
        r = miracle_test_general(j, k, [])
        print(f'  {"(no cut)":>10}: {r}')
        print(f'  {1 if r.startswith("MIRACLE") else 0} miracle(s) found '
              f'among 1 case')
        return
    if ncuts == 1:
        hits = 0
        for lab, g in items:
            r = miracle_test_general(j, k, [(g, 1)])
            if r.startswith('MIRACLE'):
                hits += 1
            print(f'  {lab:>10}: {r}')
        print(f'  {hits} miracle(s) found among {len(items)} classes')
    else:
        hits = 0
        for combo in combinations_with_replacement(items, ncuts):
            for vals in [(1,) * ncuts, (1,) + (0,) * (ncuts - 1)]:
                cuts = [(g, v) for (l, g), v in zip(combo, vals)]
                r = miracle_test_general(j, k, cuts)
                if r.startswith('MIRACLE'):
                    hits += 1
                    labs = '+'.join(l for l, _ in combo)
                    print(f'  {labs} values{vals}: {r}')
        print(f'  {hits} miracle(s) found among all class combinations')


if __name__ == '__main__':
    if len(sys.argv) >= 3:
        j, k = int(sys.argv[1]), int(sys.argv[2])
        n = int(sys.argv[3]) if len(sys.argv) > 3 else 1
        setup(j, k)
        print()
        scan_general(j, k, n)
    else:
        check_11()
        print()
        print('=' * 68)
        print('POSITIVE CONTROL: (1,2) must show the miracle at (2,1) only.')
        print('=' * 68)
        scan_general(1, 2, 1)
        print()
        print('=' * 68)
        print('(1,3) one cut -> C^4 (degree 4); (1,3) two cuts -> C^3;')
        print('(2,2) one cut -> C^4 (degree 6).')
        print('=' * 68)
        scan_general(1, 3, 1)
        print()
        scan_general(1, 3, 2)
        print()
        scan_general(2, 2, 1)
        print()
        print('TO EXTEND FURTHER:')
        print('  * degree 5+: add entries to PARTITIONS using')
        print('    partition_gammas((3,1,1)) etc.')
        print('  * more cuts: scan_general(j, k, ncuts)')
        print('  * command line: python3 symmult.py J K [NCUTS]')

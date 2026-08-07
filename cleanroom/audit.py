"""
audit.py — global certificates for the negative results, and an honest
classification of which claims are PROOF and which are EVIDENCE.

THE GAP (found by adversarial review, 2026-08-06). Every fiber-based
negative in this project rests on the bridge

    V isomorphic to C^n  ==>  V ∩ {a0 = 0} is affine space,

and that implication is FALSE for general pairs (variety, function): a
hypersurface {g = 0} inside C^n can be disconnected, singular, or punctured
while C^n is still C^n ({xy = 1} in C^2 is the standard example). The bridge
was asserted, never proven. Consequently:

PROOF (unconditional):
  * The positive direction: reconstruct.py builds an explicit counterexample
    from the (1,2)/(2,1) configuration -- constructive, verified.
  * Covering degree 2 is impossible: {Res=1} = SL_2(C) retracts to S^3, not
    contractible, so it is not C^3. GLOBAL argument.
  * The two-cut k=2 configurations are not counterexamples: a1 is a global
    nonconstant unit on V (Groebner certificate + two points). GLOBAL.
  * NEW (this file): every configuration whose cuts include the full-
    derivative class (k+1) is excluded: the cut reads (k+1)! a0 c0 = 1, so
    a0 is a global unit, and it is nonconstant (values 1 and 2 attained,
    certified by Groebner). GLOBAL.
  * Res|_{a0=0} = a1^k c0: provable for all k from the classical evaluation
    identity Res(a0 z + a1 w, Q) = Q(a1, -a0); at a0 = 0 this is c0 a1^k.
    (Verified symbolically for k = 2..6 in multicut.py; the identity itself
    is standard.)

PROOF, CONDITIONAL on the fiber bridge:
  * The j = 1 uniqueness "theorem" (multicut.py): unit argument on V_0 +
    left-kernel criterion + dimension count. Rigorous FROM the premise that
    a counterexample forces V_0 affine. Without the bridge it proves:
    "(1,2)/(2,1) is the unique configuration whose degenerate fiber is
    affine space" -- i.e., uniqueness within the degenerate-fiber METHOD,
    not yet uniqueness of counterexamples.

EVIDENCE (mod-p computation, no global certificate yet):
  * All remaining sweep negatives: BRANCHING and PUNCTURE verdicts certify
    that V_0 is not affine, which -- absent the bridge -- does not by itself
    preclude V being C^n.

UPGRADE PATHS, in order of expected yield:
  1. Reread Tao's post for his treatment of the FAILING cubic classes (3)
     and (1,1,1): if he proves those varieties are not C^3, his argument may
     generalize and supply the bridge (or replace it).
  2. Extend this audit beyond coordinate units: products a_i c_j, the cut
     functionals themselves, low-degree combinations.
  3. Smoothness: if V is singular it is not C^n; the Jacobian-criterion
     check is mechanical (heavier than the unit test but decidable).
  4. Topological invariants (Euler characteristic via degrees) where the
     unit test fails.

Usage:
  python3 audit.py            # run the audit over the small configurations
"""
import sympy as sp

from symmult import form, resultant, apply_D, PARTITIONS

z, w = sp.symbols('z w')


def global_unit_audit(j, k, cut_labels, vals, candidates=None, verbose=True):
    """Search for f with V ∩ {f=0} empty (1 in I+(f), Nullstellensatz) and f
    nonconstant on V (two attained values). Such an f is a nonconstant unit
    on V, and affine space has only constant units: V is NOT any C^n,
    UNCONDITIONALLY."""
    L, Lc = form(j, 'a')
    Q, Qc = form(k, 'c')
    prod = sp.expand(L * Q)
    R = resultant(L, Q)
    table = PARTITIONS[j + k]
    I = [sp.expand(R - 1)]
    for lab, v in zip(cut_labels, vals):
        d = apply_D(table[lab], prod)
        if d == 0:
            return 'cut annihilates the product'
        I.append(sp.expand(d - v))
    vars_ = Lc + Qc
    cands = candidates if candidates is not None else list(vars_)
    for f in cands:
        G = sp.groebner(I + [f], *vars_, order='grevlex')
        if list(G.exprs) != [sp.Integer(1)]:
            continue
        attained = []
        for t in (1, 2, 3, 5, 7):
            Gt = sp.groebner(I + [f - t], *vars_, order='grevlex')
            if list(Gt.exprs) != [sp.Integer(1)]:
                attained.append(t)
            if len(attained) >= 2:
                break
        if len(attained) >= 2:
            return (f'UNCONDITIONAL: {f} is a nonconstant unit on V '
                    f'(values {attained} attained) -- V is not any C^n')
        return f'{f} vanishes nowhere but may be constant -- inconclusive'
    return 'no certificate among the candidates -- verdict remains CONDITIONAL'


def main():
    tests = [
        (1, 2, ['(3)'], (1,), 'k=2 (3)       full-derivative cut'),
        (1, 2, ['(1,1,1)'], (1,), 'k=2 (1,1,1)'),
        (1, 3, ['(4)'], (1,), 'k=3 (4)       full-derivative cut'),
        (1, 3, ['(3,1)'], (1,), 'k=3 (3,1)     branching'),
        (1, 3, ['(2,2)'], (1,), 'k=3 (2,2)     punctured'),
        (1, 4, ['(5)'], (1,), 'k=4 (5)       full-derivative cut'),
        (1, 4, ['(3,2)'], (1,), 'k=4 (3,2)     punctured'),
    ]
    print('global-unit audit (unconditional certificates where found):')
    for j, k, labs, vals, desc in tests:
        print(f'  {desc}:')
        print(f'    {global_unit_audit(j, k, labs, vals)}')
    print()
    print('pattern: the full-derivative class (k+1) is excluded for every k')
    print('(cut = (k+1)! a0 c0 = 1 makes a0 a global unit). The branching and')
    print('punctured classes have no coordinate certificate; those negatives')
    print('remain conditional on the fiber bridge.')


if __name__ == '__main__':
    main()

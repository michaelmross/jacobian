# jacobian

Two notes on the Jacobian conjecture counterexample of Levent Alpöge (with
Claude Fable 5), analyzed through the construction described in Terence
Tao's [digestion of the
counterexample](https://terrytao.wordpress.com/2026/07/21/a-digestion-of-the-jacobian-conjecture-counterexample/).

### [trichotomy.md](trichotomy.md) — A cubic corollary

For one cut in the (j,k) = (1,2) construction, the three classes of cubic
operator behave as follows: **(3)** carries a nonconstant-unit obstruction,
**(2,1)** gives the affine miracle and Alpöge's counterexample, and
**(1,1,1)** carries an Euler-characteristic obstruction. Three short,
independent arguments using three different invariants.

### [uniqueness.md](uniqueness.md) — Uniqueness for j = 1

For Sym¹ × Sym^k with any number of linearly independent cuts and target
dimension at least three, the compactly supported Euler characteristic
satisfies χ(V) = (k−1)N. Since a counterexample requires χ(V) = 1, this
forces k = 2, one cut, class (2,1), with nonzero cut value — Alpöge's
configuration, uniquely.

### /cleanroom

Verification code and logs for both notes. Each note's own verification
section explains what the scripts establish, and — as importantly — what
they do not.

---

Michael M. Ross · michaelmross@cantab.net

# Uniqueness of the Jacobian counterexample within Tao's construction
## (j = 1, target dimension at least three)

*Self-contained. Version 3; see the revision history in §9.*

---

## 1. Background and statement

In July 2026 Levent Alpöge, with Claude Fable 5, produced a counterexample to
the Jacobian conjecture: a polynomial map F : C^3 -> C^3 of degree 7 with
det DF = -2 constant, 3-to-1 rather than injective. Terence Tao's subsequent
digestion explained it structurally, as follows.

Let Sym^d be the space of binary forms of degree d in (z, w), so
dim Sym^d = d + 1. Multiplication

        mu : Sym^j x Sym^k -> Sym^{j+k},     (L, Q) |-> L Q

has a positive-dimensional symmetry, the scaling
(L, Q) |-> (lam L, lam^{-1} Q), which fixes the product. One imposes the
resultant normalization Res(L, Q) = 1, and then m further affine hyperplane
conditions in the target. What remains is a variety V carrying a self-map.
**For k >= 2, if V is isomorphic to affine space C^n**, transporting mu
through such an isomorphism gives a polynomial self-map of C^n that is étale
— the Jacobian determinant is a nonzero constant — and finite-to-one of
degree greater than one: a Jacobian counterexample. Tao calls V ≅ C^n an
"affine miracle". The restriction to k >= 2 matters: as the scaling
computation just below shows, for k = 1 the normalized map is not generically
finite, so no degree language applies there at all.

**The covering degree, correctly.** A general form of degree j + k splits
into factors of degrees j and k in C(j+k, j) ways, but that counts splittings
before normalization. Take j = 1. Since Q is homogeneous of degree k and the
root of L is unchanged by scaling,

        Res(lam L, lam^{-1} Q) = lam^{k-1} Res(L, Q).

Hence for k >= 2 the normalization leaves k - 1 scalings (the (k-1)-st roots
of unity) in each orbit, and the normalized map has generic degree

        (k + 1)(k - 1) = k^2 - 1,

not k + 1. At k = 2 both equal 3, which is Tao's case and Alpöge's example.
For k = 1 the resultant is scaling-INVARIANT, {Res = 1} retains the whole
C^*-orbit, and the normalized map is not generically finite at all — so
there is no finite covering of degree two in this construction to rule out.

Tao states that for k = 2 with one cut the affine miracle occurs precisely in
the case of two identical roots, and verifies that case himself by an
explicit coordinate computation; his post contains no argument for the
negative half of "precisely", which his construction does not require. Those
two cases are proved in the companion note [trichotomy.md](trichotomy.md).)*, by arguments
independent of everything below. This note proves the corresponding
uniqueness statement for all of j = 1.

> **Theorem.** Fix k >= 1 and m >= 0 cuts whose defining functionals are
> LINEARLY INDEPENDENT, with arbitrary classes (including moduli, §2) and
> arbitrary values, and suppose n = k + 2 - m >= 3. If V is isomorphic to
> C^n then
>
>        k = 2,   m = 1,   the cut class is (2,1),   and its value is
>        nonzero.
>
> Any nonzero value normalizes to 1 by rescaling the functional, so this is
> the Alpöge–Tao configuration up to the coordinate and scaling actions.

The invariant is the compactly supported Euler characteristic chi_c. For
complex algebraic varieties chi_c agrees with the topological chi, is a
homeomorphism invariant, is additive over decompositions into constructible
pieces, and satisfies chi_c(C^n) = 1, chi_c(C^*) = 0. Exhibiting
chi_c(V) != 1 rules out V ≅ C^n unconditionally.

**Scope.** The hypothesis n >= 3 is used throughout, and is equivalent to
m <= k - 1, which the finiteness argument of §5 requires. The case n = 2
(m = k) is NOT covered here: a hypothetical affine miracle in dimension two
would bear on the two-dimensional Jacobian conjecture and is not excluded by
anything below.

---

## 2. Notation

Throughout j = 1. Write

        L = t z + s w
        Q = c_0 z^k + c_1 z^{k-1} w + ... + c_k w^k

so V sits in C^{k+3} with coordinates (t, s, c_0, ..., c_k).

**Resultant.** For L linear, Res(L, Q) = Q(s, -t): evaluation of Q at the
root of L. The first defining equation is

        (R)      sum_i c_i s^{k-i} (-t)^i = 1.

**Cuts.** A hyperplane condition on Sym^{k+1} is given by a constant-
coefficient differential operator D of order k+1, imposed as D(LQ) = v. Any
such D factors into k+1 directional derivatives,

        D = prod_{j=0}^{k} d_{g_j},      d_g = d_z - g d_w,

with d_z the case g = 0 and d_w the case g = infinity. The multiset of roots
{g_j} is the CLASS of the cut, recorded as the partition of k+1 of root
multiplicities. PGL_2(C) acts 3-transitively on P^1, so three roots normalize
and no more: a single class with four or more distinct roots carries moduli
(cross-ratios), and a COLLECTION of cuts carries joint moduli even when each
individual class has at most three distinct roots, since one PGL_2 action
cannot normalize them all independently. See the remark in §7.

**Independence.** The functionals D_1, ..., D_m are assumed linearly
independent on Sym^{k+1}, and m denotes their rank. This is forced by the
dimension hypothesis: a dependent collection either makes V empty or fails to
cut dimension m. It is also necessary for Lemma 2 below, which is false
without it (take G_1 = G_2 = G with l_p^k not dividing G: the two cut rows
coincide, so rank M <= m automatically, while span{G_1, G_2} = CG contains no
multiple of l_p^k).

**Linear structure.** For (t, s) FIXED, (R) and every cut are linear in
c = (c_0, ..., c_k). Write the system as

        M(t, s) c = v,       v = (1, v_1, ..., v_m)^T,

with M of size (1+m) x (k+1): row 0 from (R), row r from cut r. Write
v_cut = (v_1, ..., v_m) for the cut values alone; the two vectors are used in
different places below and should not be conflated. Then

        (t, s) CONSISTENT  <=>  rank M(t,s) = rank [M | v](t,s).

**Reduction to a plane computation.** Stratify C^2 into the locally closed
loci on which rank M is constant, and refine that stratification into
finitely many locally closed pieces on each of which a SPECIFIED maximal
minor is nonzero — a constant-rank locus need not admit one globally fixed
nonzero minor. Over each such piece, Gaussian elimination on that minor
presents the consistent locus as an
algebraically locally trivial bundle whose fibers are affine spaces (cosets
of ker M, of possibly varying dimension), and the inconsistent locus as
carrying empty fibers. Since chi_c of any affine space is 1, additivity gives

        chi_c(V) = chi_c(S),   S = {(t,s) : the system is consistent}.

S and its complement Inc are constructible, not in general locally closed;
chi_c is additive over constructible sets, so with chi_c(C^2) = 1,

        (*)     chi_c(V) = 1 - chi_c(Inc).

---

## 3. Apolarity dictionary

Functionals on Sym^k correspond to degree-k forms in dual variables (Z, W)
under the apolarity pairing <z^a w^b, Z^a W^b> = a! b!. Two entries suffice.

* **Row 0 is a pure power.** Evaluation Q |-> Q(alpha, beta) corresponds to
  (alpha Z + beta W)^k, since <Q, (alpha Z + beta W)^k> = k! Q(alpha, beta).
  By (R), row 0 corresponds to

        l_p^k,      l_p = s Z - t W,

  where p = [s : -t] in P^1 is **the root of L**.

* **Row r is a derivative.** Let G_r be the dual form of D_r — the product of
  the linear forms of its roots, of degree k+1. Then Q |-> D_r(L Q)
  corresponds to

        L(d) G_r,       L(d) = t d_Z + s d_W.

  (Leibniz: deg L = 1 and ord D_r = k+1, so exactly one derivative can fall
  on L.)

**Lemma 1.** L(d) l_p = t s + s (-t) = 0.

The operator defining the cut rows annihilates the form defining the Res row.

---

## 4. The degeneracy criterion

**Lemma 2.** Assume the D_r linearly independent. For (t, s) != (0, 0),

        rank M(t,s) <= m   <=>   span{G_1, ..., G_m} contains a nonzero
                                 element divisible by l_p^k.

*Proof.* Since (t,s) != 0, set u = l_p = sZ - tW and choose y = aZ + bW with
t a + s b = 1; this is possible because (t,s) != 0, and then

        det [[s, -t], [a, b]] = s b + t a = 1,

so (u, y) is a basis of the linear forms, and L(d) y = t a + s b = 1 while
L(d) u = 0 by Lemma 1. In these coordinates L(d) = d_y. Expanding a
degree-(k+1) form as G = sum_i b_i u^{k+1-i} y^i,

        L(d) G = sum_i i b_i u^{k+1-i} y^{i-1},

whose only u^k term is b_1 u^k. Row 0 is u^k. So L(d)G is a scalar multiple
of u^k exactly when b_i = 0 for all i >= 2, i.e. when
G = u^k (b_0 u + b_1 y), i.e. when u^k | G. Applying this to a combination
sum_r beta_r G_r and using independence of the D_r (so that a nonzero
coefficient vector gives a nonzero form) yields the statement for m > 1. QED

**Consequence.** deg G = k+1, so u^k | G forces p to be a root of the
combination of multiplicity >= k. For a SINGLE cut the class must be (k+1) or
(k,1). For m > 1 no such restriction on the individual classes follows: a
combination can acquire a k-fold root when no individual G_r has one. See §8.

**Lemma 3 (geometry).** p = [s : -t] depends only on the ratio, so degeneracy
is scaling-invariant: each degenerate p is a punctured LINE through the
origin of the (t,s)-plane, the line where the corresponding pairing vanishes.
With d_g = d_z - g d_w the pairings are

        g = 0        : <L, .> = t,      so the line is t = 0
        g = infinity : <L, .> = s,      so the line is s = 0
        g finite     : <L, .> = t - g s,  so the line is t = g s.

The origin itself is always inconsistent: there M = 0 while (R) demands 1.

---

## 5. Finiteness of the degenerate directions

The Euler count in §6 sums over degenerate lines, so their number must be
finite. This is not automatic and requires proof.

**Lemma 4 (finiteness).** Assume the D_r linearly independent and
m <= k - 1. Then only finitely many p in P^1 are degenerate.

*Proof.* If m = 0 there are no cuts and hence no degenerate directions, so
assume m >= 1 (otherwise the Wronskian below would be taken of an empty
collection). Work on an affine chart of P^1 and dehomogenize, writing
G_1(x), ..., G_m(x) for the corresponding polynomials. Suppose some nonzero
G = sum_r beta_r G_r has a zero of multiplicity at least k at x = a. Since
m <= k - 1, the first m derivatives of G vanish at a, so the columns of the
m x m Wronskian matrix (d^i G_r / d x^i)(a), 0 <= i <= m-1, admit the
nontrivial relation beta. Hence

        Wr(G_1, ..., G_m)(a) = 0.

In characteristic zero the Wronskian of linearly independent polynomials is
not identically zero, so it has finitely many roots. The point at infinity is
handled by repeating the argument in the other chart. QED

---

## 6. The count on a degenerate line, and the theorem

**The case k = 1, disposed of first.** The hypothesis n = 3 - m >= 3 forces
m = 0, so there are no cuts and hence no degenerate lines: Inc = {0}, and
chi_c(V) = 1 - 1 = 0 != 1. **We may therefore assume k >= 2 for the rest of
this section**, which is what makes sigma^{k-1} a nonconstant function of
sigma in Lemma 5.

**Lemma 5.** Let k >= 2 and let L_p be a degenerate line, parametrized by
sigma via (t, s) = sigma (t_0, s_0). The inconsistent part of L_p \ {0}
has

        chi_c = 0     or     chi_c = -(k-1).

*Proof.* Row 0 is homogeneous of degree k in (t, s); each cut row is
homogeneous of degree 1. Writing R_0(1), R_r(1) for the rows at sigma = 1, a
dependency beta_0 R_0 + sum_r beta_r R_r = 0 becomes, after dividing by
sigma,

        beta_0 sigma^{k-1} R_0(1) + sum_r beta_r R_r(1) = 0.

Let K = {beta' in C^m : sum_r beta'_r R_r(1) = 0}, independent of sigma.

*Case A: R_0(1) not in span{R_r(1)}.* Then beta_0 = 0 is forced and

        ker M(sigma)^T = {(0, beta') : beta' in K},

independent of sigma. Consistency is the condition that every left-kernel
vector annihilate the right-hand side, here beta' . v_cut = 0 for all
beta' in K — also independent of sigma. So the whole punctured line is
consistent, or the whole of it is inconsistent; in the first case the
inconsistent part is empty and in the second it is C^*. Either way
chi_c = 0.

*Case B: R_0(1) = sum_r gamma_r R_r(1).* Then

        ker M(sigma)^T = span{(1, -sigma^{k-1} gamma)} + {(0, beta') :
                              beta' in K},

(gamma is determined only modulo K, but gamma . v_cut is well defined once
condition (i) below holds). Consistency requires

    (i)  beta' . v_cut = 0 for every beta' in K      [independent of sigma],
    (ii) sigma^{k-1} (gamma . v_cut) = 1.

If (i) fails, or gamma . v_cut = 0, no point of the line is consistent: the
inconsistent part is C^*, chi_c = 0. Otherwise, since k >= 2 the exponent
k - 1 is at least 1, so (ii) holds at exactly k-1 nonzero values of sigma,
and the inconsistent part is C^* minus k-1 points, chi_c = -(k-1). QED

(For k = 1 the exponent vanishes and (ii) reads gamma . v_cut = 1, which is
independent of sigma: the inconsistent part is then empty or all of C^*, of
chi_c = 0 either way. That case was disposed of above and is not needed.)

**Proof of the theorem.** By Lemmas 3 and 4, Inc is the disjoint union of the
origin and the inconsistent parts of FINITELY many degenerate lines. By
additivity of chi_c and Lemma 5,

        chi_c(Inc) = 1 - (k-1) N,

with N the number of degenerate lines in Case B satisfying (i) with
gamma . v_cut != 0. By (*),

        chi_c(V) = (k - 1) N.

Suppose V ≅ C^n, so chi_c(V) = 1. Then (k-1) N = 1, forcing

        k = 2   and   N = 1.

With k = 2 the hypothesis n = k + 2 - m >= 3 gives m <= 1. If m = 0 there are
no cut rows, hence no degenerate lines, N = 0 and chi_c(V) = 0 — excluded. So
m = 1. By Lemma 2 the class is (3) or (2,1). For (3) the dual form is a cube,
so l_p^{k+1} divides it, L(d)G vanishes identically on the degenerate line,
and the line falls in Case A: N = 0. Hence the class is (2,1). For that
class the dual form has a unique root of multiplicity >= k = 2 — its double
root — so by Lemma 2 there is exactly one degenerate direction, and hence at
most one line to count. Finally, condition (ii) requires
gamma . v_cut = gamma_1 v_1 != 0, so the cut value
v_1 must be nonzero; when it is, (ii) has exactly k - 1 = 1 solution, giving
N = 1. QED

**Remark (k = 1).** Handled before Lemma 5: n = 3 - m >= 3 forces m = 0, so
Inc = {0} and chi_c(V) = 0. Note this is not a statement about "covering
degree 2": by §1 the normalized map is not generically finite when k = 1.

---

## 7. Moduli

The argument is uniform in all root positions and in the joint moduli of a
collection of cuts. What the root data affect is WHICH lines are degenerate
and the value of the integer N; what they do not affect is the divisibility

        chi_c(V) in (k-1) Z,

which is what forces k = 2. This is weaker than the claim in version 1 of
this note, which asserted that only root multiplicities matter. That is true
for a single cut, where Lemma 2 restricts the class to (k+1) or (k,1), but
false for m > 1, where Lemma 2 constrains the linear SPAN rather than the
individual classes.

---

## 8. A worked example: the span condition is not about individual classes

For each k there is a pair of classes, neither of which is (k,1) or (k+1),
whose span nevertheless degenerates:

        k = 3:  (2,2) + (2,1,1)
        k = 4:  (3,2) + (3,1,1)
        k = 5:  (4,2) + (4,1,1)

all of the form (k-1,2) + (k-1,1,1). Concretely, with the shared root of
multiplicity k-1 taken to be Z,

        G_1 = W^2 Z^{k-1},      G_2 = -W Z^{k-1} (W - Z),

neither divisible by Z^k, while

        G_1 + G_2 = W Z^k

is. So the pair degenerates at p = [Z]. At that direction R_1(1) + R_2(1) is
a nonzero scalar multiple of R_0(1): with the apolarity normalization of §3,
in which the functional Q |-> Q(alpha, beta) corresponds to
(1/k!)(alpha Z + beta W)^k rather than to the unscaled pure power, one has

        R_1(1) + R_2(1) = k! R_0(1)      (verified for k = 3, 4, 5),

so one may take gamma = (1/k!)(1, 1). By Lemma 5, Case B is productive
precisely when

        gamma . v_cut = (v_1 + v_2)/k! != 0,   i.e.   v_1 + v_2 != 0.

The scalar k! cancels out of the condition, so the classification below does
not depend on the normalization convention.

Hence

        chi_c(V) = k - 1   if v_1 + v_2 != 0,
        chi_c(V) = 0       if v_1 + v_2 = 0.

The value vectors enumerated in §10, namely (1,1) and (1,0), both have
nonzero sum, which is why the tables report k - 1 for these pairs. Checked
directly at k = 3, 4, 5 with the value vectors (1,1), (1,0), (3,1) — all
giving k - 1 — and (1,-1), (2,-2) — both giving 0.

The example therefore illustrates two things at once: that Lemma 2 constrains
the SPAN rather than the individual classes, and that the affine cut VALUES
enter the classification independently of the classes.

---

## 9. Revision history

### First review

Version 1 was reviewed adversarially; the following were corrected.

1. §1: the normalized covering degree is k^2 - 1, not C(j+k, j); and for
   k = 1 the map is not generically finite, so the earlier "covering degree
   two is unattainable" framing was wrong.
2. §2, §4: linear independence of the cut functionals added as a standing
   hypothesis. Lemma 2 is false without it.
3. §5: the finiteness of the degenerate directions, previously assumed, is
   now proved (Wronskian).
4. §2: the reduction to the plane restated with chi_c and rank strata; S and
   Inc are constructible, not locally closed.
5. §6: the conclusion now includes the condition v_1 != 0, without which
   N = 0 and chi_c(V) = 0. Also corrected: for class (k+1) with cut value
   zero the degenerate line is wholly CONSISTENT, not wholly inconsistent;
   either way N = 0.
6. §1: scope restricted to n >= 3 in the title, statement and text; n = 2 is
   not covered.
7. §7: the moduli remark weakened, as it was true only for a single cut.
8. §4: the root labels corrected. With d_g = d_z - g d_w, the operator d_z is
   g = 0 (line t = 0) and d_w is g = infinity (line s = 0); version 1 had
   these reversed.
9. §4, §6: the coordinate-change determinant and the left-kernel identities
   displayed explicitly, rather than left to be inferred.

---

### Second review

10. Lemma 5 assumed k >= 2 implicitly: for k = 1 condition (ii) is
    sigma-independent and does not have "exactly k-1 solutions". The case
    k = 1 is now disposed of before the lemma (n >= 3 forces m = 0), and the
    lemma states the hypothesis.
11. §8: the worked example depends on the cut VALUES as well as the span.
    At the degenerate direction R_0 = R_1 + R_2, so gamma = (1,1) and the
    pair is productive iff v_1 + v_2 != 0. Stated, with verification at
    value vectors of both kinds.
12. §1: the "étale and finite-to-one" sentence now carries the hypothesis
    k >= 2, since the same section later notes that k = 1 is not generically
    finite.
13. §2: the rank stratification refined so that a specified maximal minor is
    nonzero on each piece; a constant-rank locus need not admit a single
    global one.
14. §5: Lemma 4 now begins by excluding m = 0, which would otherwise take a
    Wronskian of an empty collection.

### Third review

15. §8: the scalar relating the rows at the degenerate direction is k!, not
    1, under the apolarity normalization of §3; gamma = (1/k!)(1,1) rather
    than (1,1). The productivity condition v_1 + v_2 != 0 is unchanged, since
    k! cancels, and no other statement is affected.
16. §6: made explicit that the class (2,1) has exactly one degenerate
    direction, its double root, so the final step counts one line.
17. Masthead and §9 relabelled: this is version 3, and the two earlier
    change lists are now subsections of a single revision history.

---

## 10. Verification

### 10.1 What supports the theorem

`verify_chi_tables.py` computes chi_c(V) independently of §§3–6, by rank
comparison over the (t,s)-plane, and asserts the outcome. For every
enumerated configuration with n >= 3 and linearly independent cuts:

* k = 3, m = 2: 20 configurations, chi_c in {0, 2};
* k = 4, m = 2 and 3: 108 configurations, chi_c in {0, 3};
* k = 5, m = 2: 110 configurations, chi_c in {0, 4};

all multiples of k-1, as §6 requires, with chi_c(V) = 1 occurring nowhere.
(The k = 5 row is pairs only; m = 3, 4 at k = 5 are not enumerated and are
not claimed.) Run with `--full` the script prints all 238 configurations, so
the aggregate assertions can be checked line by line rather than on trust.
The controls report chi_c(V) = 1 for (1,2) class (2,1) with nonzero value,
and `reconstruct.py` builds the counterexample explicitly from that
configuration, verifying a constant nonzero Jacobian determinant and a
tested fiber with three distinct points.

### 10.2 The clean-room runbook

`checkall.sh` runs the whole package in tiers, stopping at the first failure,
and writes `checkall.log`. Tiers must be read as having DIFFERENT evidential
weight; the runner labels each accordingly.

| tier | what runs | status |
|---|---|---|
| 0 | environment: python3, sympy, Singular, file inventory | preconditions |
| 1 | `miracle.py 1 2 1` -> `reconstruct.py`: the machinery must find and rebuild Alpöge's counterexample | positive control |
| 2 | `verify_chi_tables.py --full`, `moduli.py` controls, `euler.py`, `audit.py`, `multicut.py` | see below |
| 3 | `miracle.py` over covering degrees <= 10 | EXPLORATORY SCREENING |
| 4 | `sweep.sh` over j+k <= 6 (`--sweep`, hours) | EXPLORATORY SCREENING |

Within tier 2 the components differ in kind:

* **`verify_chi_tables.py`** — the verification of this paper's theorem, as
  in §10.1. Exact over Q; no modulus, no hypothesis.
* **`euler.py`** — the k = 2 Euler computation done by hand, an independent
  cross-check of `chi_plane.py` on the classes (1,1,1) and (2,1).
* **`moduli.py`** — symbolic cross-ratio analysis for single-cut classes,
  confirming that the single-cut classification is uniform in the moduli.
* **`audit.py`** — global-unit certificates. These prove a DIFFERENT
  statement from the theorem: that certain V carry a nonconstant unit and so
  are not affine space of any dimension. Independent of the Euler argument.
* **`multicut.py`** — **DIAGNOSTIC ONLY.** It checks a left-kernel
  classification of the degenerate fiber V ∩ {t = 0} against Gröbner
  computation (92/92 agreement). It does NOT bear on whether V itself is
  affine space: an earlier version of this project inferred global
  non-affineness from the fiber, which is invalid (see §11), and this script
  is retained only as a consistency check on the fiber-level bookkeeping.

### 10.3 What the screening tiers do and do not establish

Tiers 3 and 4 use `miracle.py`, whose test asks whether the lexicographic
Gröbner basis of a degenerate fiber is triangular. That is a SUFFICIENT
certificate that the fiber is affine space, not an equivalence: x^2 - y = 0
defines A^1 yet presents as "branching" when x leads. Consequently a null
result from these tiers is screening — evidence that no candidate surfaced —
and NOT rigorous exclusion. They are also mod p (default 32003).

Their standing value is for j >= 2, which is outside the reach of §§3–6: the
apolarity argument uses deg L = 1 twice (the resultant is an evaluation, and
exactly one derivative falls on L), and for j >= 2 the resultant is no longer
linear in c_0. The recorded sweep covers 32 cells with 0 candidates, 27 cells
completed with no unresolved combinations, and 4 cells incomplete at the
per-cell time cap — of which (3,3) with 3 cuts is separately known to
complete cleanly (330 combinations, 0 unresolved) from a longer run of
`resume.sh` — a manual follow-up, not part of `checkall.sh` — leaving three
cells at covering degree 15 with 31 unresolved combinations, all enumerated
by name in `sweep_results.txt`. None of this is used by the theorem.

### 10.4 Files

Proof-bearing: `verify_chi_tables.py`, `chi_plane.py`, `euler.py`,
`reconstruct.py`, `audit.py`, `moduli.py`, and `trichotomy.py` (the
self-contained k = 2 companion). Diagnostic: `multicut.py`, `hcheck.py` (an
earlier chi computation requiring a fiber-constancy hypothesis, superseded by
`chi_plane.py` and retained because the two agree on every shared case).
Screening: `miracle.py`, `sweep.sh`, `resume.sh`. Shared helpers:
`symmult.py`. Runner: `checkall.sh`.

Logs, and which run produces each:

* `checkall.log` — written by `checkall.sh` on every run (tiers 0–3).
* `sweep_results.txt` and `sweep_miracles.txt` — written by `sweep.sh`, i.e.
  only when `checkall.sh` is invoked with `--sweep` (tier 4).
* `resume.log` — written by `resume.sh`, which `checkall.sh` never calls.
  It is a manual follow-up for cells that hit the per-cell time cap during a
  sweep; it exists only if that follow-up was run. The (3,3) three-cut result
  cited in §10.3 comes from such a run.

Singular is required; a sympy fallback previously advertised in `miracle.py`
was removed rather than repaired, after it was found to fail on the positive
control.

---

## 11. Provenance

Framework: T. Tao, "A digestion of the Jacobian conjecture counterexample"
(July 2026), on the counterexample of L. Alpöge with Claude Fable 5.

Developed in adversarial human-AI collaboration. Two intermediate claims did
not survive scrutiny and are recorded because they shaped the final proof: a
deformation direction that proved to be a formal jet which does not
integrate, and — more consequentially — an assertion that V ≅ C^n implies the
fiber V ∩ {t = 0} is affine space. The latter is false in general (a
hypersurface in C^n can be C^* x C^{n-2}), and its collapse is what forced
the replacement of fiber arguments by the global invariants — units, then
Euler characteristic — used here and in `TRICHOTOMY.md`.

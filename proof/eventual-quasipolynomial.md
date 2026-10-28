# Eventual quasipolynomial complexity for every star

**Review status:** the finite-description, representative and counting argument was independently reviewed in a separate context on 24 October 2026. It does not depend on the supplementary explicit-coefficient construction. The general counting tool is an existing theorem; its application here is the inclusion of all star components and the identification of their common leading coefficient.

## 1. Statement and counting theorem

Let $`\theta`$ be any star as defined in [objects](objects.md), and let $`P`$ be a fixed positive-area rational convex polygon. Write $`S_n=nP\cap\mathbb Z^2`$ and let $`\Lambda(P)`$ be the direction-band union area from that page.

**Eventual quasipolynomial theorem.** There are positive integers $`M,N`$ and rational constants $`\beta_r,\gamma_r`$, for $`0\le r\lt M`$, such that

```math
P_\theta(S_n)=\Lambda(P)n^2+\beta_{n\bmod M}n+\gamma_{n\bmod M}
\qquad(n\ge N).
```

There is no restriction to Case B, nonzero line support, equal tails, a unimodular pair of directions, or a zero periodic background. The modulus can depend on the entire finite star specification. This theorem alone does not give its value or assert $`M=1`$ for integral $`P`$.

The counting input is Bogart–Goodrick–Woods, *Parametric Presburger Arithmetic: Logic, Combinatorics, and Quasi-polynomial Behavior*, arXiv:1608.08520v2, Theorem 1.15, printed p. 6, and Property 2, p. 5. For a parametric Presburger family of subsets of $`\mathbb Z^d`$, it supplies an eventual quasipolynomial agreeing with the cardinality on every finite fiber. Its Definition 1.3, p. 2, allows fixed formulas built from integer linear operations, Boolean operations, integer quantifiers and multiplication by the distinguished parameter. It does not allow multiplication of two other variables or quantification of that parameter. Definitions 3.1 and 3.3, p. 14, give the formal language. The extra nonnegative-coordinate assumption in Theorem 1.15 concerns its generating-function property, not the counting property used here.

For fixed rational polygon dilations, the formula below is even ordinary Presburger arithmetic when $`n`$ is treated as a free variable. Older fixed-normal counting results, such as Woods's Theorem 1.25 and Theorem 3.5(b) in arXiv:1308.4694v2, suffice for this special case. We use the later theorem for a uniform citation; we do not claim a new Presburger counting theorem.

## 2. Finite descriptions of all star colors

For every doubly periodic tail or background $`G`$, its period lattice has finite index, so it contains $`q\mathbb Z^2`$ for a positive integer $`q`$. A common multiple works for all these fields. Their colors are therefore specified by the finite tables $`G(b)`$, $`b\in\{0,\ldots,q-1\}^2`$. For example the predicate for color $`a`$ is

```math
\bigvee_{\substack{b\in\{0,\ldots,q-1\}^2\\G(b)=a}}
 (z_1\equiv b_1\pmod q\ \wedge\ z_2\equiv b_2\pmod q).
```

A congruence means $`\exists h\ (z_j-b_j=qh)`$, so it uses only a fixed integer coefficient.

Complete each primitive $`v_i`$ to an integer basis $`(v_i,e_i)`$ of determinant one. Define

```math
t_i(z)=\det(v_i,z),\qquad s_i(z)=\det(z,e_i),
\qquad z=s_i(z)v_i+t_i(z)e_i.
```

For $`\ell_i\le t\le r_i`$ and $`0\le b\lt k_i`$, record $`c_i(t,b)=F_i(bv_i+te_i)`$. The period $`k_iv_i`$ shows that this determines the whole transition strip. Define $`\Phi_{i,a}(z)`$ as the disjunction of: the left-tail color predicate with $`t_i(z)<\ell_i`$; the right-tail predicate with $`t_i(z)>r_i`$; and the finitely many conditions $`t_i(z)=t,\ s_i(z)\equiv b\pmod{k_i}`$ for which $`c_i(t,b)=a`$.

These regions partition the integer plane, including the empty-strip case $`\ell_i=r_i+1`$. Thus $`\Phi_{i,a}`$ means exactly $`F_i(z)=a`$. A finite disjunction over color tuples gives

```math
\Theta_a(z)=
\bigvee_{\substack{(a_1,\ldots,a_m)\in\mathbb F_p^m\\\sum_i a_i=a}}
 \bigwedge_i\Phi_{i,a_i}(z).
```

The finite-field sum is a finite truth table, not multiplication or quantification over a second infinite arithmetic domain. It includes all cancellations. A periodic background can be included as one more finite table, or absorbed into a component as in CN Remark 0.1. Exactly one $`\Theta_a`$ holds at every point. No characteristic-zero annihilator or Case B assumption was used.

## 3. A finite definable transversal of the patterns

After clearing positive denominators, write $`P=\{x:A_jx\le b_j\ (1\le j\le J)\}`$, with integral $`A_j,b_j`$. For $`n\ge1`$, let $`\mathrm{Win}(n,s)`$ be the conjunction $`A_js\le b_jn`$. At $`n=0`$ use the singleton window $`\{0\}`$, an irrelevant finite prefix.

Pattern equality is the fixed formula

```math
E_n(u,v)\ \Longleftrightarrow\ 
\forall s\in\mathbb Z^2\ 
\left[\mathrm{Win}(n,s)\Longrightarrow
\bigvee_a(\Theta_a(u+s)\wedge\Theta_a(v+s))\right].
```

Outside the window the implication is vacuous. No list of variables of length depending on $`n`$ is used, and no pattern is quantified as a function. Universal site quantifiers are legal.

Order anchors by

```math
u\prec v
\quad\Longleftrightarrow\quad
(\|u\|_1,u_1,u_2)<_{\mathrm{lex}}(\|v\|_1,v_1,v_2).
```

Absolute value is definable by its two sign cases; the order is a fixed Presburger formula. A nonempty subset of the integer plane has a least attained nonnegative norm, and that norm shell is finite. Consequently every pattern class, bounded or not, has exactly one least anchor. Define

```math
R_n=\{u:\neg\exists v\,[E_n(u,v)\wedge v\prec u]\}.
```

The map from $`R_n`$ to the realized patterns is a bijection. Since the window and alphabet are finite, $`R_n`$ is finite. This does not assume a guessed enclosing box or a polynomial bound on representative coordinates. The empty-window pattern, if encountered, has representative zero.

The family $`R_n`$ therefore satisfies the quoted counting theorem. There are $`M,N`$ and rational polynomials $`q_r`$ with $`P_\theta(S_n)=q_{n\bmod M}(n)`$ for $`n\ge N`$. The same construction proves the broader principle: a fixed Presburger coloring on $`\mathbb Z^d`$, on any one-parameter Presburger family of finite windows, has eventual quasipolynomial translated-pattern counts. That broader principle alone imposes no quadratic degree.

## 4. Identification of the degree and coefficients

The [all-star area theorem](all-star-extension.md) applies to this same fixed star and polygon:

```math
P_\theta(S_n)=\Lambda(P)n^2+O(n).
```

On any one infinite residue progression modulo $`M`$, the polynomial $`q_r(n)-\Lambda(P)n^2`$ is $`O(n)`$. A nonzero coefficient of degree greater than one would violate that bound. Thus every $`q_r`$ is $`\Lambda(P)n^2+\beta_rn+\gamma_r`$. The band intersections are rational polygons, so $`\Lambda(P)`$ is rational; it is positive because $`P`$ lies in every band. This proves the theorem, including a common quadratic coefficient in every residue.

## 5. Compatibility with the direct line coefficients

The [full-support periodic-line theorem](full-support-lines.md) gives a periodic linear term with bounded remainder using a vertex-denominator modulus $`M_0`$. Refine its modulus and the preceding $`M`$ to their least common multiple. On each refined residue, a different linear coefficient would produce an unbounded difference. Hence the coefficients agree and the bounded remainder is eventually periodic.

In particular the binary lines $`y-2ix=i`$, $`i=0,1,2`$, on the standard triangle satisfy

```math
P_\theta(S_n)=\frac{31}{4}n^2+\frac{11}{4}n+\gamma_{n\bmod M}
```

eventually. No numerical value of the full period or constant table is asserted from a finite data fit.

For the full $`\mathbb F_3`$ cross on the half-width rectangle, its exact count is

```math
(\lfloor n/2\rfloor+2)(n+2)
=\frac12n^2+\left(3-\frac r2\right)n+4-r,\qquad r=n\bmod2.
```

The sparse binary cross has instead

```math
\frac12n^2+\left(\frac92-\frac r2\right)n+6-\frac{3r}{2}
\qquad(n\ge8).
```

Their complete representation proofs are in the full-support module, Section 2; the sparse example is a counterexample rather than a member of that module's full-support positive class.

## 6. Why rationality cannot simply be removed

Keep the full $`\mathbb F_3`$ cross and take $`P=[0,\alpha]\times[0,1]`$ with positive irrational $`\alpha`$. For sufficiently large $`n`$, the window has at least two rows and columns. A pattern is determined uniquely by an absent-or-present row and an absent-or-present column: a two-line pattern has a unique color-two crossing; pure row, pure column and zero patterns are distinct. All choices are realized by translations. Therefore

```math
P_\theta(S_n)=(\lfloor\alpha n\rfloor+2)(n+2)
=\alpha n^2+(2\alpha+2-\{\alpha n\})n+4-2\{\alpha n\}.
```

If a periodic linear coefficient left a bounded remainder, on $`n=Mk+r`$ the numbers $`\{\alpha(Mk+r)\}`$ would converge. Their successive differences are either $`\{M\alpha\}`$ or $`\{M\alpha\}-1`$, both bounded away from zero because $`M\alpha`$ is irrational. A convergent sequence has successive differences tending to zero, a contradiction. Thus even this full-support example has no periodic linear coefficient with bounded remainder on this irrational window.

## 7. Verification and status

[The finite-description diagnostics](../scripts/presburger_examples_audit.py) compare direct fields against finite color tables on 2,401 points, including nonunimodular directions, unequal nonconstant tails, zero entries, an empty transition strip and a nonzero periodic background. They also globally enumerate the specified line examples by all multiple-hit anchors, all realizable pure phases and zero. They do not infer a period from finite values.

The [review-status record](../sources/review-status.md) identifies the source and independent review. The separate [explicit-coefficient supplement](explicit-coefficients.md) passed one fresh-context GPT-6.1 Sol adversarial review on 25 October 2026. None of the arguments in Sections 1–5 relies on that supplement.

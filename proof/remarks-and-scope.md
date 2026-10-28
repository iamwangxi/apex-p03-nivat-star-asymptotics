# Numerical comparisons, attribution and scope

## 1. An explicit Case B bound on the three-line triangle

For the binary three-line example in [Theorem 3](three-line-second-order.md), the global strip profiles are, by transverse height,

```math
f_0:\{-6:1,-4:-1,-2:-1,0:1\},\qquad
f_1:\{-1:-1,1:2,3:-1\},\qquad
f_2:\{2:1,4:-1,6:-1,8:1\}.
```

Apply the [Case B anchor lemma](signals-and-anchors.md) to $`S_{30}`$. Its three difference erosions are $`S_{22},S_{24},S_{26}`$. Since all tangential periods are one, their projected row sets give

```math
C_0=[-22,-6]\cap\mathbb Z,
```

```math
C_1=([-21,45]\cap\mathbb Z)\cup\{47\},\qquad
C_2=([-18,94]\cap\mathbb Z)\cup\{96,97,98,102\}.
```

For example the projections of the last two erosions are

```math
\pi_1(S_{24})=\{-48\}\cup([-46,24]\cap\mathbb Z),
```

```math
\pi_2(S_{26})=
\{-104,-100,-99,-96,-95,-94\}\cup([-92,26]\cap\mathbb Z).
```

The three pairwise good-anchor counts, with integer congruences retained, are 577, 499 and 3975. The congruences are respectively $`t_1\equiv y\pmod2`$, $`t_2\equiv y\pmod4`$ and $`t_1\equiv t_2\pmod2`$. The triple count is 496; hence

```math
|\mathcal A(S_{30})|=577+499+3975-2\cdot496=4059.
```

On this same 496-site window, CN Theorem T (p. 4; Theorem 7.3, p. 18) gives 497. Exact global enumeration gives 7037. The first coefficient is $`31/4`$, whereas CN's $`|S_n|+1`$ gives leading lower-bound term $`n^2/2`$.

For the Case A quadrant star, [Theorem 2's example](finite-window-bound.md#2-a-case-a-numerical-example) gives 786 on that same $`S_{30}`$, compared with CN's 497 and the exact global count 1022. These are two different configurations on one window; their pattern sets are not combined. Neither new finite bound is asserted to dominate CN on every small window.

The complete global enumeration is justified in the next two subsections. The scripts use exact integer masks to decide pattern equality; no finite sampling box is taken as evidence of global coverage. See also [the reproducibility instructions](../README.md#reproducing-the-computations).

### 1.1. Complete line enumeration

For the triangle $`S_n`$, a nonempty section in direction $`i\in\{0,1,2\}`$ at an anchor $`(X,Y)`$ has parameter

```math
t_i=i-Y+2iX\in[-2in,n].
```

The intervals are enclosing ranges; the missing layers listed in [the section calculation](three-line-second-order.md#1-sections-and-representations) do not enlarge them. Solving the two parameter equations for each pair gives the following bounds. For hits in directions 0 and 1,

```math
Y=-t_0\in[-n,0],\qquad
X=(t_1-t_0-1)/2\in[-(3n+1)/2,(n-1)/2].
```

For hits in directions 0 and 2,

```math
Y=-t_0\in[-n,0],\qquad
X=(t_2-t_0-2)/4\in[-(5n+2)/4,(n-2)/4].
```

For hits in directions 1 and 2,

```math
X=(t_2-t_1-1)/2\in[-(5n+1)/2,(3n-1)/2],\qquad
Y=t_2-2t_1\in[-6n,5n].
```

Thus every multiple-hit integer anchor lies in

```math
[-3n-2,2n+2]\times[-6n-2,5n+2],
```

exactly the inclusive rectangle scanned by `three_line_enumeration.py` and `area_and_anchor_audit.py`. The auxiliary concurrent-line check has $`t_i=-Y+2iX`$: removing the constants 1 and 2 in the inverse formulas gives the same enclosing rectangle. The finite good-anchor count is covered too. In this disjoint-line example the integer encoding is the sum of the three line indicators; $`Q_i`$ kills the other two components. A nonzero observed $`Q_i`$ value therefore requires the line in direction $`i`$ to hit the window. Two good directions imply two line hits.

Every pure section $`Q_i(t)`$ is realized at the concrete anchor

```math
(X,Y)=(L,i+2iL-t),\qquad L=20(n+1).
```

At this anchor the other parameter is $`t_j=t+(j-i)(2L+1)`$. Since $`-4n\le t\le n`$, it lies above $`n`$ for $`j>i`$ and below $`-2jn`$ for $`j\lt i`$. Thus exactly the chosen section is present. In the concurrent check use $`Y=2iL-t`$; then $`t_j=t+2(j-i)L`$ gives the same conclusion. The anchor $`(L,-10L)`$ realizes zero in both checks, since all its parameters exceed $`n`$.

Every anchor outside the scanned rectangle has at most one hit, so it contributes only a pure section or zero already in the representation list. Deduplicating that complete list gives every global pattern, including all collisions. This proves the coverage used for the count 7037 and for each enumerated scale $`3\le n\le164`$; it makes no claim that those finite checks prove the general asymptotic theorem.

### 1.2. Complete quadrant enumeration

For the quadrant star $`\theta(x,y)=\mathbf1[y\ge0]+\mathbf1[x\ge0]\pmod2`$ and a site $`s=(x_s,y_s)\in S`$, its value at the anchor $`(X,Y)`$ is

```math
\theta(X+x_s,Y+y_s)
=\mathbf1[Y\ge-y_s]+\mathbf1[X\ge-x_s]\pmod2.
```

It changes only when $`X`$ crosses a cut $`-x_s`$ or $`Y`$ crosses a cut $`-y_s`$. If the distinct cuts in one coordinate are $`c_1<\cdots<c_r`$, the disjoint integer intervals are

```math
(-\infty,c_1-1]\cap\mathbb Z,\quad
[c_j,c_{j+1}-1]\cap\mathbb Z\ (1\le j\lt r),\quad
[c_r,\infty)\cap\mathbb Z.
```

Their union is every integer in that coordinate, including both unbounded ends. The product cells therefore partition every integer anchor in the plane, and every site value is constant on each cell. One representative per cell gives every global pattern, even for a window with holes. This is exactly `threshold_cells` followed by `quadrant` in `finite_window_audit.py`.

The program retains all cells giving each mask, not just their representatives. A pattern occurs at exactly one integer anchor if and only if its cells consist of a single product of singleton intervals. Any unbounded cell contains infinitely many anchors; any other nonsingleton cell or second cell gives at least two. The counted central and distant anchors are checked against precisely one singleton cell each. Hence the count 1022 and the uniqueness checks are global, rather than confined to a sampling box.

## 2. Remark: classical threshold sharpness

The elementary single-defect configuration $`\xi=\mathbf1[\{0\}]`$ is nonperiodic and satisfies, for every nonempty finite $`S`$,

```math
P_\xi(S)=|S|+1.
```

Indeed its patterns are zero and a single one at each possible site of $`S`$. Consequently, for a positive-area lattice polygon $`P`$ and any $`a>\mathrm{area}(P)`$,

```math
P_\xi(nP\cap\mathbb Z^2)\lt an^2
```

for every sufficiently large $`n`$, because lattice counting gives $`|nP\cap\mathbb Z^2|=\mathrm{area}(P)n^2+O_P(n)`$.

This is a direct observation expressing classical sharpness, not a new theorem. Szabados, *Nivat's conjecture holds for sums of two periodic configurations*, arXiv:1710.05360v1, Introduction, p. 2, explicitly recalls nonperiodic configurations with rectangular complexity $`MN+1`$ and cites Cassaigne's 1999 classification. Cassaigne's original classification was not directly checked for this package; the historical attribution is through Szabados.

The larger star coefficient cannot therefore be substituted as a periodicity threshold for arbitrary configurations. The single defect is outside the star class: in $`\mathbb F_p`$, a nonzero Laurent polynomial acting on the single defect cannot vanish, while a star is annihilated by $`D`$ in that field. The proof above does not rely on CN's general reduction.

## 3. Contribution and dependence

CN supplies the star hypotheses and finite-difference geometry (§0, pp. 4–5; Lemma 1.1, pp. 5–6), the Case B nonzero-strip argument (Lemma 4.1, p. 13), and the product of two shifted strip signals used to make a finite marker (Lemma 4.2, pp. 13–14). The construction of that product is inherited. Here it is combined with common transverse extrema at both ends, bounded displacement of any equal-pattern anchor, and observation of the complete marker support. Together with fixed erosion losses and matching active-band counts, these give the exact coefficient for all stars and arbitrary fixed convex shapes. The finite-window theorem makes the core correction explicit. The example's linear term comes from a complete representation and collision count.

Quadratic complexity growth by itself is not claimed as new. Kari–Szabados, *An Algebraic Geometric Approach to Nivat's Conjecture*, arXiv:1605.05929v1, Theorem 44, p. 30, proves $`P_c(M,N)>MN`$ for every nonperiodic configuration outside finitely many rectangle sizes. Its supporting Lemmas 38–40 (pp. 25–27) and Lemma 42 (pp. 28–29) already yield quadratic increments in suitable rectangular regimes. The present quantitative result identifies the exact star coefficient for arbitrary fixed convex shapes and supplies a matching upper bound. This comparison is with the specified sources; no comprehensive literature search or priority certification has been carried out.

The proofs do not use the Kari–Szabados ideal-structure theorem, a spectral polynomial, or a scalar ideal bridge. This is a dependency statement about these quantitative star results. It does not remove those ingredients from CN's separate proof of a lower bound for every individual convex window.

The reviewed [eventual-quasipolynomial theorem](eventual-quasipolynomial.md) and [full-support line formula](full-support-lines.md) extend the asymptotic information. The [general coefficient supplement](explicit-coefficients.md) passed one fresh-context GPT-6.1 Sol adversarial review on 25 October 2026. The [second supplement](minimal-periods-and-constant-term.md), on exact constant terms and minimal periods, passed one on 26 October 2026. No compact original-boundary formula for every star, uniform configuration-independent boundary constant, or higher-dimensional area theorem is claimed. Computations can check instances and expose errors; they do not replace the infinite-range proofs. No human expert has certified this work.

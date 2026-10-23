# Leading term: geometry and the matching upper bound

Prerequisites: [objects](objects.md) and [the Case B anchor lemma](signals-and-anchors.md). The all-star lower bound is completed in [the extension](all-star-extension.md).

## 1. Uniform interior rows

We give the geometric boundary estimate needed for the general asymptotic statement. No smoothness, rationality or polygon assumption is hidden here.

**Lemma.** Let $`v`$ be primitive, $`\pi =\det(v,\cdot )`$, let $`F\subseteq \mathbb Z^2`$ be finite and nonempty, and fix $`k\ge 1`$. If $`P`$ is compact, convex and has nonempty interior, there exist constants $`C,N`$ such that, for every integer $`n\ge N`$ and every integer
```math
n\min_P\pi+C\le t\le n\max_P\pi-C,
\qquad\text{(7)}
```
the section at height $`t`$ of
```math
\{z\in\mathbb Z^2:z+F\subseteq nP\}
```
contains $`k`$ consecutive points in direction $`v`$.

**Proof.** Complete $`v`$ to a lattice basis $`(v,e)`$ with $`\det(v,e)=1`$. Measure chord lengths in units of $`v`$. Write $`l=\min_P\pi`$, $`r=\max_P\pi`$, choose $`q\in \mathrm{int} P`$, and let $`a=\pi (q)`$, so $`l\lt a\lt r`$. The chord-length function $`h(t)`$ is concave: convex combinations of two chord segments are contained in the intermediate chord, and their lengths combine linearly. With $`h_0=h(a)>0`$, concavity and nonnegative endpoint lengths give
```math
h(t)\ge c\min(t-l,r-t),\qquad
c=\frac{h_0}{\max(a-l,r-a)}>0.
\qquad\text{(8)}
```

Put $`G=F+\{0,v,\ldots,(k-1)v\}`$. Choose a fixed $`\lambda >0`$ so large that $`q+g/\lambda \in P`$ for every $`g\in G`$. For $`n>\lambda`$, convexity gives
```math
K_n:=(n-\lambda)P+\lambda q
\ \subseteq\ \bigcap_{g\in G}(nP-g).
\qquad\text{(9)}
```
Indeed, for $`x=(n-\lambda )y+\lambda q`$, $`x+g=(n-\lambda )y+\lambda (q+g/\lambda )\in nP`$.

The projection endpoints of $`K_n`$ are
```math
l_n=nl+\lambda(a-l),\qquad
r_n=nr-\lambda(r-a).
```
Scaling (8) gives chord length at least $`c\cdot \min(t-l_n,r_n-t)`$ on $`K_n`$. Choose
```math
C\ge\max\{\lambda(a-l),\lambda(r-a)\}+1/c.
```
Every height in (7) then has a chord in $`K_n`$ of length at least one in $`v`$ units. An integer-height line contains a lattice point because $`v`$ is primitive; its lattice points are spaced one $`v`$ unit apart. A closed segment of length at least one therefore contains a lattice point $`z`$. By (9), $`z+G\subseteq nP`$, so $`z,z+v,\ldots ,z+(k-1)v`$ belong to the erosion by $`F`$. This proves the lemma, with any integer $`N>\lambda`$. ∎

Apply the lemma for each $`i`$, using $`F=\mathrm{supp} Q_i`$ and $`k=\kappa _i`$. There is a fixed $`C_i^0`$ such that $`L_i(S_n)`$ contains every integer between $`np_i^-+C_i^0`$ and $`np_i^+-C_i^0`$. Hence $`C_i(S_n)`$ contains all integers in
```math
I_i^-(n)=
[\mu_i^+-np_i^++C_i^0,\ \mu_i^- -np_i^- -C_i^0].
\qquad\text{(10)}
```
The losses at both ends are constants, not a fixed proportion of the window.

## 2. Counting band intersections

We use an elementary lattice estimate. Fix pairwise independent-in-pairs linear forms $`\pi _i`$, real numbers $`a_i\lt b_i`$, and fixed endpoint shifts $`c_i,d_i`$. For the bands
```math
N_i(n)=\{x:na_i+c_i\le\pi_i(x)\le nb_i+d_i\},
```
let $`N(n)=\bigcup _{i\lt j}(N_i(n)\cap N_j(n))`$, and let $`N`$ be the analogous union with unscaled unshifted endpoints. Then
```math
|N(n)\cap\mathbb Z^2|=n^2\mathrm{area}(N)+O(n).
\qquad\text{(11)}
```

For completeness, each pair parallelogram has vertices $`n`$ times the corresponding fixed vertices plus $`O(1)`$: solve the same two independent linear equations at its four corners. Its symmetric difference from the unshifted dilate lies in a bounded-width neighborhood of a polygonal boundary of length $`O(n)`$, inside a ball of radius $`O(n)`$. That neighborhood has area and lattice-point count $`O(n)`$, by covering it with unit lattice squares. Summing over finitely many pairs gives the same estimate for the union. Finally, the difference between area and lattice-point count for $`nN`$ is $`O(n)`$: only unit squares intersecting the finitely many boundary segments can contribute a discrepancy, and there are $`O(n)`$ of them. This proves (11), including real, nonrational endpoints. Constants can depend on the fixed directions and shifts.

The bands in (10) have unshifted limits $`-B_i(P)`$. By the Case B anchor lemma and (11),
```math
P_\theta(S_n)\ge \Lambda(P)n^2-O(n).
\qquad\text{(12)}
```

## 3. The matching upper bound

This part uses only the star hypotheses, not Case B. For an anchor $`u`$, call direction $`i`$ active if
```math
\pi_i(u)\in I_i^+(n):=[\ell_i-np_i^+,\ r_i-np_i^-].
\qquad\text{(13)}
```
If it is inactive, all of $`u+S_n`$ lies in a single tail half-plane of component $`i`$: below the interval in (13) every point has $`\pi _i<\ell _i`$, and above it every point has $`\pi _i>r_i`$. This definition may overcount genuinely active components, which is harmless.

Anchors with at least two active directions belong to the union of pairwise intersections of the bands (13). Formula (11) bounds their number by
```math
\Lambda(P)n^2+O(n).
\qquad\text{(14)}
```
Their number is an upper bound on the number of patterns contributed by these anchors.

For exactly one active direction $`i`$, fix $`t=\pi _i(u)`$ and fix the left/right tail choices of the other components. There are $`O(n)`$ possible integers $`t`$, and at most $`2^{m-1}`$ tail choices. Two anchors at the same height differ by an integer multiple of primitive $`v_i`$. Modulo $`\kappa _i`$ such multiples, their patterns agree, since $`H_i`$ preserves $`F_i`$ and all the chosen tails. There are therefore at most $`\kappa _i`$ patterns for each fixed height and set of choices. Summing over $`i`$ gives $`O(n)`$ patterns.

With no active direction, the pattern comes from a sum of pure tails. The common period lattice
```math
\Gamma=\bigcap_i\bigl(\mathrm{Per}(L_i)\cap\mathrm{Per}(R_i)\bigr)
```
has finite index. At most $`2^m[\mathbb Z^2:\Gamma ]`$ patterns arise. Counting unrealized tail choices merely overcounts; no identity between their backgrounds is asserted.

Combining these cases with (14),
```math
P_\theta(S_n)\le\Lambda(P)n^2+O(n).
\qquad\text{(15)}
```
Equations (12) and (15) prove the area formula under Case B. All constants are for the same fixed configuration, periods, directions and shape. No conclusion for Case A is inferred from the Case B lower bound. ∎


## 4. Evaluating the area coefficient

Put $`W_i=p_i^+-p_i^-`$. The linear map $`x\mapsto(\pi_i(x),\pi_j(x))`$ has absolute determinant $`|\det(v_i,v_j)|`$ and sends $`B_i(P)\cap B_j(P)`$ onto the rectangle $`[p_i^-,p_i^+]\times[p_j^-,p_j^+]`$. The linear change-of-variables formula therefore gives

```math
\mathrm{area}(B_i(P)\cap B_j(P))
=\frac{W_iW_j}{|\det(v_i,v_j)|}.
\qquad\text{(16)}
```

More generally, with $`[m]=\{1,\ldots,m\}`$,

```math
\Lambda(P)=
\sum_{\substack{I\subseteq[m]\\ |I|\ge2}}
(-1)^{|I|}(|I|-1)
\mathrm{area}\left(\bigcap_{i\in I}B_i(P)\right).
\qquad\text{(17)}
```

Indeed, at a point belonging to exactly $`k`$ bands, the corresponding sum of indicators is zero when $`k=0,1`$. For $`k\ge2`$, the binomial identities give

```math
\sum_{r=0}^k(-1)^r\binom{k}{r}=0,\qquad
\sum_{r=0}^k(-1)^r r\binom{k}{r}=0,
```

and hence

```math
\sum_{r=2}^k(-1)^r(r-1)\binom{k}{r}=1.
```

This is the indicator of membership in at least two bands. Every intersection in (17) contains a pair of nonparallel bands and is bounded; integrating the finite indicator identity proves (17).

If $`P`$ is a lattice polygon, each $`p_i^\pm`$ is an integer, because the integer linear form $`\pi_i`$ attains its extrema at lattice vertices. All band boundary lines therefore have integer coefficients. Each vertex of an intersection polygon solves two independent such equations, so its coordinates are rational. The area is rational by the polygon area formula. The boundary vertices of the finite union also come from these line intersections, and its area is rational by (17). Thus the general coefficient can be computed by finite exact rational operations in this case. No rationality assumption is needed for (16) or (17) themselves.

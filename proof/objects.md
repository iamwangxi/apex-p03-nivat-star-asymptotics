# Objects and main results

Throughout, CN denotes the 12 September 2026 contest version of *The Convex Nivat Conjecture: A Complexity Lower Bound for Star Configurations, and a Reduction from Low Convex Complexity to Star Configurations*, Apex Intelligence. Its star hypotheses are in §0.1, p. 4. The exact source fingerprint and checked locators are in [the bibliography](../sources/bibliography.md).

## 1. Stars, encodings and observed rows

Let $`p`$ be a prime and $`m\ge2`$. Let $`v_1,\ldots,v_m\in\mathbb Z^2`$ be primitive, pairwise nonparallel vectors, and set $`\pi_i(z)=\det(v_i,z)`$. A CN star is

```math
\theta=\sum_{i=1}^m F_i:\mathbb Z^2\longrightarrow\mathbb F_p,
```

where $`F_i`$ has period $`k_iv_i`$, is not doubly periodic, and has doubly periodic tails $`L_i,R_i`$ satisfying

```math
F_i=L_i\text{ on }\{\pi_i<\ell_i\},\qquad
F_i=R_i\text{ on }\{\pi_i>r_i\},\qquad \ell_i\le r_i+1.
```

The transition strip may be empty. The tails need not be equal or constant, and the tangential periods need not be one. Doubly periodic means that the period group contains two independent vectors. A doubly periodic background can be absorbed into one component, as in CN Remark 0.1, p. 4.

Choose $`H_i=\kappa_i v_i`$ to preserve $`F_i`$ and every tail field; common multiples exist by the finite index of the tail period lattices (CN §0.3, p. 5). Put

```math
(T^u f)(z)=f(z+u),\qquad
D=\prod_{i=1}^m(T^{H_i}-1),\qquad
Q_i=\prod_{j\ne i}(T^{H_j}-1).
```

Choose any injective integer encoding $`w:\mathbb F_p\to\mathbb Z`$, and put $`\eta=w(\theta)`$, $`f_i=Q_i\eta`$. Only the original component sum is in $`\mathbb F_p`$; differences of encoded fields take place in characteristic zero. No additivity of $`w`$ is assumed.

Case B means the full color-indicator condition

```math
D\mathbf1[\theta=a]=0\quad(a\in\mathbb F_p).\qquad\text{(B)}
```

It implies $`D\eta=0`$. Case A means that at least one of these indicator differences is nonzero; a nonzero $`D\eta`$ certifies Case A. The main results cover both cases.

For a nonempty finite $`S\subset\mathbb Z^2`$, its pattern complexity is

```math
P_\theta(S)=\#\{(s\mapsto\theta(u+s)):u\in\mathbb Z^2\}.
```

An anchor is the translation $`u`$. Its pattern is globally unique if no other anchor anywhere in $`\mathbb Z^2`$ has that pattern. Define the difference erosion

```math
E_i(S)=\{s\in\mathbb Z^2:s+\mathrm{supp}Q_i\subseteq S\}.
```

Let $`\mathcal L_i(S)`$ consist of the integer heights at which $`E_i(S)`$ contains $`\kappa_i`$ consecutive lattice points in direction $`v_i`$. In the other proof files this row set is written $`L_i(S)`$; the argument $`S`$ distinguishes it from the tail field $`L_i`$.

For any specified nonzero extreme rows $`\mu_i^-\le\mu_i^+`$, put

```math
C_i(S)=(\mu_i^- -\mathcal L_i(S))\cap(\mu_i^+-\mathcal L_i(S)),
\qquad
\mathcal A(S)=\{u\in\mathbb Z^2:
\pi_i(u)\in C_i(S)\text{ for at least two distinct }i\}.
```

In Case B use the global support extrema of $`f_i`$. For arbitrary stars use the common extrema of the two end fields, and write $`C_i^{\mathrm{end}}`$, $`\mathcal A^{\mathrm{end}}`$. The end fields and their common extrema are proved to exist in [the all-star extension](all-star-extension.md). Empty row and anchor sets are allowed.

The [Case B anchor lemma](signals-and-anchors.md) states that every anchor in $`\mathcal A(S)`$ is globally unique, hence $`P_\theta(S)\ge|\mathcal A(S)|`$. Its proof also supplies the extremal-row recovery used for general stars.

## 2. Main results

Fix any compact convex $`P\subset\mathbb R^2`$ with nonempty interior, independently of the positive integer $`n`$. Define

```math
S_n=nP\cap\mathbb Z^2,\qquad
p_i^-=\min_{x\in P}\pi_i(x),\qquad p_i^+=\max_{x\in P}\pi_i(x),
```

```math
B_i(P)=\{x:p_i^-\le\pi_i(x)\le p_i^+\},\qquad
\mathcal M(P)=\bigcup_{i\lt j}(B_i(P)\cap B_j(P)),\qquad
\Lambda(P)=\mathrm{area}(\mathcal M(P)).
```

Each pair intersection is bounded, because the directions are nonparallel. Area is Euclidean, with the fundamental lattice square of area one.

**Theorem 1 (exact leading term for every CN star).** For the chosen decomposition and every such fixed $`P`$,

```math
P_\theta(S_n)=\Lambda(P)n^2+O_{\theta,P}(n).
```

The coefficient depends only on the direction set and the window shape. The fixed decomposition, periods, tails and encoding may enter the error constant. Neither rationality nor a polygonal or smooth boundary is required. The proof combines [the geometric estimates and matching upper bound](leading-term.md) with [the all-star lower bound](all-star-extension.md).

**Theorem 2 (explicit all-star finite-window bound).** For every nonempty finite $`S`$, including nonconvex windows, the fixed certificate data and marker constructed in [the finite-window proof](finite-window-bound.md) give finite sets

```math
\mathcal J(S)=\bigcap_{f\in F_0}(f-S),\qquad
\mathcal O_R(S)=([-R,R]^2-\mathrm{conv}S)\cap\mathbb Z^2,
```

and

```math
P_\theta(S)\ge
\left|\mathcal J(S)\cup
\left(\mathcal A^{\mathrm{end}}(S)\setminus\mathcal O_R(S)\right)\right|.
```

Every counted anchor has a globally unique pattern. All constants are computable from a specified star decomposition with finite tail and transition data; their automatic construction is included in that proof. For a lattice-convex window one may take the maximum with CN Theorem T (p. 4; Theorem 7.3, p. 18), $`|S|+1`$. The displayed new bound may be weaker on small windows.

**Theorem 3 (one three-line second-order formula).** For the binary star

```math
c(x,y)=\sum_{i=0}^2\mathbf1[y-2ix=i],\qquad
S_n=\{(x,y)\in\mathbb Z_{\ge0}^2:x+y\le n\},
```

```math
P_c(S_n)=\frac{31}{4}n^2+\frac{11}{4}n+O(1).
```

The lines have disjoint lattice supports, so the sum is binary. The [proof](three-line-second-order.md) establishes the linear term by a collision classification, without a numerical fit. It asserts this example's linear coefficient, not a general second-order formula or an eventual quasipolynomial.

## 3. Observation centers and Newton polygons

The Laurent polynomial ring over characteristic zero is an integral domain, and Newton polygons add under nonzero multiplication (the standard fact used in CN §2.2, equation (2.3), p. 8). Thus

```math
\mathrm{conv}(\mathrm{supp}Q_i)=\sum_{j\ne i}[0,H_j],\qquad
0\in\mathrm{conv}(\mathrm{supp}Q_i).
```

Taking convex hulls of $`s+\mathrm{supp}Q_i\subseteq S`$ gives $`s\in\mathrm{conv}S`$. Consequently

```math
E_i(S)\subseteq\mathrm{conv}S.
```

This is the observation-center fact needed in both the all-star asymptotic proof and Theorem 2. It does not assert that $`E_i(S)\subseteq S`$ for nonconvex windows.

For example

```math
Q=(T^{(1,0)}-1)(T^{(0,1)}-1)(T^{(-1,-1)}-1)
```

has zero coefficient at the origin: the empty subset and the full subset cancel there. Taking $`S=\mathrm{supp}Q`$ gives $`0\in E(S)`$ but $`0\notin S`$. The convex hull still contains zero. Formal exponent sets and actual polynomial supports must therefore be distinguished.

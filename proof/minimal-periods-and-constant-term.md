# Exact constant terms and minimal periods, with an obstruction to additive boundary formulas

[Overview](../README.md) · [Coefficient supplement](explicit-coefficients.md) · [Scope](remarks-and-scope.md)

**Status: established in one fresh-context AI adversarial review.** No fatal or proof-repair issue was found; the wording correction on primitive tangents in §3.1 is incorporated. All computations below are finite cross-checks, not premises of the proofs. The general construction uses the previously established finite representative formula. It does not claim a new quantifier-elimination theorem, a new vertex-cone identity, or a literature priority result.

## 1. Results and precise scope

Let $`\theta:\mathbb Z^2\to\mathbb F_p`$ be a star with a specified finite description: at least two primitive pairwise nonparallel directions $`v_i`$, directional periods, two doubly periodic tails for each component, and finitely many transition-row tables. Each component is not doubly periodic. A doubly periodic background is permitted. Let $`P`$ be a compact rational convex polygon with nonempty interior, and write

```math
f(n)=P_\theta(nP\cap\mathbb Z^2),\qquad
\pi_i(x)=\det(v_i,x),
```

```math
B_i(P)=\{x:\min_P\pi_i\le\pi_i(x)\le\max_P\pi_i\},\qquad
\Lambda=\mathrm{area}\left(\bigcup_{i\lt j}(B_i(P)\cap B_j(P))\right).
```

The previously proved input is the area theorem $`f(n)=\Lambda n^2+O(n)`$ and a finite quantifier-free integer-linear formula selecting exactly one anchor per pattern. Its explicit construction is in [the coefficient supplement](explicit-coefficients.md), Sections 2–3, which received its own fresh-context review (see the [review-status record](../sources/review-status.md)). This page uses those mathematical results, not numerical observations or review labels as substitutes for proofs.

The following are established here.

**Theorem A — vertex-velocity period and explicit coefficients.** From the representative formula one can construct a positive integer $`q`$, finitely many stable polygon, segment and point families, and the least positive integer $`E`$ making every eventual vertex velocity integral. Then $`M=qE`$ is a sufficient period, and

```math
f(n)=\Lambda n^2+\beta_{n\bmod M}n+\gamma_{n\bmod M}
\qquad(n\text{ sufficiently large}).
```

Sections 3–5 give finite rational sums for both coefficients, including an explicit formula for $`\gamma`$. They involve the lattice points of fixed half-open fundamental parallelograms of vertex cones. The number of summands is independent of $`n`$. There is no interpolation, evaluation of unknown complexity values, or unproved stabilization test.

Applied to the already normalized cells in the old supplement, the same theorem replaces the old full period $`M_\beta\Delta J`$ by $`M_\beta E_{\mathrm{old}}`$, with $`E_{\mathrm{old}}\mid\Delta`$. Thus the factor $`J`$ is unnecessary for sufficiency. A second, often smaller, construction in Section 2 starts directly with the coordinate-congruence modulus $`q`$, without first imposing $`G`$.

**Theorem B — exact minimal periods.** The rational coefficient tables furnished by Theorem A determine, by exact cyclotomic divisibility, the minimal period of $`\beta`$, the minimal period of $`\gamma`$, and the minimal period of the entire eventual quasipolynomial. The three answers are specified in Section 6. This is a terminating algorithm on the full star input; no period guessed from a finite run is used.

**Theorem C — a precise obstruction to two natural simpler formulas.** On the fixed integral triangle and the fixed three affine lattice lines

```math
P=\mathrm{conv}\{(0,0),(1,0),(0,1)\},\qquad
\theta_c(x,y)=\sum_{i=0}^2 c_i\mathbf1[y-2ix=i],
\qquad c_i\in\mathbb F_p^\times,\quad p\ge3,
```

the direct coefficient formula is

```math
\boxed{\ \beta(\theta_c,P)=3-\frac14\mathbf1[c_0=c_2]\ }.
```

All directions, affine supports, widths, minimal directional periods and tails are unchanged when the nonzero colors vary. Consequently:

1. these uncolored geometric and period data do not determine $`\beta`$;
2. even allowing arbitrary dependence on each individual colored component, a sum of independent component contributions cannot determine $`\beta`$.

The second assertion follows from a nonzero mixed difference, not just a comparison of two examples. It also rules out a sum of independent per-component boundary corrections. It does **not** rule out a compact formula that explicitly includes interactions between colored profiles.

This gives the requested effective characterization of the all-star minimal period and a compiled explicit formula for the constant term. The unrestricted direct formula in the original star and window boundary data remains open. The negative result has exactly the natural but limited meaning just stated.

## 2. Cells without primitive-normal rounding

### 2.1 The established representative input

For clarity, the finite-description input has the following content. A color predicate for $`\theta(z)`$ is a finite Boolean combination of integer linear inequalities and divisibilities with fixed moduli: use the two tail tables off the transition rows, the finite tangential tables on those rows, and the finite color-sum table. Equality of two window patterns is expressed by quantifying over the integer sites satisfying the facet inequalities of $`nP`$.

Order anchors by $`(|u_1|+|u_2|,u_1,u_2)`$ lexicographically. Every nonempty pattern class has a least anchor, since the minimum norm exists and its shell is finite. The representative predicate says that no smaller anchor has the same pattern. The finite elimination proved in the input supplement gives an equivalent quantifier-free formula $`F(n,u)`$, with inequality atoms

```math
a_h\cdot u\le b_hn+c_h
```

and divisibility atoms of fixed positive moduli. All coefficients are integers, and

```math
\#\{u\in\mathbb Z^2:F(n,u)\}=f(n)<\infty
\qquad(n\ge1).
```

No coordinate bound on the original pattern classes is assumed.

### 2.2 Fixing congruences

Let $`q`$ be the least common multiple of the divisibility moduli, with empty least common multiple equal to one. Fix $`0\le r\lt q`$ and $`\rho\in\{0,\ldots,q-1\}^2`$, and set

```math
n=qk+r,\qquad u=qz+\rho.
```

All divisibility atoms now have constant truth values. Enumerate the truth assignments to the inequality atoms accepted by $`F`$. Replace a false inequality by its exact integer complement; thus all chosen conditions again have the form $`a\cdot u\le bn+c`$. For integer $`z`$, this is equivalent to

```math
\boxed{\ a\cdot z\le bk+
   \left\lfloor\frac{br+c-a\cdot\rho}{q}\right\rfloor\ }.
\qquad\text{(1)}
```

The reason is that $`a\cdot z`$ and $`bk`$ are integers. **Do not divide $`a`$ by its coordinate gcd here.** No primitive support-line assumption is needed in the cone argument.

An atom with $`a=0`$ has an effectively determined eventual truth value. Each remaining accepted truth assignment defines a closed rational polyhedral family

```math
C(k)=\{z\in\mathbb R^2:a_j\cdot z\le b_jk+d_j\}.
\qquad\text{(2)}
```

These cells need be disjoint only on their integer points, which follows from the truth assignments and coordinate residues.

If the recession cone of a cell is nonzero, it contains a nonzero integer vector. Any integer point in that cell would generate infinitely many by recession translation, contradicting the established finite representative fiber. Such cells have no integer points at the relevant positive scales and can be discarded. This does not claim their real fibers are empty.

All other cells are bounded or empty. Candidate vertices solve two independent boundary equations, so have the form

```math
p_v(k)=k\alpha_v+\delta_v,\qquad
\alpha_v,\delta_v\in\mathbb Q^2.
\qquad\text{(3)}
```

Their feasibility is determined by finitely many affine signs in $`k`$. Eventual feasibility, coincident vertices, cyclic order and dimension can therefore be determined exactly. One may additionally use root bounds on the finitely many determinant polynomials of degree at most two to fix the order and dimension. Taking an integer above all these bounds gives a certified common stable threshold. This is an algebraic threshold, not one inferred by sampling counts.

Retain the actual eventual vertices of every nonempty cell, including endpoints of segments and points. Set

```math
E=\mathop{\mathrm{lcm}}\{\mathrm{den}(\alpha_{v,1}),
                       \mathrm{den}(\alpha_{v,2})
                       :v\text{ retained}\},\qquad M=qE.
\qquad\text{(4)}
```

Here $`\mathrm{den}`$ is the positive reduced denominator; an empty list gives $`E=1`$. All lists include all $`r,\rho`$ and accepted cells. Cramer's rule gives the optional coarse bound

```math
E\mid
\mathop{\mathrm{lcm}}\bigl(\{1\}\cup
 \{|\det(a,a')|:\det(a,a')\ne0\}\bigr),
\qquad\text{(5)}
```

using the integer, not necessarily primitive, cell normals.

For $`0\le s\lt M`$, write uniquely $`s=qh+r`$, with $`0\le h\lt E`$ and $`0\le r\lt q`$. Substitution $`k=Et+h`$ changes each eventual vertex to

```math
p_v(Et+h)=t a_v+b_v,\qquad
a_v=E\alpha_v\in\mathbb Z^2,\qquad
b_v=h\alpha_v+\delta_v\in\mathbb Q^2.
\qquad\text{(6)}
```

The constants $`b_v`$ need not be vertices of a valid cell at parameter zero. They are affine extrapolations of the stable vertex branches. The formulas below are designed to allow precisely this situation.

## 3. A self-contained vertex-cone coefficient formula

### 3.1 The exact finite cone identity

Let $`K`$ be a bounded rational two-dimensional polygon, with redundant collinear vertices removed. At vertex $`b_v`$, let $`w_{v,1},w_{v,2}`$ be the two primitive integer vectors along its incident edges, both pointing into the edges from the vertex. Define the half-open parallelogram and its lattice set

```math
\Pi_v=b_v+[0,1)w_{v,1}+[0,1)w_{v,2},\qquad
Z_v=\Pi_v\cap\mathbb Z^2.
\qquad\text{(7)}
```

The two rays are independent. Every lattice point in the tangent cone at $`b_v`$ is uniquely

```math
z+m_1w_{v,1}+m_2w_{v,2},
\qquad z\in Z_v,\quad m_1,m_2\in\mathbb Z_{\ge0}.
```

Indeed, take the integer and fractional parts of its two real ray coordinates. In particular $`|Z_v|=|\det(w_{v,1},w_{v,2})|`$: the parallelogram is a translated fundamental domain for the sublattice generated by the two rays.

The rational generating-function identity is

```math
\sum_{x\in K\cap\mathbb Z^2}X^x
 =
\sum_v
 \frac{\sum_{z\in Z_v}X^z}
 {(1-X^{w_{v,1}})(1-X^{w_{v,2}})}.
\qquad\text{(8)}
```

Here $`X^{(a,b)}=X_1^aX_2^b`$; the right side is interpreted as a rational function, not as simultaneously convergent infinite sums.

For completeness, (8) has the following elementary proof in dimension two. Let $`T_v`$ be the closed tangent cone at a vertex and $`H_e`$ the closed supporting half-plane at an edge containing $`K`$. The pointwise identity is

```math
\mathbf1_K=\sum_v\mathbf1_{T_v}-\sum_e\mathbf1_{H_e}
                         +\mathbf1_{\mathbb R^2}.
\qquad\text{(9)}
```

For a point in $`K`$, every term is present and the alternating count is one. For a point outside $`K`$, the violated facets form a nonempty proper consecutive visible arc between the two supporting tangents from that point. Thus at least one facet is satisfied. If $`j`$ facets are satisfied, their satisfied adjacent pairs number $`j-1`$. The right side is therefore $`(j-1)-j+1=0`$. The same argument includes points on facet extensions, with violation taken strictly.

To turn (9) into a rigorous finite generating-function calculation, choose one **primitive** integer tangent $`w_e`$ along each edge and apply the product of lattice difference operators corresponding to

```math
D(X)=\prod_e(1-X^{w_e}).
```

The factor $`1-X^w`$ means the pointwise operator $`g(z)\mapsto g(z)-g(z-w)`$. Each half-plane term is killed by its edge-tangent factor, and the whole-plane term is killed as well. For a vertex cone, its two incident factors, after changing ray signs if necessary, reduce its lattice indicator to a signed translate of the finite parallelogram indicator in (7). All remaining differences preserve finite support. Thus every generating-function sum after applying $`D`$ is finite. The parallelogram decomposition proves that the resulting equality is $`D`$ times (8). Divide by the nonzero Laurent polynomial $`D`$. This proves (8) without assigning a divergent series a value.

### 3.2 Extracting the constant coefficient

Choose an integer vector $`c`$ not orthogonal to any edge ray. It can be found by trying $`(1,L)`$ for positive integers $`L`$: only finitely many are excluded. Put

```math
d_1=c\cdot w_{v,1},\qquad d_2=c\cdot w_{v,2},
\qquad d_1d_2\ne0,
```

and define the rational quadratic polynomial

```math
\boxed{\ H_{d_1,d_2}(x)
 =\frac{x^2}{2d_1d_2}
 -\frac{x}{2}\left(\frac1{d_1}+\frac1{d_2}\right)
 +\frac14+\frac1{12}\left(\frac{d_1}{d_2}+\frac{d_2}{d_1}\right)\ }.
\qquad\text{(10)}
```

Substitute $`X^z=\exp(\varepsilon\,c\cdot z)`$ in (8) and take the constant Laurent coefficient at $`\varepsilon=0`$. The expansion

```math
\frac1{1-e^{d\varepsilon}}
 =-\frac1{d\varepsilon}+\frac12-\frac{d\varepsilon}{12}
  +O(\varepsilon^3)
```

shows that

```math
[\varepsilon^0]\,
\frac{e^{x\varepsilon}}
 {(1-e^{d_1\varepsilon})(1-e^{d_2\varepsilon})}
 =H_{d_1,d_2}(x).
```

The left side of (8) is finite, and its constant coefficient is the number of lattice points. Therefore

```math
|K\cap\mathbb Z^2|
 =\sum_v\sum_{z\in Z_v}H_{d_{v,1},d_{v,2}}(c\cdot z).
\qquad\text{(11)}
```

Negative values of $`d_1`$ or $`d_2`$ are allowed and essential; no absolute values are inserted in (10). The resulting total is independent of the admissible choice of $`c`$, because (11) equals the exact lattice count.

### 3.3 Moving vertices, including thin polygons

Now let $`K(t)`$ be any stable two-dimensional cell of (6). Its facet directions and tangent rays are fixed. At vertex $`t a_v+b_v`$, the tangent cone is the integer translate by $`t a_v`$ of the cone at the extrapolated point $`b_v`$. Its parallelogram lattice set is exactly $`t a_v+Z_v`$, where (7) is defined at $`b_v`$.

Consequently, for all sufficiently large integers $`t`$,

```math
|K(t)\cap\mathbb Z^2|
 =
\sum_v\sum_{z\in Z_v}
 H_{d_{v,1},d_{v,2}}\bigl(t\eta_v+c\cdot z\bigr),
\qquad \eta_v=c\cdot a_v.
\qquad\text{(12)}
```

This is an exact quadratic polynomial in $`t`$. It applies to polygons whose area is only linear, polygons of bounded diameter, and polygons with no lattice points. Cancellation in the finite sum handles these cases automatically.

It is not necessary for the extrapolated collection $`(b_v)`$ to form a polygon. Identity (8) is applied only at large $`t`$, where the stable polygon exists; each of its individual cones is then translated back by its own integer vector. In particular, one must not replace the constant term of (12) by a lattice count of a potentially nonexistent polygon at $`t=0`$.

For example, the cell $`0\le x\le t-7,\ 0\le y\le t-11`$ has eventual count $`(t-6)(t-10)=t^2-16t+60`$, although its fiber at zero is empty. The constant coefficient is sixty, not that empty fiber's count.

For a convenient fully expanded version put

```math
K_v=|Z_v|,\qquad T_v=\sum_{z\in Z_v}c\cdot z,\qquad
U_v=\sum_{z\in Z_v}(c\cdot z)^2.
```

The contributions of this vertex to the coefficients of $`t^2,t,1`$ are respectively

```math
A_v=\frac{K_v\eta_v^2}{2d_1d_2},
\qquad\text{(13)}
```

```math
B_v=\frac{\eta_v T_v}{d_1d_2}
 -\frac{K_v\eta_v}{2}\left(\frac1{d_1}+\frac1{d_2}\right),
\qquad\text{(14)}
```

```math
\boxed{\ C_v=
 \frac{U_v}{2d_1d_2}
 -\frac{T_v}{2}\left(\frac1{d_1}+\frac1{d_2}\right)
 +\frac{K_v}{4}
 +\frac{K_v}{12}\left(\frac{d_1}{d_2}+\frac{d_2}{d_1}\right)\ }.
\qquad\text{(15)}
```

These moments can be enumerated exactly in a fixed rational bounding box for (7). The amount of enumeration may be large, but does not grow with the dilation parameter.

## 4. Lower-dimensional cells and the complete coefficient sum

### 4.1 Segments

Let a stable segment have endpoints

```math
p_-(t)=t a_-+b_-,\qquad
p_+(t)=t a_++b_+,
```

with $`a_\pm\in\mathbb Z^2`$, and let its eventual orientation have primitive tangent $`w`$. Take the primitive normal $`\nu=(-w_2,w_1)`$. Its containing line is

```math
\nu\cdot x=t\,\nu\cdot a_-+\nu\cdot b_-.
```

If $`\nu\cdot b_-\notin\mathbb Z`$, the segment has no integer points for any integer $`t`$. If this constant is integral, find $`x_0\in\mathbb Z^2`$ on the line $`\nu\cdot x=\nu\cdot b_-`$ by Bezout's identity. Write

```math
b_\pm=x_0+\tau_\pm w,\qquad
a_+-a_-=Lw.
```

Primitivity of $`w`$ and integrality of $`a_\pm`$ imply $`L\in\mathbb Z`$. Eventual orientation implies $`L\ge0`$. Translation by $`-t a_-`$ reduces the integer parameters on the segment to

```math
\tau_-\le j\le Lt+\tau_+,\qquad j\in\mathbb Z.
```

Thus its coefficient triple is

```math
(A_C,B_C,C_C)
 =\left(0,\ L,\ \lfloor\tau_+\rfloor-\lceil\tau_-\rceil+1\right).
\qquad\text{(16)}
```

This equality is eventual if the extrapolated endpoints are reversed at zero. When $`L=0`$ the actual segment length is constant and nonnegative; the displayed count is already nonnegative, possibly zero. This explicitly treats nonintegral support lines, which were avoided by primitive-normal rounding in the old construction.

### 4.2 Points and empty cells

For a point $`p(t)=ta+b`$ with $`a\in\mathbb Z^2`$, the coefficient triple is

```math
(A_C,B_C,C_C)=\left(0,0,\mathbf1[b\in\mathbb Z^2]\right).
\qquad\text{(17)}
```

Empty cells contribute zero. The empty unbounded cells discarded in Section 2 contribute zero as well.

### 4.3 Returning to the original variable

For each residue $`s=qh+r`$ modulo $`M=qE`$, add (13)–(15), (16), and (17) over the retained cells and coordinate residues belonging to $`r`$. Denote the resulting total by $`(A_s,B_s,C_s)`$. Disjointness on integer points gives, exactly for large $`t`$,

```math
f(Mt+s)=A_st^2+B_st+C_s.
\qquad\text{(18)}
```

Comparing its quadratic coefficient with the proved area theorem yields $`A_s=M^2\Lambda`$. Substituting $`t=(n-s)/M`$ gives the explicit formulas

```math
\boxed{\ \beta_s=\frac{B_s}{M}-2s\Lambda,\qquad
\gamma_s=C_s-\frac{sB_s}{M}+s^2\Lambda\ }.
\qquad\text{(19)}
```

Equivalently $`\gamma_s=C_s-\Lambda s^2-\beta_s s`$. Both formulas use the original scale $`n`$, rather than silently identifying a coefficient in $`t`$ with one in $`n`$.

Everything in (19) is a finite rational expression in the constructed cells and their cone parallelograms. This proves Theorem A and supplies the promised general explicit constant term. The formula remains compiled from representatives; it is not advertised as a direct formula in the original star boundary profiles.

## 5. Why vertex velocities suffice, and a worked constant term

An intercept denominator does not itself force a dilation period. Once a residue is fixed and every vertex moves by an integer vector, each fixed tangent cone has a lattice numerator translated by that vector; (12) is polynomial. This is why no extra slicing modulus occurs.

In the old primitive-normal construction, a vertex velocity solves two equations with primitive integer normals and integer right-side slopes. Its coordinate denominators divide their determinant. Consequently $`E_{\mathrm{old}}\mid\Delta`$ and

```math
M_{\mathrm{old,cone}}=M_\beta E_{\mathrm{old}}
\quad\text{is sufficient, with}\quad
M_{\mathrm{old,cone}}\mid M_\beta\Delta.
\qquad\text{(20)}
```

This statement concerns the same old cells. The alternative unrounded construction (1) may use a different subdivision and a different sufficient modulus; no unproved comparison of all such moduli is needed. Minimal-period reduction will remove any surplus in either construction.

For example, consider the elementary cell

```math
K(k)=\{(x,y):x\ge0,\ y\ge0,\ 2x+3y\le k\}.
```

Its vertex velocities are $`(0,0),(1/2,0),(0,1/3)`$, so $`E=6`$. The three tangent-ray pairs are

```math
((1,0),(0,1)),\qquad
((-1,0),(-3,2)),\qquad
((0,-1),(3,-2)).
```

Take $`c=(1,1)`$, which is nonorthogonal to all six rays. For residue $`h`$, use bases $`(0,0),(h/2,0),(0,h/3)`$ in (7). The three parallelograms contain respectively $`1,2,3`$ lattice points. Applying (13)–(19) gives

```math
|K(k)\cap\mathbb Z^2|
=\frac{k^2}{12}+\frac{k}{2}+\gamma_{k\bmod6},
```

```math
(\gamma_0,\ldots,\gamma_5)
=\left(1,\frac5{12},\frac23,\frac34,\frac23,\frac5{12}\right).
\qquad\text{(21)}
```

These entries are evaluations of (15) and (19), not a fit to six counts. The coefficient of $`k`$ has minimal period one. The constant term has minimal period six: it has neither period one, two nor three. The old $`\Delta J`$ prescription gives eighteen for this cell.

The vertex-velocity period can still be nonminimal. A moving box with a nonintegral constant-width transverse strip may have no lattice points at all, in which case all cone contributions cancel and the count has period one regardless of vertex denominators. In a disjoint sum of cells, cancellation can occur between cells as well. This is the reason to reduce the final coefficient tables, rather than declare a geometric denominator minimal.

## 6. Exact minimal periods and a terminating algorithm

Let $`M`$ be any proved sufficient period, and let the tables in (19) be indexed by $`0\le s\lt M`$. Form

```math
\mathcal B(X)=\sum_{s=0}^{M-1}\beta_sX^s,\qquad
\mathcal G(X)=\sum_{s=0}^{M-1}\gamma_sX^s.
\qquad\text{(22)}
```

Let $`\Phi_d(X)`$ denote the $`d`$-th cyclotomic polynomial. Define

```math
\mathcal D_\beta=\{d:d\mid M,\ \Phi_d\nmid\mathcal B\},\qquad
\mathcal D_\gamma=\{d:d\mid M,\ \Phi_d\nmid\mathcal G\}.
```

Divisibility here is exact in $`\mathbb Q[X]`$; every polynomial divides the zero polynomial. Empty least common multiples are one.

**Theorem B, explicit form.** The three minimal periods are

```math
\boxed{\ 
p_\beta=\mathop{\mathrm{lcm}}\mathcal D_\beta,\qquad
p_\gamma=\mathop{\mathrm{lcm}}\mathcal D_\gamma,\qquad
p_f=\mathrm{lcm}(p_\beta,p_\gamma).
\ }
\qquad\text{(23)}
```

By $`p_f`$ we mean the least positive period for a polynomial formula on each residue class valid for all sufficiently large $`n`$. A finite irregular initial segment is irrelevant.

**Proof of the coefficient-period assertions.** A periodic table $`b_s`$ has generating function

```math
\sum_{n\ge0}b_{n\bmod M}X^n
 =\frac{\sum_{s=0}^{M-1}b_sX^s}{1-X^M}.
\qquad\text{(24)}
```

The denominator is squarefree. Cancellation at the primitive $`d`$-th roots is equivalent, over the rationals, to divisibility of the numerator by $`\Phi_d`$. Therefore its reduced denominator has precisely the factors indexed by the uncanceled $`d`$.

A period $`h`$ implies that this reduced denominator divides $`1-X^h`$, so every such $`d`$ divides $`h`$. Conversely, if $`h`$ is the least common multiple of these orders, the reduced denominator divides $`1-X^h`$. The rational function in (24) is proper; multiplying by $`1-X^h`$ gives a polynomial of degree less than $`h`$. Expansion then gives an $`h`$-periodic sequence from index zero. This proves minimality, including the identically zero sequence.

**Proof of the full-period assertion.** Period $`p_\beta`$ and period $`p_\gamma`$ immediately give period $`\mathrm{lcm}(p_\beta,p_\gamma)`$ for the eventual quasipolynomial. Conversely, suppose $`h`$ is any eventual quasipolynomial period of $`f`$. Its residue polynomials have degree at most two because $`f(n)=O(n^2)`$, and their quadratic coefficients must be $`\Lambda`$ by the area theorem. Refine the $`h`$- and $`M`$-residue classes to their common refinement. On each infinite progression, equality with (19) is equality of polynomials, so its linear and constant coefficients must equal the corresponding $`\beta_s,\gamma_s`$. Hence $`h`$ is a period of both coefficient sequences. Their minimal periods divide $`h`$, proving (23). This also proves that the answer does not depend on the chosen sufficient modulus, representative order, or cell subdivision. ∎

For direct implementation, one can equivalently test all divisors $`h\mid M`$ and retain the least $`h`$ satisfying

```math
\beta_s=\beta_{s+h},\qquad \gamma_s=\gamma_{s+h}
\quad\text{for every }s\bmod M.
\qquad\text{(25)}
```

The cyclotomic formulation identifies exactly which frequencies survive cancellation; the cyclic-word formulation needs only rational comparisons.

There is also a single rational-function characterization. Up to a polynomial accounting for the finite initial segment,

```math
\sum_{n\ge1}f(n)X^n
\equiv
\Lambda\,\frac{X(1+X)}{(1-X)^3}
+X\frac{d}{dX}\left(\frac{\mathcal B(X)}{1-X^M}\right)
+\frac{\mathcal G(X)}{1-X^M}
\pmod{\mathbb Q[X]}.
\qquad\text{(26)}
```

At a nontrivial root where $`\mathcal B`$ does not vanish, the derivative term has a genuine double pole, which the last term cannot cancel. Where $`\mathcal B`$ vanishes, that entire derivative term is analytic, and a pole remains exactly when $`\mathcal G`$ does not vanish. Thus $`p_f`$ is also the least common multiple of the orders of the roots of unity that remain poles. The root one contributes order one to this least common multiple.

A complete decision algorithm is now explicit:

1. Construct the finite color predicates and representative formula, and apply the previously proved terminating elimination.
2. Form (1)–(3), discard integer-empty unbounded cells, and compute exact eventual vertex branches and a stable threshold.
3. Compute $`E,M`$, enumerate (7), and evaluate (13)–(19), including (16)–(17).
4. Use either rational polynomial division in (23) or exact comparisons in (25).
5. Output $`p_\beta,p_\gamma,p_f`$ and the corresponding reduced coefficient tables.

Every loop is finite and every test is rational or integral arithmetic. No claim of practical polynomial running time, a small elimination output, or a complete implemented compiler is made.

## 7. A direct colored example and the obstruction proof

This section proves Theorem C from complete pattern representations. It does not use elimination or Theorem A.

### 7.1 Raw representations

For $`i=0,1,2`$ let

```math
Q_i(t)=\{(x,y)\in\mathbb Z_{\ge0}^2:x+y\le n,\ y-2ix=t\},
\qquad \mathcal Q_i=\{t:Q_i(t)\ne\varnothing\}.
```

An anchor $`(X,Y)`$ produces section parameters

```math
t_i=i-Y+2iX=t_0+i d,\qquad d=1+2X\text{ odd}.
\qquad\text{(27)}
```

The three fixed real lines meet at $`(-1/2,0)`$ and have disjoint lattice supports. Thus the support of every pattern is the disjoint union of its active sections, with colors $`c_i\ne0`$; no cancellation occurs.

For $`n\ge3`$, elementary enumeration gives

```math
\mathcal Q_0=[0,n]\cap\mathbb Z,
```

```math
\mathcal Q_1=([-2n,n]\cap\mathbb Z)\setminus\{-2n+1\},
```

```math
\mathcal Q_2=([-4n,n]\cap\mathbb Z)
 \setminus\{-4n+1,-4n+2,-4n+3,-4n+6,-4n+7,-4n+11\}.
\qquad\text{(28)}
```

To verify the gaps, write $`x=n-r`$. The possible heights for slope $`2i`$ at that column form the integer interval $`[-2in+2ir,-2in+(2i+1)r]`$. Their union gives exactly (28). In particular

```math
|\mathcal Q_0|=n+1,\qquad|\mathcal Q_1|=3n,\qquad
|\mathcal Q_2|=5n-5.
\qquad\text{(29)}
```

Let $`N_{ij}`$ count anchors having both sections $`i,j`$. Their height congruences from (27) are: $`t_1-t_0`$ odd, $`t_2-t_0\equiv2\pmod4`$, and $`t_2-t_1`$ odd, respectively. For two intervals with boundedly many holes, residue counts are their sizes divided by the fixed modulus plus bounded deviations summing to zero. Multiplying and summing over matching residues cancels both linear cross terms. Therefore

```math
N_{01}=\frac32n^2+\frac32n+O(1),\quad
N_{02}=\frac54n^2+O(1),\quad
N_{12}=\frac{15}{2}n^2-\frac{15}{2}n+O(1).
\qquad\text{(30)}
```

Every compatible pair $`t_0,t_2`$ forces $`t_1=(t_0+t_2)/2`$. It is in the full interval $`[-2n,n]`$ and fails (28) only when $`t_0+t_2=-4n+2`$, allowing at most three pairs. Hence $`N_{012}=N_{02}+O(1)`$. Inclusion-exclusion gives the number of anchors with at least two sections:

```math
A_n=N_{01}+N_{02}+N_{12}-2N_{012}
    =\frac{31}{4}n^2-6n+O(1).
```

A complete finite representation list contains these anchors, one pure section for every $`i,t\in\mathcal Q_i`$, and zero. Every pure section is realized by moving its anchor far along its own direction, and zero is realized outside all three visibility strips. Its length is

```math
R_n=A_n+\sum_i|\mathcal Q_i|+1
    =\frac{31}{4}n^2+3n+O(1).
\qquad\text{(31)}
```

This raw count is independent of the nonzero colors. Only repeated representations remain to be analyzed.

### 7.2 Completeness of the collision classification

Call a section long if it has at least three lattice sites. A long relative line is intrinsic to its pattern: a competing representation with a different parallel line could cover at most two of these sites by its other two directions. Applying this argument to both representations recovers all long lines. Two long lines recover the anchor uniquely. A pure long pattern cannot equal a multi-line pattern, whose other nonempty sections are disjoint from its long line.

Each section with one or two sites lies in a bounded neighborhood of $`U=(n,0)`$ or $`V=(0,n)`$, with a bound independent of $`n`$. Explicitly, for slope $`k=2`$ or $`4`$ and nonnegative height $`t`$, such a section starts at $`(0,t)`$ and $`n-t<2(k+1)`$. For negative height its first point has $`0\le y_0\lt k`$, and absence of a third point implies $`n-x_0\lt y_0+2(k+1)`$. Horizontal short sections lie in the top two rows. Radius twenty in the sup norm is more than enough for all these sites.

Thus patterns with no long direction have only boundedly many representations: each short height belongs to a bounded list and two heights determine an anchor. Pure short patterns also have boundedly many entries.

Suppose exactly one direction $`i`$ is long and let $`E_{\mathrm{off}}`$ be the sites off its intrinsic line. In two distinct representations every off-line site must change its assigned direction: an unchanged assignment would give two fixed nonparallel lines and hence the same anchor. There are only two remaining directions. If there are at least three off-line sites, two belonged to the same direction and would both have to move to the other direction, impossible for two distinct sites on nonparallel lines. If there are two sites, they must exchange their assignments. They lie in the bounded vertex neighborhoods, so their choices and the possible anchors form a bounded list. These lists also include all their collision partners, by the same argument.

Consequently every unbounded collision family has exactly one long line and exactly one off-line site, to be assigned alternatively to the two remaining directions. Colors can only remove such support coincidences; they cannot introduce a new support coincidence.

### 7.3 The one linear family, and its color condition

Write the off-line site as $`(a,b)`$ and the long height as $`t`$. For long direction zero, assigning this site to direction one forces $`b-t`$ odd, whereas assignment to direction two forces $`b-t\equiv2\pmod4`$. These are incompatible. For long direction two the same two incompatible conditions occur on $`t-b`$.

For long direction one, the two assignments require $`t-b`$ odd. One alternative is a horizontal singleton, whose only location is $`V=(0,n)`$. The other is the direction-two singleton at the same site. Their parameter triples are

```math
(n,t,2t-n),\qquad (2t-n,t,n).
\qquad\text{(32)}
```

They are distinct integer anchors exactly when $`t-n`$ is odd. Both have no additional active section exactly when $`2t-n\notin\mathcal Q_0\cup\mathcal Q_2`$. The long-section requirement is $`|Q_1(t)|\ge3`$.

The two patterns in (32) agree **if and only if $`c_0=c_2`$**. Their long line already has the same color $`c_1`$; their sole off-line site has color $`c_0`$ or $`c_2`$. There is no third representation, since this site has only the two remaining direction assignments.

The number of these two-entry classes, before the color test, is

```math
F_n=\#\{t\in\mathcal Q_1:|Q_1(t)|\ge3,\ t-n\text{ odd},
                          \ 2t-n\notin\mathcal Q_0\cup\mathcal Q_2\}.
```

Apart from the six holes of $`\mathcal Q_2`$, the main range is $`[-2n,-3n/2)`$. The long-section condition removes only a bounded number of heights next to $`-2n`$, and parity retains half of this interval. Hence

```math
F_n=\frac n4+O(1).
\qquad\text{(33)}
```

All other representation excess is bounded by Section 7.2. Subtracting it and (33) from (31) proves

```math
f_c(n)=\frac{31}{4}n^2+
 \left(3-\frac14\mathbf1[c_0=c_2]\right)n+O(1).
\qquad\text{(34)}
```

This proves the direct coefficient formula in Theorem C for every nonzero $`c_1`$ and every prime $`p\ge3`$.

### 7.4 What simpler formulas are impossible

Keep $`c_1=1`$ and allow $`c_0,c_2\in\{1,2\}`$ in $`\mathbb F_3`$. The coefficient matrix is

| $`c_0\backslash c_2`$ | $`1`$ | $`2`$ |
| --- | ---: | ---: |
| $`1`$ | $`11/4`$ | $`3`$ |
| $`2`$ | $`3`$ | $`11/4`$ |

All supporting lines, lattice intersections, directional periods, zero tails, transition heights and the polygon are identical across this matrix. Therefore any formula using only those uncolored data gives the same answer in all four cases and is false.

More strongly, suppose on this fixed geometry that a formula had the component-separable form

```math
\beta(\theta_c,P)=G+\sum_{i=0}^2\Psi_i(F_i,P),
\qquad\text{(35)}
```

where each $`\Psi_i`$ may use the entire individual colored component, but no other component's color data, and $`G`$ is fixed by the uncolored geometry. The mixed difference in $`c_0,c_2`$ of (35) is zero. The actual mixed difference is

```math
\beta_{11}+\beta_{22}-\beta_{12}-\beta_{21}
=\frac{11}{4}+\frac{11}{4}-3-3=-\frac12.
\qquad\text{(36)}
```

This contradiction proves the second obstruction. Allowing multiple separate edge or vertex terms for each component does not help, since they can be absorbed into $`\Psi_i`$.

The obstruction is the equality relation between colors of different components at a reassignable boundary singleton. The formula (34) itself shows that a color-aware interaction term can be short. Accordingly (36) is not a proof that every general compact formula is impossible, nor a lower bound on the complexity of an algorithm that retains all color tables.

## 8. Verification, dependencies and remaining obstacles

The diagnostic `scripts/minimal_period_audit.py` computes coefficients from (13)–(17) before comparing them with independently counted integer vertical sections. It includes nonprimitive normals, rational vertex offsets, constant-width and empty thin cells, moving bounded polygons, nonintegral segment support lines, moving points, and extrapolated parameter-zero vertices that are not an actual polygon. Its deterministic output also compares cyclotomic reduction with cyclic-table equality.

For the colored three-line family it enumerates every multi-line anchor by solving all compatible pairs of section heights, then adds all globally realizable pure entries and zero. It does not call a finite translation box exhaustive. It checks all eight assignments of colors one and two on the specified range. These tests can discover mistakes but cannot prove an asymptotic, minimal period, or absence of further collisions; those conclusions follow from Sections 2–7.

The mathematical dependencies and limits are:

- The only all-star inputs are the established finite representative construction and area theorem. The new exact coefficients and period statements follow from the cone argument, including its proof here.
- The cone identity and cyclotomic cancellation are classical mechanisms. The contribution here is their explicit assembly for all representative cells, including degeneracies, an explicit constant-term expression, the removal of the slicing factor, and exact minimal-period reduction.
- This representative-cell construction does not itself give a formula directly in the original star's finite boundary profiles. A separate [integral-window formula](integral-windows.md) does so on integral convex polygons. The obstruction proves that uncolored geometry and independent-component terms are insufficient, not that all interaction-based formulas fail.
- This page alone does not decide all-star constant $`\beta`$ on integral polygons. The [integral-window theorem](integral-windows.md) proves it for every integral convex polygon; the example in Section 7 is consistent with that result.
- The algorithm is a decision procedure, not a claim of a practical implementation of the full star-to-formula compiler. Elimination size and parallelogram determinants can be large.
- One fresh-context AI adversarial review found no fatal or repairable issue. No external mathematician or proof assistant has certified this page.

For rational polygons beyond the integral case, a direct all-star boundary formula would require replacing the eliminated representative cells by a complete, directly computable classification of boundary-profile equivalences. The [integral-window proof](integral-windows.md), Appendix B, supplies such a classification in its stated integral scope. Such a classification must retain cross-component color and phase matching, and must prove that all omitted collisions have only bounded total excess, including their partners. The three-line singleton family is a concrete necessary test for any proposed classification.

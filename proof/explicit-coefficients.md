# Explicit second coefficients and sufficient periods from representative polygons

**Status: established in one fresh-context AI adversarial review.** Separately reviewed in a fresh GPT-6.1 Sol context on 25 October 2026; no fatal or proof-repair issue was found. This is an AI adversarial review, not human certification or proof-assistant verification. This remains a separate supplementary result; the convex area and finite-window bounds do not use it as a premise. The construction below uses elementary Presburger elimination and planar lattice geometry; no new general counting principle or priority claim is made.

The current status and the fingerprints of the reviewed proof and review are in the [review-status record](../sources/review-status.md).

## 1. What the formula supplies

Fix a star $`\theta:\mathbb Z^2\to\mathbb F_p`$ in the sense of the contest paper, with a finite specification of its directional periods, its two periodic tails for each component, and its finite transition-row tables. A doubly periodic background is allowed. Fix a compact rational convex polygon $`P`$ with nonempty interior, and set $`f(n)=P_\theta(nP\cap\mathbb Z^2)`$ for positive integers $`n`$. Put

```math
\pi_i(x)=\det(v_i,x),\qquad
B_i(P)=\{x:\min_P\pi_i\le\pi_i(x)\le\max_P\pi_i\},
\qquad
\Lambda=\mathrm{area}\left(\bigcup_{i\lt j}(B_i(P)\cap B_j(P))\right).
```

The already reviewed area theorem gives $`f(n)=\Lambda n^2+O(n)`$. The already reviewed definable-representative argument and Bogart–Goodrick–Woods give eventual quasipolynomiality. The present result makes a sufficient period and the linear coefficients effective by an explicit finite geometric construction.

**Theorem.** Construct the finite integer data and polygons in Sections 2–4. They give explicit positive integers $`M_\beta,\Delta,J`$, with $`M=M_\beta\Delta J`$, and rational numbers $`A_r,B_r`$, for $`0\le r\lt M_\beta`$, such that

```math
A_r=M_\beta^2\Lambda,\qquad
\boxed{\ \beta_r=\frac{B_r}{M_\beta}-2r\Lambda\ },
```

and, for all sufficiently large integers $`n`$,

```math
f(n)=\Lambda n^2+\beta_{n\bmod M_\beta}n+\gamma_{n\bmod M}
```

for rational constants $`\gamma_s`$. Thus $`M_\beta`$ is a sufficient period of the linear coefficient; $`M`$ is a sufficient period of the whole eventual quasipolynomial. Neither is claimed minimal.

Here $`B_r`$ is a finite sum of area-displacement coefficients, half the coefficients of lattice boundary length, and the lattice lengths of one-dimensional cells. All polygons, sums and least common multiples are specified by terminating operations on the input tables and facet inequalities. No unknown complexity values, empirical stabilization, or fitting enter the construction.

This is a **compiled geometric formula**: its polygons are those of a definable set of pattern representatives, after explicit quantifier elimination. It is not a short formula using only the original direction bands or the original window facets. Quantifier elimination may make it impractically large. This theorem alone does not prove that an integral input polygon always has constant $`\beta`$, nor that its vertex-denominator modulus suffices in general. Constancy on integral convex polygons is proved separately in [the integral-window theorem](integral-windows.md).

## 2. An explicit formula for one representative per pattern

Choose an integer basis $`(v_i,e_i)`$ with determinant one. Write $`t_i(z)=\det(v_i,z)`$, $`s_i(z)=\det(z,e_i)`$, so $`z=s_i(z)v_i+t_i(z)e_i`$. On the finitely many transition heights, the value of component $`i`$ is specified by $`t_i(z)`$ and $`s_i(z)\bmod k_i`$. Each tail is specified by $`z\bmod q_0`$ for some common positive integer $`q_0`$: a finite-index period lattice contains $`q_0\mathbb Z^2`$. The same holds for a periodic background.

Consequently each color predicate $`\Theta_a(z)`$, meaning $`\theta(z)=a`$, is a finite Boolean combination of integer linear inequalities and divisibilities by fixed positive integers. Concretely, use the left table when $`t_i<\ell_i`$, the right table when $`t_i>r_i`$, and the transition table otherwise; then take the finite disjunction over color tuples whose sum is $`a`$ in $`\mathbb F_p`$. Empty transition strips and cancellations cause no exception. Each point satisfies exactly one color predicate.

Write the window as $`A_js\le b_jn`$, with integer $`A_j,b_j`$, by clearing fixed positive denominators. Define

```math
E(n,u,v)\ \Longleftrightarrow\ 
\forall s\in\mathbb Z^2\ 
\left[\left(\bigwedge_j A_js\le b_jn\right)
\ \Longrightarrow\ 
\bigvee_{a\in\mathbb F_p}
 \bigl(\Theta_a(u+s)\wedge\Theta_a(v+s)\bigr)\right].
```

All coordinates are integer variables and $`n`$ is a free variable. There is no product of variables: $`b_jn`$ is multiplication by a fixed integer. Order anchors by

```math
u\prec v
\quad\Longleftrightarrow\quad
(|u_1|+|u_2|,u_1,u_2)
 <_{\mathrm{lex}}
(|v_1|+|v_2|,v_1,v_2).
```

Absolute values and this order expand into finite linear formulas. Every nonempty subset of $`\mathbb Z^2`$ has a least element for this order: first minimize the nonnegative integer norm, then minimize in the finite norm shell. Put

```math
\mathcal R(n,u)\ \Longleftrightarrow\ 
n\ge1\ \wedge\ 
\neg\exists v\,[E(n,u,v)\wedge v\prec u].
```

For every positive $`n`$, this selects exactly one anchor from each pattern class, even if that class is infinite. Thus

```math
R_n=\{u:\mathcal R(n,u)\},\qquad |R_n|=f(n)<\infty.
```

This argument proves finiteness before any polyhedral decomposition. It also covers empty windows, whose unique class has least anchor zero.

## 3. Finite elimination, with an explicit algorithm

We give the elementary elimination step to make “compute the representative cells” an actual finite prescription.

Reduce universal quantifiers by negation and eliminate innermost existential quantifiers successively. Put a quantifier-free formula under one existential quantifier in disjunctive normal form. An individual conjunction consists of integer linear inequalities, divisibilities and nondivisibilities. Equalities are two inequalities; the negation of $`L\le c`$ is $`-L\le-c-1`$. Clauses independent of the variable remain untouched.

For a variable $`x`$ to eliminate, let $`L`$ be the least common multiple of the absolute nonzero coefficients of $`x`$ in the clause, taking $`L=1`$ if there are none. Set $`X=Lx`$ and add $`L\mid X`$. In a term with coefficient $`a\ne0`$, multiply by $`L/|a|`$. An inequality becomes a lower or upper bound on $`X`$ with coefficient one. A divisibility transforms as

```math
d\mid (ax+t(y))
\quad\Longleftrightarrow\quad
d\,\frac{L}{|a|}
 \mid \left(\mathrm{sgn}(a)X+\frac{L}{|a|}t(y)\right),
\quad\text{subject to }L\mid X.
```

This equivalence holds because $`d\mid z`$ if and only if $`dc\mid cz`$ for positive integer $`c`$. The same equivalence applies to nondivisibility.

Let $`Q`$ be the least common multiple of all moduli in the transformed clause, including $`L`$. At fixed remaining variables, every divisibility literal is invariant under $`X\mapsto X+Q`$. If lower bounds occur, split into the finitely many choices of a largest lower bound $`l_j(y)`$, including the inequalities $`l_j\ge l_h`$ for all $`h`$. Test the finitely many values

```math
X=l_j(y)+h,\qquad h=0,\ldots,Q-1.
```

If a feasible $`X`$ exists, the smallest number in its residue class above $`l_j`$ is in this list; it is no larger than the original witness, so it respects every upper bound. Conversely any successful test supplies a witness. If there are no lower bounds but there are upper bounds, choose a smallest upper bound $`u_j`$ and test $`X=u_j-h`$ for the same range of $`h`$. If there are no bounds, test $`X=h`$.

Replace the existential clause by the corresponding finite disjunction of tests. The result again uses only integer linear inequalities and divisibilities. All substitutions multiply variables only by constants. Repeating this finite procedure eliminates every bound variable of $`\mathcal R`$. No bound on a pattern class or a representative coordinate is assumed. A coefficient-free clause simply retains its independent conditions; its auxiliary constraint has $`L=1`$.

Denote the resulting quantifier-free Boolean formula by $`F(n,u)`$. Write its inequality atoms as

```math
a_h\cdot u\le b_hn+c_h,
\qquad a_h\in\mathbb Z^2,\quad b_h,c_h\in\mathbb Z,
```

and its divisibility atoms as $`d\mid(\alpha\cdot u+bn+c)`$. Its equivalence to $`\mathcal R`$ has just been proved constructively. Different elimination orders can give different, equally valid, periods and cells.

## 4. The cells and the two explicit moduli

Let $`q`$ be the least common multiple of all divisibility moduli in $`F`$, or one if there are none. For every nonzero inequality normal $`a_h`$, set

```math
g_h=\gcd(|a_{h,1}|,|a_{h,2}|),\qquad
G=\mathrm{lcm}_h g_h,\qquad
M_\beta=qG,
```

using $`G=1`$ if the list is empty. Fix $`0\le r\lt M_\beta`$, a coordinate residue $`\rho\in\{0,\ldots,q-1\}^2`$, and substitute

```math
n=M_\beta k+r,\qquad u=qz+\rho,\qquad z\in\mathbb Z^2.
```

Every divisibility atom now has a constant truth value. An inequality with zero anchor normal also has a constant truth value for all sufficiently large $`k`$; determine it from its affine right side. For the remaining inequality atoms, enumerate all truth assignments $`\sigma`$ for which $`F`$ is true. A false inequality is replaced by its integer complement, so write every chosen condition again as

```math
a\cdot u\le bn+c.
```

The normal of a complement is $`-a_h`$, with the same $`g_h`$. With $`g=\gcd(|a_1|,|a_2|)`$, this condition on integer $`z`$ is exactly

```math
\nu\cdot z\le s k+t,\qquad
\nu=\frac ag,\quad s=\frac{bG}{g}\in\mathbb Z,\quad
t=\left\lfloor\frac{br+c-a\cdot\rho}{qg}\right\rfloor\in\mathbb Z.
```

Let $`C_{r,\rho,\sigma}(k)`$ be the closed real polygonal set defined by these conditions. Its normals are primitive integral and all its support constants are integers. Cells are used only for their integer points. Different accepted truth assignments cannot count the same integer $`z`$, and different $`\rho`$'s cannot count the same anchor.

**Discarding unbounded cells is justified.** Their horizontal recession cone is

```math
K=\{w\in\mathbb R^2:\nu\cdot w\le0\text{ for every cell inequality}\},
```

independent of $`k`$. If $`K\ne\{0\}`$, it contains a nonzero integer vector, since it is a rational polyhedral cone. Any integer point of any such cell would then generate infinitely many integer points by adding nonnegative multiples of that vector. All would belong to $`R_n`$, contrary to its proved finiteness. Such a cell has no integer points at every relevant positive $`n`$, and can be discarded even if its real set is nonempty. No feasibility sampling is used in this step.

Every remaining nonempty fiber is a bounded rational polytope. For every pair of independent normals solve the two boundary equations. Its candidate vertex has the form

```math
p_j(k)=k\alpha_j+\delta_j,\qquad
\alpha_j,\delta_j\in\mathbb Q^2.
```

Feasibility of each candidate against every other inequality is an affine inequality in $`k`$. Hence feasibility, vertex coincidences, cyclic order and dimension stabilize after an effectively computable threshold. One can calculate that threshold from the finitely many rational crossings; all statements below concern the stable range. Cells may be two-dimensional polygons, segments, points or empty. They must not all be assumed to have positive-area limits.

Let $`\mathcal V`$ be the list of primitive normals of the retained conditions, taking a union over all residues and cells. Define the entirely explicit integers

```math
\Delta=\mathrm{lcm}\left(\{1\}\cup
\{|\det(\nu,\nu')|:\nu,\nu'\in\mathcal V,\ \det(\nu,\nu')\ne0\}\right),
```

```math
J=\mathrm{lcm}\left(\{1\}\cup
\{|\nu_2|:\nu=(\nu_1,\nu_2)\in\mathcal V,\ \nu_2\ne0\}\right),
\qquad \boxed{\,M=M_\beta\Delta J\,}.
```

An empty retained list gives $`\Delta=J=1`$. These are sufficient periods, not denominators of the original polygon alone. They are finite functions of its facet data and the full finite star tables, through the explicit elimination.

## 5. The geometric coefficient of a cell

For a two-dimensional cell in its stable range, list its vertices counterclockwise as $`p_j(k)=k\alpha_j+\delta_j`$, with cyclic subscripts. Define

```math
a_2(C)=\frac12\sum_j\det(\alpha_j,\alpha_{j+1}),
```

```math
a_1(C)=\frac12\sum_j
 \left(\det(\alpha_j,\delta_{j+1})
       +\det(\delta_j,\alpha_{j+1})\right).
```

These are exactly the quadratic and linear coefficients of its area, by the shoelace formula. For every edge choose its primitive integral tangent $`w_j`$, oriented from $`p_j`$ to $`p_{j+1}`$, and write

```math
p_{j+1}(k)-p_j(k)=(\lambda_jk+\mu_j)w_j,\qquad
\lambda_j\ge0.
```

The numbers are rational and found by dividing any nonzero coordinate of $`w_j`$. Edges of bounded length have $`\lambda_j=0`$. Put

```math
b(C)=a_1(C)+\frac12\sum_j\lambda_j
\qquad\text{for a two-dimensional cell}.
```

For a one-dimensional cell use its primitive tangent $`w`$, orient its endpoints so their difference is $`(\lambda k+\mu)w`$ with $`\lambda\ge0`$, and put $`a_2(C)=0,\ b(C)=\lambda`$. For points and empty cells put both numbers zero. Finally define

```math
A_r=\sum_{\rho,\sigma}a_2(C_{r,\rho,\sigma}),\qquad
B_r=\sum_{\rho,\sigma}b(C_{r,\rho,\sigma}),
```

with discarded cells omitted. These finite rational sums are the promised geometric formula.

### 5.1 Integer support lines, including thin cells

We prove

```math
|C(k)\cap\mathbb Z^2|
 =a_2(C)k^2+b(C)k+O(1).
```

For a two-dimensional cell, each long facet lies on a line with primitive integer normal and integer support constant. That line contains lattice points spaced by its primitive tangent. Select lattice points on each long facet within one primitive step of each endpoint. The cell has fixed facet directions. Consecutive bounded facets and the ends of adjacent long facets group into corner neighborhoods of bounded diameter: their vertices have equal $`\alpha_j`$'s, so their mutual differences are bounded. Removing those corner neighborhoods leaves precisely the long facets, up to bounded endpoint pieces.

The convex hull of the selected lattice points loses only bounded area in those neighborhoods. Consequently, if $`H_k`$ is the integer hull of $`C(k)`$, then

```math
\mathrm{area}(H_k)=\mathrm{area}(C(k))+O(1),\qquad
|\partial H_k\cap\mathbb Z^2|
 =\sum_j(\lambda_jk+\mu_j)+O(1).
```

If the limiting polygon has positive area, $`H_k`$ is two-dimensional for large $`k`$. If its limit is a segment and some edges grow, the two long sides are on distinct parallel integer support lines; their separation in their common primitive-normal coordinate is a positive integer. The selected points on both sides again give a two-dimensional hull for large $`k`$. If there are no long edges, the whole cell has bounded diameter and its lattice count, area and perimeter are $`O(1)`$, which already proves the desired assertion. These observations cover cells of constant transverse width; omitting them would in general lose part of $`\beta`$.

In the two-dimensional hull cases, Pick's theorem gives

```math
|C(k)\cap\mathbb Z^2|
=\mathrm{area}(C(k))
 +\frac12\sum_j(\lambda_jk+\mu_j)+O(1).
```

Here the integer hull has exactly the same lattice points as the original convex cell. The coefficient is therefore $`b(C)`$.

For a segment, its containing line is an active primitive integer support line, hence contains a full one-dimensional lattice. An interval of lattice length $`\lambda k+\mu`$ contains that length plus $`O(1)`$ lattice points, independent of its phase. For points the count is at most one. This proves every case of the formula.

### 5.2 Summing and returning to the original scale

The cells enumerate the representative set exactly on the relevant residue. Their lattice sets are disjoint, so

```math
f(M_\beta k+r)=A_rk^2+B_rk+O(1).
```

The reviewed area theorem on this same arithmetic progression gives $`A_r=M_\beta^2\Lambda`$. Substitution of $`k=(n-r)/M_\beta`$ gives the displayed formula for $`\beta_r`$. This identifies every residue separately, not an average over residues.

## 6. Why the displayed full modulus works

We prove eventual exact polynomiality on each residue modulo $`\Delta J`$ in the variable $`k`$. This proves the stated full period without using an unspecified modulus from a counting theorem.

For every candidate vertex, Cramer's rule shows $`\Delta\alpha_j\in\mathbb Z^2`$. Fix $`h`$ modulo $`\Delta J`$ and set $`k=\Delta Jt+h`$. Each vertex coordinate is then affine in $`t`$ with integer slope divisible by $`J`$. Its floor and ceiling are affine integer functions of $`t`$, with slopes divisible by $`J`$.

Partition the real horizontal projection of a cell at the distinct vertex first coordinates. Their order has stabilized. On a slab use the half-open convention $`[x_l(k),x_{l+1}(k))`$, and handle the last endpoint separately. The integer columns of a slab range from $`\lceil x_l(k)\rceil`$ to $`\lceil x_{l+1}(k)\rceil-1`$. At every such column, its lower and upper boundaries are fixed affine rational functions supplied by two fixed facets. Its count is

```math
\lfloor U_k(x)\rfloor-\lceil L_k(x)\rceil+1.
```

Because $`L_k(x)\le U_k(x)`$ on the slab, this expression is nonnegative even in a thin cell or on a nonvertical segment. A boundary function has the form

```math
\frac{sk+t_0-\nu_1x}{\nu_2},
\qquad \nu_2\ne0,
```

with all coefficients integral. For the fixed progression in $`k`$, $`sk`$ has constant residue modulo every $`|\nu_2|`$, since $`|\nu_2|\mid J`$. Split the columns into their residues modulo $`J`$. The two endpoints of each resulting arithmetic progression are affine in $`t`$, their residues being fixed. After writing $`x=Jz+a`$, each floor or ceiling is an integer affine function of $`z,t`$. Summing it over the affine integer bounds for $`z`$ is an exact polynomial in $`t`$ of degree at most two. Slabs with no columns are either eventually empty or handled by the same zero-length sum.

At the final endpoint, the vertex first coordinate either is integral throughout this progression or never is: its slope is integral and its fractional part constant. If integral, the one-column formula is affine in $`t`$ because its coefficients of $`t`$ are divisible by all relevant denominators. The same argument treats a vertical segment. A point either is always integral on the progression or never is; both its coordinates have integral slopes.

These arguments also apply to one-dimensional and bounded cells. The finitely many slabs give an exact polynomial for each cell on each progression, for all sufficiently large $`t`$. Summing cells gives such a polynomial for $`f(n)`$ on every residue modulo $`M_\beta\Delta J`$. Section 5 identifies its quadratic and linear coefficients. Subtracting those leaves a constant on each residue, which is the asserted $`\gamma`$. This completes the proof.

A stable threshold can be obtained from the vertex-feasibility and slab-order tests in Section 4. Exact residue polynomials, and hence the constants, can also be calculated by the column sums above. Neither requires pattern enumeration at sampled scales.

## 7. Sanity checks and exact limits of the result

For the full-support cross $`\theta(x,y)=\mathbf1[y=0]+\mathbf1[x=0]`$ over $`\mathbb F_3`$ and $`P=[0,1/2]\times[0,1]`$, every pattern is specified by an absent-or-present horizontal row and an absent-or-present vertical column. For large $`n`$, the single color-two crossing distinguishes the two-present case. Thus a bijective descriptor set is the integer rectangle

```math
-1\le x,\quad 2x\le n,\qquad -1\le y\le n.
```

This representative family need not be the norm-minimal one; any explicitly proved definable transversal or bijective descriptor set can replace Section 2. Here $`q=1,G=2,M_\beta=2`$. With $`n=2k+r`$, the real rectangle has side lengths $`k+1`$ and $`2k+r+1`$. Its area linear coefficient is $`r+3`$, and half its perimeter coefficient is $`3`$. Thus $`A_r=2,B_r=r+6`$ and

```math
\beta_r=(r+6)/2-r=3-r/2.
```

This agrees with the exact count $`(\lfloor n/2\rfloor+2)(n+2)`$. In this example $`\Delta=J=1`$, so $`M=2`$.

For the sparse binary cross with supports $`y=0,\ x\equiv0\pmod2`$ and $`x=1,\ y\equiv0\pmod2`$, the exact count is $`AB+2A+2B+1`$, where $`A=\lfloor n/2\rfloor+1,\ B=n+1`$. A disjoint descriptor decomposition consists of the $`A`$-by-$`B`$ bulk rectangle, two pure-row phase columns, two pure-column phase rows, and one zero descriptor. Its formula gives $`\beta_0=9/2,\ \beta_1=4`$. The constant-width two-column and two-row cells are essential in this calculation.

An integer-support triangle $`x\ge0,\ y\ge0,\ 2x+3y\le k`$ illustrates why a period of $`\beta`$ is not automatically a period of the constant remainder. Its normals are primitive, so $`M_\beta=1`$; its area coefficient is $`1/12`$, its half-boundary coefficient is $`1/2`$, and the exact bounded remainder is periodic. The prescription gives the sufficient full period $`18`$, not necessarily the minimal period $`6`$.

The diagnostics test these coefficient calculations, bounded-width cells, segments, empty unbounded cells, the finite elimination step, and exact lattice counts on further fixed rational cells. They are finite cross-checks of the proof, not its premises. No complete production implementation of the whole star-to-eliminated-formula compiler is claimed.

Only a reduced data set is obstructed: the directions and window alone cannot determine $`\beta`$, as the full and sparse crosses demonstrate. This is not an obstruction to the complete finite-table formula proved here. This elimination construction does not give a practically useful size bound or a compact formula on the original direction bands alone. The separate [exact-period construction](minimal-periods-and-constant-term.md) decides the minimal modulus; the [integral-window theorem](integral-windows.md) gives constant $`\beta`$ and a direct finite boundary-profile formula for integral convex polygons.

## 8. Review and attribution boundary

The finite color description and norm-first transversal were already present in the independently reviewed eventual-quasipolynomial argument. The present contribution is the fully specified elimination-to-polygon recipe, the treatment of thin and lower-dimensional cells, the geometric coefficient extraction, and the explicit sufficient moduli. The elimination is a classical Presburger method and Pick's theorem is classical; the result is their application to these representatives, not a claim to new logic or new lattice-point theory.

This supplementary proof passed one fresh-context GPT-6.1 Sol adversarial review on 25 October 2026, with no fatal or proof-repair issue found. No external mathematician or proof assistant has certified it. The reviewed leading-term theorem is used only to identify the already geometrically computed $`A_r/M_\beta^2`$ with the original direction-band area. The proof of the sufficient modulus itself does not invoke the unspecified period in Bogart–Goodrick–Woods.

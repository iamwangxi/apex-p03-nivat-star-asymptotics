# Second-order pattern complexity for periodic colored lattice lines

The results in this module received a separate-context adversarial review on 24 October 2026. Theorem numbers and equation numbers on this page are local to this module. The arbitrary full-support result is Theorem 9; it does not rely on the supplementary general coefficient construction, separately reviewed on 25 October 2026.

## 1. Scope and results

For a finite set $`S\subset\mathbb Z^2`$, write
```math
\mathcal P_\theta(S)=\{(\theta(u+z))_{z\in S}:u\in\mathbb Z^2\},
\qquad P_\theta(S)=|\mathcal P_\theta(S)|.
\qquad\text{(1)}
```
Throughout, $`P`$ is a positive-area rational convex polygon and $`S_n=nP\cap\mathbb Z^2`$. Integer $`n`$ tends to infinity. A modulus $`M`$ may be any positive integer making every vertex of $`MP`$ integral.

The [all-star area theorem](all-star-extension.md) gives the exact leading term
```math
P_\theta(S_n)=\mathrm{area}(\mathcal M(P))n^2+O(n),
\quad
\mathcal M(P)=\bigcup_{i\lt j}\bigl(B_i(P)\cap B_j(P)\bigr),
\qquad\text{(2)}
```
where $`B_i(P)`$ is the strip with the same projection range as $`P`$ in the normal to the primitive direction $`v_i`$. We use (2) only for comparison; the new results below have independent counting proofs.

This note proves four assertions.

**Periodic-coefficient counterexample.** A Case B star on a rational rectangle has different even and odd linear coefficients. A constant $`\beta`$ is therefore impossible for the rational-window problem as stated.
**Two-line result.** For **two full-support, nonzero periodically colored primitive affine lattice lines**, arbitrary rational polygons admit
   ```math
P_\theta(S_n)=\Lambda(P)n^2+\beta_{n\bmod M}(\theta,P)n+O(1).
   \qquad\text{(3)}
```
   The coefficients are computed from finitely many lattice-feasibility tests in vertex tangent cones. All pattern collisions reduce to a uniformly bounded, explicitly enumerable list; no collision-finiteness hypothesis is used.
**Long-extreme-face result.** For **any finite number of such lines**, (3) holds whenever both extreme faces of $`P`$ parallel to each line direction are nondegenerate edges. In this class the complete representation count is eventually collision-free, and $`\beta_r`$ has a finite formula in facet lattice lengths and support displacements. There is no restriction on pair determinants or on multiple intersections of the fixed lines.
**Arbitrary-window result.** For **any finite number of full-support periodically colored lines**, the long-extreme-face window restriction in the preceding result can be removed, including intersecting supports. For every rational polygon, a complete classification reduces all linear-sized collisions to singleton reassignments near vertices. A finite interval-and-phase integral computes their exact linear correction. The proof first treats pairwise disjoint lattice supports and then removes that assumption using a bounded list of short-section intersection anchors. This class contains the earlier three-line Case B configuration and allows arbitrary many nonconcurrent lines, nontrivial color periods, multiple intersections and cancellation.

All three positive results give a constant linear coefficient for integral $`P`$. They do **not** prove this for all CN stars. Zero entries in a component word, finite-width strips, finite defects or different tails in a component, and nonzero periodic backgrounds remain outside the positive theorems. The residue counterexample itself is allowed to have zero entries.

## 2. A rational-window Case B counterexample

Consider the binary configuration
```math
c(x,y)=\mathbf1_{\{y=0,\ x\equiv0\pmod2\}}
       +\mathbf1_{\{x=1,\ y\equiv0\pmod2\}}
       \pmod2.
\qquad\text{(4)}
```
The two supports are disjoint. Each is periodic along its supporting primitive lattice line and has zero tails. The binary color indicators are $`c`$ and $`1-c`$, both annihilated in characteristic zero by the product of the two tangential differences, so this is a CN Case B star. The underlying configuration is an existing sparse orthogonal example; the rational-window consequence below is the new use of it.

**Theorem 1 (a genuinely periodic linear coefficient).** Put
$`P=[0,\tfrac12]\times[0,1]`$. For every $`n\ge8`$, with $`r=n\bmod2`$,
```math
P_c(S_n)=\frac12n^2+
 \left(\frac92-\frac r2\right)n+6-\frac{3r}{2}.
\qquad\text{(5)}
```
Here $`\mathrm{area}(\mathcal M(P))=1/2`$. Thus the even and odd linear coefficients are $`9/2`$ and $`4`$, respectively.

**Proof.** First use the rectangular site set
$`R_{A,B}=\{0,\ldots,A-1\}\times\{0,\ldots,B-1\}`$, with $`A,B\ge5`$. A translated horizontal component, if visible, has at least two sites on one row; similarly a visible vertical component has at least two sites on one column. In a two-component pattern, that row and column are intrinsic: every other row or column has at most one nonzero site. They determine the translation uniquely. There are exactly $`AB`$ such translations.

Pure horizontal patterns have $`B`$ possible rows and two phases; pure vertical patterns have $`A`$ possible columns and two phases. Every phase is realizable by translating sufficiently far along the corresponding line, while keeping its parity. These patterns, the two-component patterns, and the zero pattern are disjoint. Hence
```math
P_c(R_{A,B})=AB+2A+2B+1.
\qquad\text{(6)}
```
For $`S_n`$, take $`A=\lfloor n/2\rfloor+1`$, $`B=n+1`$, and expand (6). The two normal strips intersect in the rectangle $`P`$, giving the asserted leading area. ∎

In particular, neither choice of a constant $`\beta`$ can make the error in (5) bounded on both parities. This settles the constant-versus-periodic alternative for **rational** polygons. It does not provide an integral-polygon counterexample.

There is also a counterexample inside the full-support class used below. Over $`\mathbb F_3`$, let
$`\theta(x,y)=\mathbf1_{y=0}+\mathbf1_{x=0}`$. For a rectangular site set with at least two rows and columns, the unique color $`2`$ identifies every two-line pattern, and
```math
P_\theta(R_{A,B})=(A+1)(B+1).
```
On the same half-width polygon this gives, exactly,
```math
P_\theta(S_n)=\frac12 n^2+
 \left(3-\frac r2\right)n+4-r\qquad(n\ge2).
\qquad\text{(7)}
```

## 3. Exact projection counts from finitely many cones

Fix a primitive vector $`v\in\mathbb Z^2`$, choose $`e\in\mathbb Z^2`$ with $`\det(v,e)=1`$, and put
```math
\pi(z)=\det(v,z),\qquad s(z)=\det(z,e).
```
Thus $`z=s(z)v+\pi(z)e`$. Define
```math
\ell=\min_P\pi,\quad u=\max_P\pi,\quad W=u-\ell>0,
\quad Q_v(n)=\{\pi(z):z\in S_n\}.
\qquad\text{(8)}
```
This section applies to any rational polygon, independently of any configuration.

### 3.1 A uniform bound on short sections

Choose a rational interior point $`q_0`$ of $`P`$, put $`a=\pi(q_0)`$, and let $`h_0>0`$ be the length, in the coordinate $`s`$, of the section $`P\cap\{\pi=a\}`$. Set
```math
c_v=\frac{h_0}{\max(a-\ell,u-a)}>0,
\qquad C_v(K)=\left\lceil\frac {K+1}{c_v}\right\rceil+2
\quad(K\ge1).
\qquad\text{(9)}
```

**Lemma 2 (only finitely many short heights).** If an integer $`t`$ satisfies
```math
n\ell+C_v(K)\le t\le nu-C_v(K),
\qquad\text{(10)}
```
then $`S_n\cap\{\pi=t\}`$ has at least $`K+1`$ consecutive sites in direction $`v`$. Consequently the number of nonempty sections of length less than $`K`$, or at most $`K`$, is bounded by $`2(C_v(K)+1)`$, uniformly in $`n`$.

**Proof.** The real section-length function $`h(t)`$ on $`[\ell,u]`$ is concave. This also follows directly by taking convex combinations of the central section with a point on either extreme face. On the left and right of $`a`$, respectively,
```math
h(t)\ge h_0\frac{t-\ell}{a-\ell},\qquad
h(t)\ge h_0\frac{u-t}{u-a}.
```
In particular $`h(t)\ge c_v\min(t-\ell,u-t)`$. The section of $`nP`$ at height $`t`$ has real length at least
$`c_v\min(t-n\ell,nu-t)`$. Under (10) this is strictly greater than $`K+1`$; a closed real interval of that length contains at least $`K+1`$ integers. For integer $`t`$, the section's lattice points are exactly the consecutive integer $`s`$-coordinates. All shorter sections therefore occur within a fixed distance $`C_v(K)`$ of one of the two extreme heights. Counting integer heights in those two intervals gives the stated safe bound, including sections of length exactly $`K`$. ∎

The extra margin in (9) avoids endpoint conventions. Missing projection heights in particular occur only at bounded distance from $`n\ell`$ or $`nu`$.

### 3.2 Finite cone certificates

Write a rational half-plane description
```math
P=\{x:A_fx\le b_f\text{ for all facets }f\}.
\qquad\text{(11)}
```
Choose vertices $`q^-`$ and $`q^+`$ attaining $`\ell`$ and $`u`$. Either endpoint of an extreme edge may be chosen. Let
```math
T_q=\{w:A_fw\le0\text{ whenever }A_fq=b_f\}
\qquad\text{(12)}
```
be the tangent cone at $`q`$. In this subsection put $`C=C_v(1)`$.

For $`0\le r\lt M`$, define two finite sets of nonnegative integers:
```math
\begin{split}
D_v^-(r)=\{0\le h\le C:
 &\nexists z\in\mathbb Z^2,
 \ \pi(z)=\lceil r\ell\rceil+h,
 \ z-rq^-\in T_{q^-}\},\\ 
D_v^+(r)=\{0\le h\le C:
 &\nexists z\in\mathbb Z^2,
 \ \pi(z)=\lfloor ru\rfloor-h,
 \ z-rq^+\in T_{q^+}\}.
\end{split}
\qquad\text{(13)}
```
All tests in (13) are explicit one-variable integer interval tests. Indeed, write $`z=sv+te`$. The incident facets give
```math
(A_fv)s\le rb_f-tA_fe,\qquad s\in\mathbb Z.
\qquad\text{(14)}
```
Intersect their rational upper and lower bounds, including the possibility of an unbounded interval. Feasibility is exactly the existence of an integer in this interval.

**Lemma 3 (eventually exact projection formula).** There is a computable $`N_v`$ such that, for $`n\ge N_v`$, $`r=n\bmod M`$, $`\alpha_n=\lceil n\ell\rceil`$, and $`\omega_n=\lfloor nu\rfloor`$,
```math
Q_v(n)=([\alpha_n,\omega_n]\cap\mathbb Z)
 \setminus\bigl((\alpha_n+D_v^-(r))
             \cup(\omega_n-D_v^+(r))\bigr).
\qquad\text{(15)}
```
The two deleted endpoint lists are disjoint for this threshold. In particular
```math
|Q_v(n)|=nW+b_v(r),
\qquad\text{(16)}
```
where the following is a finite geometric formula:
```math
b_v(r)=\lfloor ru\rfloor-\lceil r\ell\rceil+1-rW
       -|D_v^-(r)|-|D_v^+(r)|.
\qquad\text{(17)}
```

**Proof.** Put $`n=Mk+r`$. For a lower endpoint height, translating a lattice point by $`-kMq^-`$ gives a lattice point at height $`\lceil r\ell\rceil+h`$, because $`Mq^-`$ is integral and $`\pi(q^-)=\ell`$. Necessity of the lower cone test follows from $`nP\subset nq^-+T_{q^-}`$. The upper endpoint is identical.

For sufficiency, choose one lattice witness $`z`$ for any feasible test in (13). Every incident facet already satisfies $`A_fz\le rb_f`$. For a nonincident facet,
```math
A_f(kMq+z)\le(Mk+r)b_f
\quad\Longleftrightarrow\quad
k\ge\frac{A_fz-rb_f}{M(b_f-A_fq)}.
\qquad\text{(18)}
```
The denominator is positive. Taking the maximum of zero and the ceilings in (18) gives an explicit stable value $`k_0`$ for this witness. There are finitely many residues, sides and offsets, so their $`k_0`$'s have a finite computable maximum. Thus every feasible cone height is realized in $`nP`$ for all sufficiently large $`n`$ in its residue class.

By Lemma 2 all other heights between the extremes occur. Increase the threshold until the endpoint lists are disjoint. For example, the maximum of $`M\max k_0`$, $`1`$, and $`\lceil(2C+4)/W\rceil+2`$ is sufficient. This proves (15); counting it and using integrality of $`M\ell,Mu`$ proves (16)–(17). ∎

No sequence stabilization is inferred from data: (18) proves it. If $`P`$ is integral, take $`M=1`$ in (17), so $`b_v`$ is constant. The sets (13) are independent of the choice of an extreme vertex in the eventual formula, since (15) is exact.

## 4. Two periodically colored full lattice lines

Let $`v_1,v_2`$ be nonparallel primitive integer vectors. For each $`i`$, choose $`e_i`$ with $`\det(v_i,e_i)=1`$, and write $`\pi_i(z)=\det(v_i,z)`$, $`s_i(z)=\det(z,e_i)`$. Fix integers $`a_i`$ and words
```math
w_i:\mathbb Z\longrightarrow\mathbb F_p^\times
\quad\text{of minimal periods }\kappa_i.
```
All word entries must be nonzero. Define
```math
F_i(z)=\begin{cases}w_i(s_i(z)),&\pi_i(z)=a_i,\\0,&\text{otherwise},\end{cases}
\qquad\theta=F_1+F_2\pmod p,
\quad D=|\det(v_1,v_2)|.
\qquad\text{(19)}
```
This is a CN star with zero tails. The affine lines may meet at a lattice point, where their colors may cancel, or may have no common lattice point. No unimodularity assumption is made.

Use $`W_i, Q_i(n),b_i(r)`$ from Section 3. Let $`N(n)`$ be the number of lattice anchors $`u`$ for which both geometric sections
```math
S_n\cap\{\pi_i=a_i-\pi_i(u)\},\quad i=1,2,
\qquad\text{(20)}
```
are nonempty. Geometric activity is defined before adding colors.

**Theorem 4 (two-line coefficient and complete finite collision reduction).** For (19) and every positive-area rational polygon $`P`$,
```math
P_\theta(S_n)=\frac{W_1W_2}{D}n^2+
 \left(\frac{W_1b_2(r)+W_2b_1(r)}{D}
             +\kappa_1W_1+\kappa_2W_2\right)n+O(1),
\qquad r=n\bmod M.
\qquad\text{(21)}
```
Here $`W_1W_2/D=\mathrm{area}(\mathcal M(P))`$. The error is uniformly bounded, and a bounded list computes its collision correction exactly. In particular integral $`P`$ has a constant linear coefficient, even when $`D>1`$, the words have nontrivial periods, or the intersection color vanishes.

### 4.1 Counting geometric two-line anchors

Complete $`v_1`$ to the basis $`(v_1,e_1)`$ and write
$`v_2=A v_1+d e_1`$, where $`d=\det(v_1,v_2)`$. Primitivity of $`v_2`$ gives $`\gcd(A,d)=1`$. The two offsets $`t_i=a_i-\pi_i(u)`$ determine at most one anchor, and this anchor is integral exactly when
```math
t_2\equiv A t_1+a_2-Aa_1\pmod D.
\qquad\text{(22)}
```
Indeed, for $`u=sv_1+te_1`$, $`\pi_1(u)=t`$ and $`\pi_2(u)=At-ds`$. Thus (22) is both necessary and sufficient.

For large $`n`$, each $`Q_i(n)`$ is an integer interval minus a bounded number of endpoint holes. Put $`T_i=|Q_i(n)|`$. Its count in any residue class modulo $`D`$ is $`T_i/D+\epsilon_i`$, where the $`\epsilon_i`$'s are uniformly bounded and sum to zero. Because multiplication by $`A`$ permutes residues, (22) gives
```math
N(n)=\frac{T_1T_2}{D}+O(1).
\qquad\text{(23)}
```
There is no unaccounted linear residue error here: the two possible cross terms vanish by their zero sums. If $`h_i=\max_r(|D_i^-(r)|+|D_i^+(r)|)`$, a safe bound for the error in (23) is $`D(h_1+1)(h_2+1)`$.

### 4.2 A long active section determines a two-line pattern uniquely

**Lemma 5 (long-section uniqueness, including cancellation).** A two-line anchor whose $`i`$-section has at least $`3\kappa_i`$ sites gives a pattern realized by no other anchor and by no pure-line representation.

**Proof.** At most one site of its $`i`$-line can be changed or canceled by the other line. Hence at least two nonzero sites remain on that line. In any competing pattern the component of the other direction can cover at most one of these sites; therefore its $`i`$-line must be the same relative affine line.

If the original pattern has any nonzero site off this line, the other relative line must pass through that site. The two nonparallel relative lines determine both transverse offsets and thus determine the anchor uniquely.

Otherwise the original other-line section is one site, at the intersection of the two relative lines. This requires the fixed affine lines to have an integral intersection $`z_*`$. Its contribution is a **nonzero additive edit** of the $`i`$-word at the relative site $`z_*-u`$, even if the resulting color at that site is zero.

A competitor with the same $`i`$-line cannot have an off-line site, so its other section is either empty or one site at its own relative intersection. Along the common $`i`$-section, compare the two unedited baseline words. Their difference is $`\kappa_i`$-periodic. Equality of the edited patterns would make this difference supported on at most two sites. A nonzero $`\kappa_i`$-periodic word on at least $`3\kappa_i`$ consecutive sites has at least three nonzero sites. Thus the baseline words agree. Their nonzero edits must then occur at the same site; a pure competitor is impossible, and $`z_*-u=z_*-u'`$ forces $`u=u'`$. ∎

### 4.3 The entire exceptional list

Use the following representation list, with multiplicities:

- one entry for each geometric two-line anchor counted by $`N(n)`$;
- for each $`i`$, each height $`t\in Q_i(n)`$, and each of the $`\kappa_i`$ tangent phases, the pure $`i`$-word on that section;
- one zero entry.

It represents every global pattern. Each pure phase is globally realizable: choose an anchor with the required transverse offset and phase and translate it by multiples of $`\kappa_i v_i`$. The other line's transverse offset moves strictly monotonically and eventually leaves its bounded visible range. A zero pattern is realizable outside both visible strips.

Let $`H=3(\kappa_1+\kappa_2)`$. Form the sublist $`E_n`$ consisting of:

1. two-line anchors with both section lengths $`<3\kappa_i`$;
2. pure-line entries whose section has at most $`H`$ sites;
3. the zero entry.

This is a complete collision list. Lemma 5 makes all other two-line entries globally unique. A pure word on a section of length at least $`\kappa_i`$ has distinct phases by minimality of the period. Pure patterns on different parallel lines cannot agree; pure patterns in different directions can agree only if their support is a singleton. Finally, a bad two-line entry has support size less than $`H`$, so any equal pure entry belongs to item 2. These observations show that every entry outside $`E_n`$ is unique and has a pattern absent from $`E_n`$.

To see finiteness explicitly, write $`R_i(K)=2(C_{v_i}(K)+1)`$. Lemma 2 gives
```math
|E_n|\le R_1(3\kappa_1)R_2(3\kappa_2)
       +\kappa_1R_1(H)+\kappa_2R_2(H)+1.
\qquad\text{(24)}
```
Each pair of short heights gives at most one anchor by (22). For enumeration, inspect these endpoint height bands, find their section intervals by (11), solve (22), and evaluate the finite words. This is an algorithm with a number of candidates bounded independently of $`n`$, not a conjectural list of collision types.

If $`\chi_n=|E_n|-|\{\text{patterns of entries of }E_n\}|`$, the following count is **exact**:
```math
P_\theta(S_n)=N(n)+\kappa_1T_1+\kappa_2T_2+1-\chi_n.
\qquad\text{(25)}
```
Formula (24) bounds $`\chi_n`$. Thus the boundary-shell anchors and the pure representations are counted exactly up to a uniformly bounded, exactly computable correction. There is no linear-sized collision family in this two-component full-support class.

### 4.4 Deriving the coefficient

Insert (16) and (23) into (25). The quadratic and linear terms are precisely (21). All remaining terms are bounded, using (24) and the finitely many $`b_i(r)`$'s. The map $`x\mapsto(\pi_1(x),\pi_2(x))`$ has determinant of absolute value $`D`$, and the intersection of the two normal strips maps to a rectangle of side lengths $`W_1,W_2`$. This proves the area statement and completes Theorem 4. ∎

### 4.5 An integral example where missing layers change the answer

Take the standard triangle $`\mathrm{conv}\{(0,0),(1,0),(0,1)\}`$, directions $`(1,2),(1,4)`$, levels $`1,2`$, constant words $`1,1`$ over $`\mathbb F_2`$. The two affine lattice lines are disjoint. The widths are $`3,5`$, and the determinant is $`2`$.

At the negative extreme, the first direction has offsets generated by $`2,3`$, with only the gap $`1`$. The second has offsets generated by $`4,5`$, with gaps $`1,2,3,6,7,11`$. At the positive extreme both semigroups include $`1`$, so there are no gaps. Equivalently, these lists follow directly from (13) at the two relevant vertices. Thus eventually
```math
T_1=3n,\qquad T_2=5n-5.
```
Theorem 4 gives
```math
P_\theta(S_n)=\frac{15}{2}n^2+\frac12n+O(1).
\qquad\text{(26)}
```
Using only continuous projection widths and omitting endpoint holes would give an incorrect linear coefficient. The checked constant remainder $`-36`$ at the sampled scales is not asserted here to hold for every large $`n`$.

## 5. Arbitrarily many directions with long extreme sections

We now allow $`m\ge2`$ pairwise nonparallel primitive vectors $`v_i`$, levels $`a_i\in\mathbb Z`$, and nonzero periodic words $`w_i`$ of minimal periods $`\kappa_i`$. Define $`F_i`$ as in (19) and put $`\theta=\sum_{i=1}^m F_i`$ over $`\mathbb F_p`$. Intersections, including simultaneous intersections and total color cancellation there, are allowed.

Assume the following geometric condition:
```math
P\cap\{\pi_i=\ell_i\}\quad\text{and}\quad
P\cap\{\pi_i=u_i\}
\text{ are nondegenerate line segments for every }i.
\qquad\text{(27)}
```
Thus the two extreme faces parallel to every component direction are edges. This includes every positive rational zonotope with a generator parallel to each $`v_i`$; it permits arbitrarily many directions and nonunimodular determinants. The polygon may have additional edge directions.

For $`I\subset\{1,\ldots,m\}`$, $`|I|\ge2`$, define
```math
K_I=\bigcap_{i\in I}\{x:-u_i\le\pi_i(x)\le-\ell_i\},
\qquad c_I=(-1)^{|I|}(|I|-1).
\qquad\text{(28)}
```
Each displayed inequality defines a nondegenerate facet of $`K_I`$: the corresponding edge of $`-P`$ lies on it and satisfies all the other inequalities. For its upper and lower facets let
```math
\tau_i^+(I)=\frac{\mathrm{length}(K_I\cap\{\pi_i=-\ell_i\})}{\|v_i\|},
\quad
\tau_i^-(I)=\frac{\mathrm{length}(K_I\cap\{\pi_i=-u_i\})}{\|v_i\|}.
\qquad\text{(29)}
```
These rational quantities are facet lengths measured in primitive lattice steps. Define the support displacements
```math
\delta_i^+(r)=a_i+r\ell_i-\lceil r\ell_i\rceil,
\qquad
\delta_i^-(r)=-a_i+\lfloor ru_i\rfloor-ru_i.
\qquad\text{(30)}
```

**Theorem 6 (a geometric coefficient for arbitrary many lines).** Under (27), (3) holds with $`\Lambda(P)=\mathrm{area}(\mathcal M(P))`$ and
```math
\begin{split}
\beta_r(\theta,P)=\sum_i\kappa_iW_i
 +\sum_{|I|\ge2}c_I\sum_{i\in I}\bigl[
 &\tau_i^+(I)(\delta_i^+(r)+\tfrac12)\\ 
 +&\tau_i^-(I)(\delta_i^-(r)+\tfrac12)\bigr].
\end{split}
\qquad\text{(31)}
```
For all sufficiently large $`n`$, there is also the exact formula
```math
P_\theta(S_n)=
 \sum_{|I|\ge2}c_I\,|K_I(n)\cap\mathbb Z^2|
 +\sum_i\kappa_i\bigl(\lfloor nu_i\rfloor-\lceil n\ell_i\rceil+1\bigr)+1,
\qquad\text{(32)}
```
where
```math
K_I(n)=\bigcap_{i\in I}
 \{x:a_i-\lfloor nu_i\rfloor\le\pi_i(x)
                         \le a_i-\lceil n\ell_i\rceil\}.
\qquad\text{(33)}
```
In particular, for integral $`P`$, $`M=1`$, $`\delta_i^+=a_i`$, $`\delta_i^-=-a_i`$, so (31) gives a constant linear coefficient. Word values apart from their minimal periods affect neither coefficient in this class.

### 5.1 Why no representation collision is omitted

Let $`\mu_i>0`$ be the smaller of the two extreme real section lengths of $`P`$ in direction $`v_i`$, measured in primitive steps. Concavity of section length gives a lower bound $`n\mu_i`$ for every real section of $`nP`$. Thus, for every integer height between its projection extremes, its lattice section has at least $`\lfloor n\mu_i\rfloor`$ sites. In particular its projection set is the full integer interval once $`n\mu_i\ge1`$.

Choose $`n`$ large enough that all these lengths are at least $`\max(2m,\kappa_i)`$. Any geometrically active line can have at most $`m-1`$ sites changed or canceled by the other lines. It therefore has more than $`m-1`$ nonzero sites remaining on it. In a competing pattern that omits the same relative $`i`$-line, the other $`m-1`$ directions can cover at most $`m-1`$ sites of that line. This is impossible. Applying the argument in both directions shows that the full set of active relative lines is intrinsic to the pattern.

When at least two directions are active, two nonparallel relative lines recover the anchor uniquely. When one direction is active, its height and tangent phase recover the pure pattern; different phases are distinct because the section is at least one full minimal period. All such phases are realizable by translating along that direction and avoiding the finitely many other visible strips. The zero pattern is distinct and realizable. Consequently all the representations counted in (32) are distinct. This proves the complete collision classification for this class, including multiple line intersections.

For an anchor with exactly $`j`$ active directions, the coefficient with which it is counted in the first sum of (32) is
```math
\sum_{k=2}^j(-1)^k(k-1)\binom jk
 =\begin{cases}0,&j=0,1,\\1,&j\ge2.\end{cases}
\qquad\text{(34)}
```
The binomial identities for $`(1-x)^j`$ and its derivative at $`x=1`$ prove (34). Thus the first sum counts exactly the anchors with at least two active lines, establishing (32).

### 5.2 A lattice-polygon counting lemma

**Lemma 7 (integer support lines).** Let $`L_n`$ be convex polygons with a fixed finite collection of primitive integral outward normals, integer support constants, and the same facet order for all sufficiently large $`n`$. Suppose every facet length is a positive multiple of $`n`$, up to a bounded error, and all angles are fixed. Write $`\mathcal B(L_n)`$ for the sum of their facet lengths measured in the corresponding primitive lattice steps. Then
```math
|L_n\cap\mathbb Z^2|=
 \mathrm{area}(L_n)+\tfrac12\mathcal B(L_n)+O(1).
\qquad\text{(35)}
```
The bound is uniform over finitely many fixed support perturbations.

**Proof.** A line with primitive integral normal and integer support constant contains a lattice point; its lattice points are spaced by the primitive tangent vector. On each sufficiently long facet choose lattice points within one primitive step of both endpoints. The convex hull of these points agrees with the facet except within bounded neighborhoods of its endpoints. Adjacent normals are fixed and nonparallel, so the corner triangles cut off by joining these selected points have bounded diameter and area.

Let $`H_n=\mathrm{conv}(L_n\cap\mathbb Z^2)`$. It contains that selected-point hull and lies in $`L_n`$. Hence $`L_n\setminus H_n`$ lies in finitely many uniformly bounded corner neighborhoods. Its area is $`O(1)`$; the number of boundary lattice points of $`H_n`$ differs from $`\mathcal B(L_n)`$ by $`O(1)`$, since along each long facet only bounded endpoint lengths change and all new boundary pieces lie in those corner neighborhoods. Pick's theorem for the lattice polygon $`H_n`$ now gives (35), using $`H_n\cap\mathbb Z^2=L_n\cap\mathbb Z^2`$. ∎

This lemma requires **integer support constants**, not integral vertices. It is the reason determinants and vertex denominators cannot silently introduce an additional linear residue term in the calculation below.

### 5.3 Computing the boundary and area terms

For a fixed residue $`r`$, (33) has outward-normal support constants
```math
-n\ell_i+\delta_i^+(r),\qquad
 nu_i+\delta_i^-(r),
\qquad\text{(36)}
```
for normals $`\pi_i`$ and $`-\pi_i`$, respectively. These constants are integers. All facets of $`K_I`$ have positive length by (27), so bounded displacements preserve their order and change their endpoints by $`O(1)`$ for large $`n`$.

Moving an outward support line by $`\delta`$ in its primitive-normal coordinate moves it a Euclidean distance $`\delta/\|v_i\|`$. Its signed area strip is therefore $`n\delta\tau+O(1)`$. The corrections at vertices have bounded area. It follows that
```math
\begin{split}
\mathrm{area}(K_I(n))
 &=n^2\mathrm{area}(K_I)
   +n\sum_{i\in I}
       \bigl(\delta_i^+\tau_i^+(I)+\delta_i^-\tau_i^-(I)\bigr)+O(1),\\ 
\mathcal B(K_I(n))
 &=n\sum_{i\in I}(\tau_i^+(I)+\tau_i^-(I))+O(1).
\end{split}
\qquad\text{(37)}
```
Lemma 7 applies to (33). Insert (37) into (35), then into (32). Its linear term is exactly (31), including the $`\sum_i\kappa_iW_i`$ pure-direction contribution.

The quadratic coefficient is $`\sum_{|I|\ge2}c_I\mathrm{area}(K_I)`$. Apply (34) pointwise to the strips of $`-P`$. This coefficient is the area where at least two such strips overlap, namely $`\mathrm{area}(\mathcal M(-P))=\mathrm{area}(\mathcal M(P))`$. The finitely many residues make all error bounds uniform. This completes the proof of Theorem 6. ∎

## 6. Full-support Case B lines and arbitrary polygons

In this section $`m\ge2`$ is arbitrary, $`P`$ is an arbitrary positive-area rational convex polygon, and the components are as in Section 5. Replace (27) by the assumption that the **lattice** lines
```math
L_i=\{z\in\mathbb Z^2:\pi_i(z)=a_i\}
\quad\text{are pairwise disjoint}.
\qquad\text{(38)}
```
This gives CN Case B: for every nonzero color, its indicator is a sum of the corresponding single-component color indicators, while the zero indicator is a constant minus the line-support indicators. Each summand is annihilated by its tangential difference, so their product annihilates every color indicator in characteristic zero. We do not identify this geometric assumption with the general CN definition of Case B. The underlying real lines may intersect at nonlattice points, including simultaneous intersections; they need not be concurrent. The primitivity and nonzero-word assumptions remain in force.

In (29) now allow $`\tau_i^\pm(I)=0`$ when the corresponding face of $`K_I`$ is a point. Put
$`h_i^\pm(r)=|D_{v_i}^\pm(r)|`$, with (13) computed for $`P`$. Define the coefficient before pattern collisions by
```math
\begin{split}
\mathfrak B_r=\sum_i\kappa_iW_i+
\sum_{|I|\ge2}c_I\sum_{i\in I}\bigl[
 &\tau_i^+(I)(\delta_i^+(r)+\tfrac12-h_i^-(r))\\ 
 +&\tau_i^-(I)(\delta_i^-(r)+\tfrac12-h_i^+(r))\bigr].
\end{split}
\qquad\text{(39)}
```
All ingredients except the collision correction have already been defined by finite rational geometry.

**Theorem 8 (all full-support periodic-line Case B stars).** Under (38),
```math
P_\theta(S_n)=\mathrm{area}(\mathcal M(P))n^2+
 (\mathfrak B_r-\Gamma_r)n+O(1),\qquad r=n\bmod M,
\qquad\text{(40)}
```
where $`\Gamma_r\ge0`$ is the finite geometric/arithmetic formula (46) below. All collision families contributing to $`\Gamma_r`$ are classified, and their finite-window counts are computable by interval endpoints and congruences. Remaining collisions involve a uniformly bounded number of representations.

For integral $`P`$ one may take $`M=1`$; hence the linear coefficient is constant. In particular pair determinants and color periods, however large, do not force a residue-dependent linear term for integral polygons in this class.

### 6.1 A complete reduction to singleton reassignments

Call a section in direction $`i`$ long if it has at least
$`L_i^*=\max(m,\kappa_i)`$ sites, and short otherwise. Because of (38), no component ever changes or cancels a site of another component. Thus a long relative $`i`$-line has at least $`m`$ nonzero sites. If a competing representation omitted that line, its other $`m-1`$ directions could cover at most $`m-1`$ of those sites. Every long relative line is therefore intrinsic to the pattern.

An anchor with two long directions is globally unique. Anchors with at least two short active directions form a uniformly bounded list: each short direction has only boundedly many heights by Lemma 2, and a pair of transverse heights fixes at most one anchor. This includes every multi-line anchor with no long direction.

Outside that bounded list, a nonunique multi-line pattern must have exactly one long direction $`i`$ and one other active direction $`j`$. If the $`j`$-section has at least two sites, these are off the relative $`i`$-line by (38). They determine the $`j`$-line. Any competing representation with only one other direction must use that same line and hence the same anchor. A competitor with two or more short directions belongs to the bounded list already removed.

This leaves exactly the following possible linear families: one long relative $`i`$-line and a **single off-line site** $`x`$, which can be assigned to one of several other directions. The long line and the site $`x`$ are intrinsic; a choice of its assigned direction fixes the anchor. Consequently a pattern in such a family has at most $`m-1`$ representations. The same multiplicity bound applies to any pattern with an intrinsic long line and an off-line site, including partners of the bounded exceptional list.

For completeness, a pattern with no intrinsic long line has bounded support if all its representations have only short sections. Any multi-line representation of that bounded support has two bounded-length sections and is on a uniformly bounded height-pair list. Pure representations of bounded support are likewise bounded in number by Lemma 2. Pure patterns with long sections have distinct phases and cannot equal a multi-line pattern with an off-line site. It follows that removing the bounded anchor list and its collision partners changes the total representation excess by only $`O(1)`$.

Thus singleton reassignments are the **entire** source of a linear collision term, not merely an observed subclass of collisions. The following finite construction counts them without assuming that their types stabilize.

### 6.2 All singleton positions from vertex cones

For a vertex $`q`$ of $`P`$ and a residue $`r`$, let
$`\mathcal C_{q,r}=rq+T_q`$. If $`q`$ is the unique minimum vertex of $`\pi_j`$, inspect heights
$`t=\lceil r\ell_j\rceil+h`$, $`0\le h\le C_{v_j}(1)`$; if it is the unique maximum vertex, inspect $`t=\lfloor ru_j\rfloor-h`$ in that range. At each height use (14) to find the integer points of $`\mathcal C_{q,r}`$. Whenever there is exactly one such point, include it in the finite set $`Z_{q,j}(r)`$. If $`q`$ is not a unique extreme vertex for $`\pi_j`$, set $`Z_{q,j}(r)=\varnothing`$.

These sets are finite: at a strict extreme vertex every inspected real cone section is bounded, and a singleton section contributes at most one lattice point. Put
```math
Z_q(r)=\bigcup_j Z_{q,j}(r),\qquad
J_{q,z}(r)=\{j:z\in Z_{q,j}(r)\}.
\qquad\text{(41)}
```
For all sufficiently large $`n=Mk+r`$, the complete list of singleton sections of $`S_n`$ is
```math
x=kMq+z,\qquad z\in Z_q(r),
\quad\text{singleton directions exactly }J_{q,z}(r).
\qquad\text{(42)}
```
Lists at distinct vertices are then disjoint.

To prove this exact statement, Lemma 2 first confines a singleton height to the inspected endpoint bands. When an extreme face is an edge, its nearby sections have length growing linearly with $`n`$, so no such singleton persists. At a strict extreme vertex, the inspected tangent-cone section is bounded. Every one of its finitely many lattice points lies in the translated polygon for all sufficiently large $`k`$, by (18); the reverse containment holds for all $`k`$. Thus the whole lattice section, and not just its feasibility, stabilizes to the cone section. The finite lists give a common effective threshold. Their bounded offsets and the positive distances between distinct vertices prove disjointness of the vertex lists for large $`n`$.

### 6.3 The finite interval-and-phase formula

Fix a long direction $`i`$. Write $`d_{ij}=\det(v_i,v_j)`$, and choose the explicit phase modulus
```math
T_i=\mathop{\mathrm{lcm}}_{j\ne i}
 \bigl(|d_{ij}|\mathrm{lcm}(\kappa_i,\kappa_j)\bigr).
\qquad\text{(43)}
```
For an integer $`t`$, the unique real point on the fixed $`j`$-line with $`\pi_i=t`$ is
```math
W_{ij}(t)=\frac{t v_j-a_jv_i}{d_{ij}}.
```
It is a permissible lattice point exactly when both coordinates are integral. Whenever it is integral, assign the label
```math
\lambda_{ij}(t)=
 \bigl(s_i(W_{ij}(t))\bmod\kappa_i,
        w_j(s_j(W_{ij}(t)))\bigr).
\qquad\text{(44)}
```
Integrality and this label are periodic in $`t`$ with period $`T_i`$: increasing $`t`$ by $`T_i`$ adds $`(T_i/d_{ij})v_j`$, and $`T_i/d_{ij}`$ is a multiple of both relevant word periods.

For a vertex $`q`$, define the rational interval
$`I_i(q)=[\pi_i(q)-u_i,\pi_i(q)-\ell_i]`$, and the following finite unions of intervals:
```math
A_{ij}(q)=\left\{\xi\in I_i(q):
 \pi_k(q)-\frac{d_{kj}}{d_{ij}}\xi\notin[\ell_k,u_k]
 \text{ for every }k\notin\{i,j\}\right\}.
\qquad\text{(45)}
```
Endpoint inclusion does not affect any integral below. All coefficients $`d_{kj}/d_{ij}`$ in these exclusions are nonzero.

For a set of directions $`J`$, a residue $`s\bmod T_i`$, and a real number $`\xi\in I_i(q)`$, let $`C_{i,J,s}(\xi;q)`$ consist of the directions $`j\in J\setminus\{i\}`$ for which $`W_{ij}(s)`$ is integral and $`\xi\in A_{ij}(q)`$. Define
```math
\begin{split}
e_{i,J,s}(\xi;q)&=|C_{i,J,s}(\xi;q)|
 -|\{\lambda_{ij}(s):j\in C_{i,J,s}(\xi;q)\}|,\\ 
\gamma_{i,J}(q)&=\frac1{T_i}\sum_{s=0}^{T_i-1}
       \int_{I_i(q)}e_{i,J,s}(\xi;q)\,d\xi,\\ 
\Gamma_r&=\sum_{q\text{ vertex}}\ \sum_{z\in Z_q(r)}\ \sum_{i=1}^m
                  \gamma_{i,J_{q,z}(r)}(q).
\end{split}
\qquad\text{(46)}
```
This is a finite rational formula. Partition $`I_i(q)`$ at the finitely many endpoints of (45); the integrand is a constant nonnegative integer on each resulting open interval, for each of the finitely many phases. In particular $`\Gamma_r\ge0`$. A type with fewer than two candidate directions contributes zero.

### 6.4 Why (46) counts all singleton collisions

Fix an actual singleton site $`x=kMq+z=nq+O(1)`$ from (42). A representation assigning this site to direction $`j\ne i`$ has
```math
u=W_{ij}(t)-x,
\quad\text{relative }i\text{-height }a_i+\pi_i(x)-t.
\qquad\text{(47)}
```
It has exactly the required active directions precisely when its $`i`$-section is long, its $`j`$-section is the singleton $`x`$, and
```math
a_k+\pi_k(x)-\pi_k(W_{ij}(t))\notin Q_{v_k}(n)
\quad\text{for every }k\notin\{i,j\}.
\qquad\text{(48)}
```
These are exact finite-window conditions. They also show that the possible $`t`$'s lie in an interval of length $`O(n)`$.

For a fixed site $`x`$, two such entries give the same pattern exactly when their integer $`t`$'s agree and their labels (44) agree. Indeed, the intrinsic long line determines $`t`$ through (47). The baseline $`i`$-word phase is
$`s_i(W_{ij}(t))-s_i(x)`$; subtracting the common $`s_i(x)`$ preserves label equivalence. The singleton color is the second label component. Minimality of the period and the long section length make equality of baseline phases necessary and sufficient. There is no cancellation or additional point on the long line because of (38).

For leading counts, put $`\xi=t/n`$. The long-section condition differs from $`\xi\in I_i(q)`$ at only $`O(1)`$ integer values of $`t`$, by Lemma 2. Since
$`\pi_k(W_{ij}(t))=(d_{kj}/d_{ij})t+O(1)`$, conditions (48) differ from $`\xi\in A_{ij}(q)`$ at only $`O(1)`$ values as well: the projection holes and endpoint roundings are bounded, and all these slopes are nonzero. The bound is uniform over the finite site types and phases.

On an interval $`(a,b)`$ in the partition of (45), the number of integers $`t\equiv s\pmod{T_i}`$ with $`t/n\in(a,b)`$ is $`n(b-a)/T_i+O(1)`$. Each label class contributes its size minus one to the representation excess. Summing these arithmetic counts proves that the singleton-family excess at this site is
$`n\gamma_{i,J_{q,z}(r)}(q)+O(1)`$. The long direction and off-line site are intrinsic, so different $`(i,x)`$ families cannot share a pattern. Summing (42) proves a total excess $`n\Gamma_r+O(1)`$, with the bounded exceptions of Section 6.1 already accounted for.

This also gives **exact** family counts without an asymptotic approximation. Use (47)–(48), the exact projection lists (15), and the finitely many short $`i`$-heights. Their inverse images in the parameter $`t`$ are finite unions of intervals and individual points. On each phase class (44) is constant. Partition at all endpoints, group equal labels, and count integers in each arithmetic progression by floor and ceiling. The number of tests and interval types is bounded in terms of the fixed data. Thus the finiteness and enumerability of all linear collision families have been proved.

### 6.5 Counting raw boundary anchors and finishing the proof

It remains to justify (39) without (27). A constraint of $`K_I`$ may now touch only a vertex. For a fixed residue, each possible intersection of two support lines of $`K_I(n)`$ is affine in $`n`$; its feasibility eventually stabilizes because each remaining inequality is linear in $`n`$. Consequently every eventual vertex is $`n`$ times a vertex of $`K_I`$, up to $`O(1)`$. Positive-length limiting facets have length $`n\tau+O(1)`$; newly appearing facets at a limiting vertex have length $`O(1)`$.

Lemma 7 extends to this case: group all bounded facets at the same limiting vertex into one bounded corner neighborhood, and use lattice points near the endpoints of the remaining long facets. The same selected-point hull and Pick argument give (35). The area displacement formula (37) also holds with $`\tau=0`$ at point faces, since all changes there have bounded area. It follows that the full-band anchors in $`K_I(n)`$ have the coefficient obtained from (31) before projection holes are removed.

By Lemma 3 a lower projection hole of direction $`i`$ removes an integer anchor line at bounded distance from the upper limiting facet $`\pi_i=-\ell_i`$. Its intersection with $`K_I(n)`$ has lattice count
$`n\tau_i^+(I)+O(1)`$, since the line has primitive integral normal and integer intercept. This remains valid with $`\tau=0`$. Upper holes similarly remove $`n\tau_i^-(I)+O(1)`$ points. Parallel hole lines are disjoint; hole lines of different directions meet in at most one point, so all overlaps alter only $`O(1)`$ anchors. Therefore the number of anchors with every direction in $`I`$ actually active is
```math
|K_I(n)\cap\mathbb Z^2|
 -n\sum_{i\in I}\bigl(h_i^-(r)\tau_i^+(I)
                         +h_i^+(r)\tau_i^-(I)\bigr)+O(1).
\qquad\text{(49)}
```
Apply the inclusion-exclusion identity (34) to (49). Add $`\kappa_i|Q_i(n)|`$ pure representations for every direction and one zero representation. Short pure sections may identify phases or singleton directions, but by Lemma 2 this changes the count by only $`O(1)`$. The result is
```math
\#\text{raw representations}
 =\mathrm{area}(\mathcal M(P))n^2+\mathfrak B_r n+O(1).
\qquad\text{(50)}
```
Section 6.1 proves completeness of the collision classification, and Sections 6.2–6.4 give the excess $`n\Gamma_r+O(1)`$. Subtract it from (50) to obtain (40). All constructions use only finitely many residues of $`M`$; color periods have been averaged over the finite moduli (43). With integral $`P`$, every geometric type is independent of $`n`$, so $`\mathfrak B_r-\Gamma_r`$ is constant. This completes the proof of Theorem 8. ∎

The theorem includes, and analytically recovers, the earlier three-line triangle coefficient: directions $`(1,0),(1,2),(1,4)`$, levels $`0,1,2`$, constant binary words give $`\mathfrak B=3`$ and $`\Gamma=1/4`$, hence $`\beta=11/4`$. If the middle word is changed to $`(1,2)`$ over $`\mathbb F_3`$, the allowed phase labels eliminate that linear collision family: $`\mathfrak B=6`$, $`\Gamma=0`$, and $`\beta=6`$. These values follow from (39) and (46), not from polynomial fitting.

### 6.6 Removing Case B for full-support lines

**Theorem 9 (arbitrary full-support periodically colored line stars).** Formula (40), with the same explicitly defined $`\mathfrak B_r`$ and $`\Gamma_r`$, holds without (38). Thus for arbitrary $`m\ge2`$, arbitrary primitive pairwise nonparallel directions and integer line levels, and arbitrary nonzero bilateral periodic words, every rational polygon has
```math
P_\theta(S_n)=\mathrm{area}(\mathcal M(P))n^2+
 \beta_{n\bmod M}(\theta,P)n+O(1),
\qquad \beta_r=\mathfrak B_r-\Gamma_r.
\qquad\text{(51)}
```
Integral polygons have a constant $`\beta`$. All collision families of linear size are precisely those in (46); all additional collisions have a uniformly bounded representation list and computable partners.

**Proof.** Increase the long thresholds to
$`\widehat L_i=\max(2m,\kappa_i)`$. At most $`m-1`$ sites of a component line can be changed or canceled by the other directions. A long relative line therefore retains more than $`m-1`$ nonzero sites. It is intrinsic by the same covering argument as in Section 5.1, so two long directions still recover an anchor globally.

Let $`X`$ be the finite set of integral intersections of pairs of the fixed affine lines. Let $`V_j(n)`$ be the set of sites of $`S_n`$ belonging to short sections in direction $`j`$, defined before adding colors. Lemma 2 gives the explicit bound
```math
|V_j(n)|\le \widehat L_j R_j(\widehat L_j),
\qquad R_j(K)=2(C_{v_j}(K)+1).
\qquad\text{(52)}
```
Consider anchors at which a short $`j`$-section contains an intersection with another active component. For some fixed $`b\in X\cap L_j`$, its relative intersection site is $`b-u\in V_j(n)`$; hence $`u=b-x`$ for one of the boundedly many $`x\in V_j(n)`$. Their number is at most
```math
\sum_j |X\cap L_j|\,\widehat L_j R_j(\widehat L_j),
\qquad\text{(53)}
```
uniformly in $`n`$. Also discard anchors with at least two short active directions; their number is bounded by
$`\sum_{i\lt j}R_i(\widehat L_i)R_j(\widehat L_j)`$.

Outside these bounded lists, every non-bulk multi-line anchor has one long direction and one short direction, and **their visible supports are disjoint**. All color interactions involving a short component have been removed. Thus a short section of size at least two fixes its line from two off-line sites, exactly as in Section 6.1; only a singleton can be reassigned. The phase labels and the exact conditions (47)–(48) apply unchanged to these singleton families. Exclude $`t=a_i`$, when the singleton lies on the relative long line; this removes only finitely many candidates per site type and does not change (46).

We must also bound collision partners of the discarded anchors; merely bounding the anchors themselves would not suffice. If such a pattern has an intrinsic long line and a nonzero site off it, a choice of one of the other $`m-1`$ directions through that site fixes its anchor, so it has at most $`m-1`$ multi-line representations. If it has an intrinsic long line and no off-line site, its pure representations have that one height and at most $`\kappa_i`$ phases. Any multi-line partner outside the discarded lists would have either two intrinsic long lines or a nonzero off-line short component, both impossible for this pattern. Its other partners are on the discarded lists.

Finally, if a discarded pattern has no intrinsic long line, all sections of that representation are short, so its support has bounded size, at most $`\sum_i\widehat L_i`$. In any competing representation a component with section length $`L`$ retains at least $`L-(m-1)`$ nonzero sites. Thus every competing section has bounded length as well. A multi-line competitor then has two bounded-length sections and belongs to a uniformly bounded height-pair list by Lemma 2; pure competitors also belong to bounded short-height lists. All collision partners are therefore uniformly bounded in number. They are enumerable from the lists just described, or from a fixed off-line site and a choice of its assigned direction.

After removing these bounded representations, the complete proof of Section 6.1 applies. Formula (46) still gives the entire linear excess; the change of long threshold removes only boundedly many heights in each singleton family. Raw geometric anchor counts, projection holes, and pure-direction leading counts do not depend on intersection colors, so (39) and (50) are unchanged. Subtracting the excess proves (51). All arguments allow multiple integral intersections and arbitrary cancellation there. ∎

The removal of Case B here is specific to **full nonzero line support**: short sections contain only boundedly many sites, while all color interaction points of distinct lines are fixed and finite. It is not an argument for general finite-width or zero-containing components.

**Exact finite-window counting procedure.** The proof gives more than a formula for the linear coefficient. Define the exact raw count
```math
R_n=\sum_{|I|\ge2}c_I\,
 \#\{u\in\mathbb Z^2:a_i-\pi_i(u)\in Q_{v_i}(n)\ (i\in I)\}
 +\sum_i\kappa_i|Q_{v_i}(n)|+1.
\qquad\text{(54)}
```
This counts all multi-line anchors, all pure phase entries and zero with multiplicity. Each intersection count in (54) is a bounded collection of rational polygons, deleted integer lines and their intersection points, by (15). It can be computed exactly, for example by integer vertical sections of each polygon and floor/ceiling counts on each deleted line. Pure entries on all long rows are distinct; every possible pure-phase or cross-direction identification is on a bounded short-row list.

Take the bounded bad-anchor lists used in Theorem 9, add the zero entry and all pure entries on rows with length at most $`H=\sum_i\widehat L_i+m`$, and close this list under pattern equality. Call the resulting representation list $`C_n`$. The collision-partner argument in that proof shows that its size is uniformly bounded and gives a finite way to enumerate it: use a long line and one off-line site to test each assigned direction, or use the enlarged bounded-height lists for a bounded-support pattern. If a pattern has a long line and no off-line site, test the finitely many pure phases and the bad-anchor lists.

Outside $`C_n`$, every repeated pattern is exactly a singleton label class of (47)–(48), with $`t\ne a_i`$. Use the exact interval/congruence procedure of Section 6.4, excluding the finitely many label classes already in $`C_n`$, and let $`F_n^{\mathrm{out}}`$ be their total excess. Then
```math
P_\theta(S_n)=R_n-
 \bigl(|C_n|-|\{\text{patterns in }C_n\}|\bigr)-F_n^{\mathrm{out}}.
\qquad\text{(55)}
```
This is an exact algorithm for every sufficiently large $`n`$, with a proved finite classification. It is not an assertion that the bounded remainder in (51) is constant or eventually quasipolynomial. It makes explicit how boundary-shell anchors, pure-direction modes and cross-representation collisions can all be counted without an unspecified collision search.

## 7. Exact checks and interpretation

The standard-library script [full_support_lines_audit.py](../scripts/full_support_lines_audit.py) computes all coefficients by cone certificates and formulas (31), (39), (46), **before** enumerating patterns. Its output is `full_support_lines_audit.json`. It does not estimate coefficients by fitting counts. The enumeration covers all lattice translations by listing geometric multi-line anchors, globally realizable pure phases, and zero; it does not sample a finite translation box and call it exhaustive.

For Theorem 4 the script also verifies that every long two-line entry is globally unique and that the entire representation excess equals the excess of the exceptional sublist $`E_n`$. For Theorem 6 it verifies that every multi-line anchor and every pure entry is distinct, including test cases with cancellation and nonunimodular directions. For Theorems 8–9 it compares the global multiplicity list with an independent enumeration of (47)–(48), checks the singleton pattern labels, and bounds the excess left on the exceptional lists. Rational cone inequalities, thresholds, witnesses, interval partitions and phase labels are stored in the JSON. Larger exact family counts avoid constructing pattern bitsets; these are consistency checks of the proved density, not fits.

Some coefficients derived by these formulas are:

| Configuration and window | $`\Lambda`$ | $`\beta_r`$ |
| --- | ---: | --- |
| Two lines, triangle, directions $`(1,0),(1,2)`$, levels $`0,1`$, words $`(1),(1,2)`$ over $`\mathbb F_3`$ | $`3/2`$ | $`17/2`$ |
| Two disjoint lines, triangle, directions $`(1,2),(1,4)`$, levels $`1,2`$, constant binary words | $`15/2`$ | $`1/2`$ |
| Two axis lines, triangle, words $`(1,4),(4,2,3)`$ over $`\mathbb F_5`$, with intersection cancellation | $`1`$ | $`7`$ |
| Two lines, rational trapezoid with vertices $`(0,0),(5/3,0),(1,4/3),(0,1)`$, directions $`(1,2),(2,1)`$, levels $`0,0`$, words $`(1,2),(3)`$ over $`\mathbb F_5`$ | $`143/27`$ | $`124/9,335/27,337/27`$ for $`r=0,1,2`$ |
| Three lines, integral zonotope $`\sum_{v\in\{(1,0),(0,1),(1,1)\}}[-v,v]`$, levels $`0,0,0`$, words $`(1),(1,4),(2,3,4)`$ over $`\mathbb F_5`$ | $`24`$ | $`36`$ |
| Three lines, one third of the preceding zonotope, levels $`0,1,2`$, words $`(1),(1,2),(3)`$ over $`\mathbb F_5`$ | $`8/3`$ | $`28/3,4,20/3`$ for $`r=0,1,2`$ |
| Four lines, integral zonotope generated by $`(1,0),(1,2),(0,1),(2,1)`$, levels $`0,1,-1,2`$, words $`(1),(1,2),(2),(3,4)`$ over $`\mathbb F_5`$ | $`168`$ | $`100`$ |

The examples establish computational consistency of the proofs. Constant remainders seen at sampled scales are not promoted to additional exact-polynomial theorems.

## 8. Relation to the all-star theorem and scope

The [eventual-quasipolynomial theorem](eventual-quasipolynomial.md) covers every star on rational polygons, including zeros, wider strips, unequal tails and periodic backgrounds. Its proof uses the finite-description and representative argument followed by Bogart–Goodrick–Woods, rather than the collision classification on this page.

For full-support lines, combine that theorem with (51): on each common refined residue a nonzero difference of linear coefficients would contradict the bounded remainder in (51). Hence the remainder after the coefficient (39)–(46) is eventually periodic. The vertex-denominator modulus used here is a sufficient period of the **linear coefficient**, not a claimed sufficient period of the entire remainder.

This page proves a direct coefficient formula only for nonzero periodic words on whole affine lattice lines with zero background. It does not extend its singleton classification to arbitrary star components. A [separate supplementary construction](explicit-coefficients.md), established in one fresh-context GPT-6.1 Sol adversarial review on 25 October 2026, addresses general stars through eliminated representative polygons. It is not a premise of any theorem on this page.

## 9. Attribution and review boundary

The star hypotheses are CN §0.1, printed p. 4. The area coefficient used for comparison is proved in [the leading-term geometry](leading-term.md) and [the all-star extension](all-star-extension.md). The cone certificates, two-line collision list, facet calculation, singleton reassignment classification and bounded intersection reduction are proved in Sections 3–6 above. Pick's theorem is the classical planar lattice-polygon theorem used in Lemma 7.

The earlier example-specific coefficient is also proved in [the three-line module](three-line-second-order.md). The present page independently recovers it from (39) and (46). Finite checks support the proofs; no experimentally constant remainder or fitted period is promoted to a theorem.

The integration preserves the mathematical scope of the separately reviewed manuscript. Its source fingerprint and the review fingerprint are in [the review-status record](../sources/review-status.md). No claim of human expert certification or comprehensive literature priority is made.

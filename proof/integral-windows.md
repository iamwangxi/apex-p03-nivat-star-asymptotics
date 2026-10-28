# Nivat stars on integral convex polygons have constant linear coefficients

**Status: established in one fresh-context AI adversarial review.** The review found no fatal, proof-repair or presentation issue. Source and review fingerprints are in the [review-status record](../sources/review-status.md). **The linear coefficient is constant for every CN star on every integral convex polygon.** The full proof and finite boundary formula are in Appendix B. Earlier sections give independently checkable special cases and an exact period-four constant-term example. Computations are finite checks, not premises of the proofs. No human certification, proof-assistant verification or literature-priority certification is claimed.

## 1. Answer and scope

For the full class in CN §0.1, Theorem 5 proves

```math
f(n)=\Lambda(P)n^2+\beta n+\gamma_{n\bmod M}
\quad(n\text{ sufficiently large}),\qquad p_\beta=1.
```

There is no full-support, equal-tail, unit-period, or noncancellation assumption. The existence of the periodic remainder uses the [all-star quasipolynomial theorem](eventual-quasipolynomial.md); constancy of the linear coefficient follows from the bounded-error formula proved here. The finite search found no counterexample, but absence in the search is not used in this conclusion.

The first subclass established in this note is as follows. Let

```math
\theta(z)=B(z)+\sum_{i=1}^m F_i(z),\qquad
F_i(z)=\mathbf1[\pi_i(z)=a_i]w_i(s_i(z)),\qquad z\in\mathbb Z^2,
```

with arithmetic in a prime field. The vectors $`v_i`$ are primitive and pairwise nonparallel, $`\pi_i(z)=\det(v_i,z)`$, and $`(v_i,e_i)`$ is an integral determinant-one basis, so $`s_i(z)=\det(z,e_i)`$. Each $`a_i`$ is integral. Each periodic word $`w_i`$ is not identically zero, and its minimal period is $`\kappa_i`$. **Zero entries are allowed in the words.** The field $`B`$ is any doubly periodic field. Affine line intersections, multiple intersections, and cancellations are allowed.

**Theorem 1.** For every positive-area integral convex polygon $`P`$, there is a rational constant $`\beta`$ such that

```math
P_\theta(nP\cap\mathbb Z^2)=\Lambda(P)n^2+\beta n+O(1).
```

Consequently, in the established eventual quasipolynomial, the minimal period of the linear coefficient is one. A finite rational formula for $`\beta`$ is given in §5. The constant remainder is eventually periodic by the [all-star quasipolynomial theorem](eventual-quasipolynomial.md); its period need not be one.

These are CN stars: each nonzero line-supported component has a tangential period and is not doubly periodic, with zero tails. Absorb $`B`$ into one component if the definition requires a sum with no separately written background. This preserves the required non-double-periodicity.

Theorem 1 includes sparse periodically colored affine lines and a global periodic background on every integral polygon. Appendix A extends this to all equal-tail stars. Appendix B then removes the equal-tail restriction by preserving complete one-direction end-field labels, rather than comparing gradients alone. These narrower results are retained as independently checkable stages of the general argument.

**Theorem 2, an exact constant-term example.** Over $`\mathbb F_2`$, put

```math
\theta(x,y)=\mathbf1[y=x,\ x\equiv0\pmod2]
             +\mathbf1[y=-x,\ x\equiv0\pmod2],\qquad
P=[0,1]^2.
```

For every $`n\ge12`$,

```math
\boxed{\quad P_\theta(\{0,\ldots,n\}^2)
       =2n^2+8n+\gamma_{n\bmod4},\qquad
       (\gamma_0,\gamma_1,\gamma_2,\gamma_3)=(2,3,2,-1).\quad}
```

Thus $`p_\beta=1`$ and $`p_\gamma=p_f=4`$. The proof in §7 is exact and does not infer a formula by fitting counts. The two primitive directions have determinant two and both words have period two. Their least common multiple, two, would miss this constant-term period.

## 2. Finite endpoint profiles

Write

```math
\ell_i=\min_P\pi_i,\qquad u_i=\max_P\pi_i,\qquad W_i=u_i-\ell_i>0.
```

These extrema are integers because the vertices of $`P`$ are integral. A geometric section is

```math
Q_i(n,t)=\{z\in nP\cap\mathbb Z^2:\pi_i(z)=t\}.
```

Its points, if any, are consecutive lattice sites in direction $`v_i`$. For any fixed integer $`R`$, all sections with fewer than $`R`$ sites occur at a uniformly bounded distance in height from one of the two extreme faces. Their eventual lists are finite and computable.

Here is a justification sufficient for all uses below. Choose a strictly interior section of $`P`$ with positive real length in primitive $`v_i`$ steps. Concavity of real section length gives a positive rational $`c_i`$ such that the real section length in $`nP`$ is at least

```math
c_i\min(t-n\ell_i,nu_i-t)
```

between the extremes. A real interval of length at least $`R+1`$ contains at least $`R`$ lattice parameters. Hence a computable constant height band contains every shorter section. If an extreme face is an edge, every fixed-height section in that band has length growing with $`n`$, so none remains short eventually. If the extreme face is the unique vertex $`q`$, translate by $`-nq`$ and use the origin-based rational tangent cone $`T_q=\mathrm{cone}(P-q)`$. For every inspected offset $`d`$, the section of $`T_q`$ is bounded and finite on the lattice. Each of its finitely many points satisfies all nonincident facets of $`n(P-q)`$ once $`n`$ exceeds an explicit rational bound. Conversely $`n(P-q)\subseteq T_q`$. This proves exact eventual equality of the sections, with a computable threshold.

In particular, at a strict minimum vertex use

```math
X_{i,q,d}=\{z\in T_q\cap\mathbb Z^2:\pi_i(z)=d\},\qquad d\ge0,
```

and use $`\pi_i(z)=-d`$ at a strict maximum vertex. The actual section is $`nq+X_{i,q,d}`$ for all sufficiently large $`n`$. The translation is integral for **every** $`n`$. This fact is central to the proof.

Let $`h_i^-`$ and $`h_i^+`$ count the eventual missing geometric heights near the minimum and maximum respectively. Then

```math
\#\{t:Q_i(n,t)\ne\varnothing\}=nW_i+1-h_i^--h_i^+
```

eventually. The [projection-cone construction](full-support-lines.md#3-exact-projection-counts-from-finitely-many-cones) computes these integers. An extreme edge contributes no such holes.

For $`L\ge1`$ define the zero-run fraction

```math
z_i(L)=\frac1{\kappa_i}
 \#\{r\bmod\kappa_i:w_i(r)=\cdots=w_i(r+L-1)=0\}.
```

It vanishes for $`L\ge\kappa_i`$. Define

```math
e_i^- =\sum_{\substack{\text{minimum endpoint profiles }X\ne\varnothing}}
                  z_i(|X|),\qquad
e_i^+ =\sum_{\substack{\text{maximum endpoint profiles }X\ne\varnothing}}
                  z_i(|X|).
```

Only profiles with $`|X|<\kappa_i`$ matter, so these are finite rational sums. They count phase-averaged invisible, geometrically nonempty sections; they are separate from geometric holes.

## 3. Separating the background and recovering long lines

Choose $`Q`$ with $`Q\mathbb Z^2\subseteq\mathrm{Per}(B)`$, taking $`Q=1`$ when $`B=0`$. Equality of two sufficiently large window patterns forces their translated backgrounds to be identical on the whole lattice.

Indeed, subtracting the two representations would express the doubly periodic difference $`B(u+\cdot)-B(u'+\cdot)`$ as a function supported on at most $`2m`$ affine lines. If that periodic difference were nonzero, it would be nonzero on some full residue class modulo $`Q`$. The polygon $`nP`$ contains a square whose side grows with $`n`$. In the chosen residue class this square contains a lattice grid of growing side. A fixed number of lines cannot cover that grid: in an $`L`$ by $`L`$ grid, each line contains at most $`L`$ points, whereas there are $`L^2`$ points. Taking $`L>2m`$ gives a contradiction. This also gives a uniform effective threshold. Consequently one may subtract a common background when analyzing collisions.

Let $`b_i`$ be the least positive integer such that $`b_i v_i\in\mathrm{Per}(B)`$, and put

```math
K_i=\mathrm{lcm}(\kappa_i,b_i).
```

For any fixed long pure $`i`$-section, its distinct joint line/background phases are precisely the $`K_i`$ tangential phases. Each is globally realizable: translate along $`K_i v_i`$ far enough to avoid all other affine lines. Thus the pure-mode contribution is

```math
\sum_i K_iW_i n+O(1).
```

Call a geometric $`i`$-section long when its length is at least $`H_i=2m\kappa_i`$; otherwise call it short. A long section has at least $`2m`$ nonzero component sites. At most $`m-1`$ sites can be changed by the other lines, so more than $`m-1`$ nonzero sites remain. A competing representation omitting the same relative $`i`$-line could cover at most $`m-1`$ sites of that line. This is impossible. Therefore a long relative line is intrinsic to the pattern, even with sparse colors and intersection cancellation. Two long directions determine the anchor uniquely.

The number of anchors with two geometrically short active directions is uniformly bounded: each short direction has a bounded list of possible heights, and two nonparallel heights determine at most one anchor. This bound is uniform in $`n`$ even though the short heights themselves move.

We also remove a uniformly bounded list where a short section meets another physical supporting line. Let $`Z`$ be the finite set of integral intersections of the physical affine lines. Every site of every short section belongs to a finite vertex list $`nq+\delta`$. A meeting at $`z_*\in Z`$ forces the anchor to be $`z_*-nq-\delta`$. There are only boundedly many such anchors. This covers multiple intersections and every cancellation involving a short section.

Outside these two bounded lists, a nonunique multi-line representation has exactly one long direction $`i`$ and exactly one other nonempty component $`j`$, with disjoint visible supports. If the latter has two nonzero sites, they determine its supporting line and hence the anchor. The only remaining collisions of linear total size have an intrinsic long $`i`$-line and a **single nonzero off-line site**, reassigned between different possible directions. Sparse geometric sections may contain several sites, but their other colors must be zero. This is the only change to the singleton classification needed here.

For completeness, removing bounded lists must not hide unboundedly many partners. A pattern with an intrinsic long line and a nonzero off-line site has at most $`m-1`$ multi-line representations, because choosing the direction through that site fixes the anchor. A pattern with an intrinsic long line and no off-line site has finitely many pure phases; every multi-line partner must belong to the removed short-intersection list. A pattern with no long line has only short-section representations: a long line in a partner would be intrinsic in the original representation as well, and the geometric length of a fixed relative line is the same in both. Such multi-line representations and their pure short partners lie on bounded lists. Thus closing the removed lists under pattern equality changes the total representation excess by only $`O(1)`$.

## 4. Counting anchors, including invisible sparse sections

For $`I\subseteq\{1,\ldots,m\}`$ with $`|I|\ge2`$, set

```math
c_I=(-1)^{|I|}(|I|-1),\qquad
K_I=\bigcap_{i\in I}\{x:-u_i\le\pi_i(x)\le-\ell_i\}.
```

Let $`\tau_i^+(I)`$ be the length, in primitive $`v_i`$ steps, of the face of $`K_I`$ on $`\pi_i=-\ell_i`$, and let $`\tau_i^-(I)`$ be the corresponding length on $`\pi_i=-u_i`$. A point face has length zero. These lengths are rational and computable by intersection of intervals on the supporting line.

The full band polygon at scale $`n`$ is

```math
K_I(n)=\bigcap_{i\in I}
\{x:a_i-nu_i\le\pi_i(x)\le a_i-n\ell_i\}.
```

Its support normals are primitive and its support constants are integers. A polygon with these fixed normals satisfies

```math
|K_I(n)\cap\mathbb Z^2|
=\mathrm{area}(K_I)n^2
+n\sum_{i\in I}\left[
\tau_i^+(I)(a_i+\tfrac12)
+\tau_i^-(I)(-a_i+\tfrac12)\right]+O(1).
```

Here rational vertices create no additional linear-period error. To see this, choose lattice points within bounded distance of the ends of every long integer support facet. Their hull cuts off only bounded corner neighborhoods. Its area differs by $`O(1)`$ and its boundary lattice count differs by $`O(1)`$ from the sum of primitive facet lengths. Pick's theorem proves the half-boundary term. Facets of bounded length at a limiting vertex are absorbed into that vertex's bounded neighborhood. Moving a support line by a fixed amount $`a_i`$ adds its signed area strip $`na_i\tau+O(1)`$, giving the area-displacement term.

A missing geometric height removes $`n\tau_i^\pm(I)+O(1)`$ anchors on an integer line. An invisible sparse section has the same kind of correction, with a phase fraction. More explicitly, at a strict endpoint vertex $`q`$, fix a nonempty profile $`X`$ and write anchors on its height as

```math
u=-nq+w_0+t v_i,\qquad t\in\mathbb Z.
```

The component values on its section are $`w_i(s_i(w_0+z)+t)`$ for $`z\in X`$. There is no $`n`$ in this expression. The other band inequalities confine $`t`$ to an interval of length $`n\tau_i^\pm(I)+O(1)`$. Exactly a fraction $`z_i(|X|)`$ of its residue classes make the entire section invisible. Counting arithmetic progressions gives the correction

```math
n\tau_i^\pm(I)z_i(|X|)+O(1).
```

Deletions from different directions overlap in only boundedly many anchors: they each prescribe a height from a bounded endpoint list, and two independent height equations fix one anchor. Therefore if $`N_I(n)`$ counts anchors where every component in $`I`$ has nonempty colored support, then

```math
\begin{split}
N_I(n)=\mathrm{area}(K_I)n^2+n\sum_{i\in I}\big[&
\tau_i^+(I)(a_i+\tfrac12-h_i^--e_i^-)\\ 
&+\tau_i^-(I)(-a_i+\tfrac12-h_i^+-e_i^+)\big]+O(1).
\end{split}
```

For an anchor with $`j`$ nonempty components,

```math
\sum_{k=2}^j(-1)^k(k-1)\binom jk=\mathbf1[j\ge2].
```

Thus $`\sum_Ic_IN_I(n)`$ counts exactly the anchors with at least two nonempty component supports. Adding the pure modes and the bounded number of backgrounds gives raw representations with quadratic coefficient $`\Lambda(P)`$ and linear coefficient

```math
\boxed{\begin{split}
\mathfrak B={}&\sum_iK_iW_i+\sum_{|I|\ge2}c_I\sum_{i\in I}\big[
\tau_i^+(I)(a_i+\tfrac12-h_i^--e_i^-)\\ 
&\hspace{42mm}+\tau_i^-(I)(-a_i+\tfrac12-h_i^+-e_i^+)\big].
\end{split}}
```

The area identity follows from the same binomial identity applied pointwise to the direction bands; reflection preserves area. All quantities in this formula are independent of $`n`$.

## 5. A finite formula for the entire linear collision correction

This section specifies the sparse analogue of the [full-support singleton integral](full-support-lines.md), including its background-phase condition.

For each vertex $`q`$, let $`Z_q`$ be the union over directions $`j`$ of all points in endpoint cone sections $`X_{j,q,d}`$ of length less than $`2\kappa_j`$. Retain only strict endpoint vertices for that direction. This is a finite set. Every possible single nonzero $`j`$-section is located at $`x=nq+\delta`$ with $`\delta\in Z_q`$: an interval of length at least $`2\kappa_j`$ contains at least two nonzero entries in every phase.

Fix a long direction $`i`$. Write $`d_{ij}=\det(v_i,v_j)`$ and set

```math
T_i=\mathrm{lcm}_{j\ne i}
       \left(|d_{ij}|\mathrm{lcm}(\kappa_i,\kappa_j,Q)\right),
\qquad
W_{ij}(t)=\frac{t v_j-a_jv_i}{d_{ij}}.
```

The point $`W_{ij}(t)`$ has $`\pi_i=t`$ and $`\pi_j=a_j`$. A representation assigning the singleton $`x=nq+\delta`$ to direction $`j`$ has anchor

```math
u=W_{ij}(t)-x.
```

It is admissible only if $`W_{ij}(t)`$ is integral and the endpoint profile of direction $`j`$ through $`\delta`$ has exactly this one nonzero site. In explicit finite form, put $`z=W_{ij}(t)`$ and let

```math
X=\{y\in T_q\cap\mathbb Z^2:\pi_j(y)=\pi_j(\delta)\}.
```

Require $`X`$ to be one of the catalogued short profiles, $`\delta\in X`$, and

```math
w_j(s_j(z))\ne0,\qquad
w_j(s_j(z-\delta+y))=0\quad(y\in X\setminus\{\delta\}).
```

Assign the label

```math
\lambda_{ij}(t)=
\big(s_i(W_{ij}(t))\bmod\kappa_i,
      w_j(s_j(W_{ij}(t))),
      [T^{W_{ij}(t)}B]\big).
```

The last entry is the finite equivalence class of the translated periodic background, computable from its table. Integrality, admissibility and this label are periodic in $`t`$ with period $`T_i`$. This follows because $`t\mapsto t+T_i`$ adds $`(T_i/d_{ij})v_j`$, a multiple of both word periods and of $`Q`$.

Let

```math
I_i(q)=[\pi_i(q)-u_i,\pi_i(q)-\ell_i],
```

```math
A_{ij}(q)=\left\{\xi\in I_i(q):
\pi_k(q)-\frac{d_{kj}}{d_{ij}}\xi\notin[\ell_k,u_k]
\text{ for every }k\notin\{i,j\}\right\}.
```

Endpoint conventions do not affect the integrals below. For $`\delta\in Z_q`$, $`r\bmod T_i`$ and $`\xi\in I_i(q)`$, let $`C_{i,q,\delta,r}(\xi)`$ consist of directions $`j\ne i`$ passing the preceding finite integrality and singleton tests at $`t=r`$, with $`\xi\in A_{ij}(q)`$. Set

```math
E_{i,q,\delta,r}(\xi)=|C_{i,q,\delta,r}(\xi)|
 -\#\{\lambda_{ij}(r):j\in C_{i,q,\delta,r}(\xi)\},
```

```math
\boxed{\quad
\Gamma=\sum_{q}\sum_{\delta\in Z_q}\sum_i
 \frac1{T_i}\sum_{r=0}^{T_i-1}
    \int_{I_i(q)}E_{i,q,\delta,r}(\xi)\,d\xi,
\qquad \beta=\mathfrak B-\Gamma.\quad}
```

This is a finite rational formula: partition the intervals at the finitely many rational endpoints in the exclusions defining $`A_{ij}`$; the integrand is constant on each resulting open interval. All sums are finite, and $`\Gamma\ge0`$.

To verify completeness and the leading count, §3 has already reduced every unbounded collision family to a fixed intrinsic long line and one off-line site $`x`$. For this fixed $`x`$, its relative long height determines the same integer $`t`$ in every representation. Equality of long phases is equality of the first label entry, because the common subtraction $`s_i(x)`$ cancels. Equality of singleton colors is the second entry. Equality of background phases is the third entry: translating both backgrounds back by $`x`$ preserves their equality relation. Hence the displayed labels are necessary and sufficient.

The long-section condition differs from $`t/n\in I_i(q)`$ at only boundedly many integer values. Every unassigned component has geometric height

```math
a_k+\pi_k(x)-\pi_k(W_{ij}(t))
=n\pi_k(q)-\frac{d_{kj}}{d_{ij}}t+O(1).
```

The coefficient of $`t`$ is nonzero. Thus the exact inactive conditions, endpoint holes and possible invisible short sections differ from $`t/n\in A_{ij}(q)`$ at only boundedly many $`t`$. The condition that $`x`$ lie off the long line excludes just $`t=a_i`$. On a partition interval $`(a,b)`$, the number of integers $`t\equiv r\pmod{T_i}`$ with $`t/n\in(a,b)`$ is $`n(b-a)/T_i+O(1)`$. Each equal-label class contributes its size minus one to the representation excess. This proves total excess $`n\Gamma+O(1)`$. Different intrinsic long-line/site families cannot share a pattern. The bounded exceptions and all their partners were controlled in §3.

Combining this with §4 proves Theorem 1. Notice exactly where integrality is used: the local site is $`nq+\delta`$ with integral $`nq`$, so every singleton and zero-run test depends on the fixed physical point $`W_{ij}(t)`$ and fixed offsets, with no remaining $`n`$-phase. Tangential phases are then averaged over long integer intervals. Vertex denominators can still affect the bounded remainder, as §7 demonstrates.

## 6. The particularly short two-direction formula

When $`m=2`$, $`\Gamma=0`$: there is only one possible direction for an off-line singleton once its long direction is fixed. Let $`D=|\det(v_1,v_2)|`$ and $`b_i^{\mathrm{proj}}=1-h_i^--h_i^+`$. The two faces in direction $`i`$ of the anchor parallelogram have primitive length $`W_j/D`$. Therefore

```math
\boxed{\quad
\beta=K_1W_1+K_2W_2+
\frac{W_2(b_1^{\mathrm{proj}}-e_1^--e_1^+)
      +W_1(b_2^{\mathrm{proj}}-e_2^--e_2^+)}{D}.
\quad}
```

The affine levels cancel from this coefficient. They may still affect the bounded remainder. For two dotted perpendicular lines on the standard integral triangle, $`\kappa_1=\kappa_2=2`$, $`W_1=W_2=D=1`$, there are no geometric holes and each endpoint singleton is invisible for half its phases. Thus $`\beta=4+1=5`$, instead of the full-support coefficient four.

## 7. Exact proof of the period-four constant term

First consider the full-support binary diagonal pair

```math
\psi(X,Y)=\mathbf1[Y=X]+\mathbf1[Y=-X]
```

on the rectangular lattice window $`R_{a,b}=\{0,\ldots,a\}\times\{0,\ldots,b\}`$ with $`a,b\ge5`$. Put $`L=a+b+1`$.

The relative line heights range over $`[-a,b]`$ for $`Y-X`$ and $`[0,a+b]`$ for $`Y+X`$. A pair of heights comes from an integer anchor exactly when the heights have equal parity. Consequently the number of geometrically two-active anchors is

```math
N(a,b)=\begin{cases}
L^2/2,&L\text{ even},\\ 
(L^2+(-1)^a)/2,&L\text{ odd}.
\end{cases}
```

All these two-active entries, the $`2L`$ nonzero pure entries, and zero are distinct for the indicated large rectangles. Here is a complete reason, including the short entries. A section with at least three sites remains nonzero at at least two sites after its possible intersection cancellation, so its relative line is intrinsic. Two long lines fix the anchor. With one long line, a nonzero off-line site fixes the other line, while a singleton at the intersection is a nonzero edit of the fixed full word and cannot be a pure pattern. A different intersection-edit location also gives a different pattern. Sections of length at most two lie within one lattice step of their strict extreme corners. The first diagonal's corners are $`(a,0),(0,b)`$, and the second's are $`(0,0),(a,b)`$. These four neighborhoods are disjoint when $`a,b\ge5`$. If both directions are short, their nonempty sections are disjoint and their direction assignments are forced by these corner neighborhoods. A competing representation cannot have a long line: that line would have to be intrinsic in the original representation, where its geometric section would then be the same short section. Thus both short lines, and hence the anchor, are recovered as well. The same corner argument excludes pure/pure and pure/multi coincidences among short entries. It follows that

```math
P_\psi(R_{a,b})=N(a,b)+2L+1. \qquad\text{(*)}
```

Now return to the sparse pair of Theorem 2. It is supported on $`(2\mathbb Z)^2`$, and on that sublattice it is exactly $`\psi`$ after division by two. An anchor residue modulo two selects one of the four residue subgrids in $`\{0,\ldots,n\}^2`$. The associated coarse rectangle has

```math
a=\left\lfloor\frac{n-r_x}{2}\right\rfloor,\qquad
b=\left\lfloor\frac{n-r_y}{2}\right\rfloor,
\qquad (r_x,r_y)\in\{0,1\}^2.
```

Every coarse anchor is realizable. Nonzero patterns from different residue subgrids cannot agree, because their supports lie in disjoint residue classes. Their zero patterns are the one common pattern. Thus, exactly,

```math
f(n)=1+\sum_{r_x,r_y\in\{0,1\}}
  \left[N\left(\left\lfloor\frac{n-r_x}{2}\right\rfloor,
                 \left\lfloor\frac{n-r_y}{2}\right\rfloor\right)
       +2\left(\left\lfloor\frac{n-r_x}{2}\right\rfloor
              +\left\lfloor\frac{n-r_y}{2}\right\rfloor+1\right)\right].
```

For $`n\ge12`$, every coarse side parameter is at least five, so $`(*)`$ applies. If $`n=2k`$, the four pairs are $`(k,k),(k,k-1),(k-1,k),(k-1,k-1)`$. Their $`N`$ values sum to $`8k^2+1`$, and their pure terms sum to $`16k`$. Hence $`f(2k)=8k^2+16k+2`$.

If $`n=2k+1`$, all four pairs are $`(k,k)`$. Formula $`(*)`$ yields

```math
f(2k+1)=8k^2+24k+11+2(-1)^k.
```

These are exactly the claimed four residue formulas. The constant table has period four and not period one or two, since $`\gamma_1=3\ne-1=\gamma_3`$. The [cyclic-table and cyclotomic algorithms](minimal-periods-and-constant-term.md) both return $`(p_\beta,p_\gamma,p_f)=(1,4,4)`$ on the exact tables. They are applied **after** a sufficient modulus and the exact formulas have been proved.

The mechanism is two successive parity effects. First an anchor chooses a residue subgrid; then the diagonal-height lattice has index two inside its pair of height coordinates. On odd fine-window sizes, all four coarse rectangles have the same side parity, so this second parity survives. On even sizes the two coarse side parities cancel. The surviving oscillation is bounded, and does not multiply $`n`$.

## 8. How the general-tail obstruction is closed

The global-background separation in §3 cannot simply be used with different sector backgrounds. Replacing components by different left/right tails changes which periodic background is present. Appendix B retains all these choices as complete global fields $`H_{i,\sigma}=F_i+\sum_{k\ne i}T_{k,\sigma_k}`$, including their tangential phases.

Three facts supply the previously missing classification. A difference along a common tail period has finite-width nonzero components and recovers strong normal heights. Pure end labels are counted by projected tail-feasibility intervals; nonparallel directions prohibit a zero-width leading fiber on an interval of normal heights. Every remaining linear-sized collision has one strong direction and a bounded vertex correction. One first compares the actual global end fields, then proves the two nonzero corrections are at the same vertex, and only then centers them by that integral vertex. This order eliminates a possible hidden phase comparison between distinct vertices.

Finally, intrinsic strong status bounds the complete partner set of each discarded representation: a second collar height and the recovered strong height fix the anchor, and pure partners have only finitely many end labels. Thus no unclassified order-$`n`$ excess is hidden in an $`O(1)`$ exceptional list. The argument is direct in the original finite tables and does not assume that individual eliminated representative cells have constant linear coefficients.

## 9. Computation and reproducibility

The packaged diagnostics use Python's standard library and import only adjacent repository modules. Run `python3 -B scripts/verify.py` from the repository root. All packaged outputs are deterministic, with no timestamps or timings. Finite tests neither establish stabilization nor enlarge the theorem's scope.

[Exact pattern tools](../scripts/star_pattern_tools.py) count all patterns by a finite exhaustive representative list at each tested scale: every anchor with two active transition strips, solved from pairs of integer heights; one representative for every tangential phase and feasible tail-sign choice when exactly one strip is active; and representatives of every open direction sector and coordinate residue when no strip is active. The common-interior threshold proves that every feasible zero-active sign choice belongs to such an open sector. No finite translation box is declared exhaustive.

The [sparse-line audit](../scripts/integral_sparse_audit.py) evaluates the finite rational formula of §§2–5 before the associated counts. It checks nine sparse/background examples at fourteen scales each, independent support-by-support counts at two scales per example, and seven coarse rectangles. It reuses the [rational vertex-cone moments](../scripts/minimal_period_audit.py) for a second derivation of Theorem 2's exact tables, then applies cyclic and cyclotomic minimal-period reduction. Its [complete output](../scripts/integral_sparse_audit.json) reduces proved coefficient tables, rather than inferring a theorem from finite differences.

The [direct boundary-formula module](../scripts/integral_boundary_coefficients.py) implements Appendix B with exact rational arithmetic and full colored-cluster/end-field equality. The [general-tail audit](../scripts/integral_boundary_audit.py) compares the formula and all auxiliary finite data against 157 frozen source cases in [the packaged count fixture](../scripts/integral_boundary_cases.json). It compares these coefficients with the saved exact finite count tables, and freshly recounts a fixed coverage set spanning every polygon, component type, input period, direction set, and the nonzero-collision case. Its [complete output](../scripts/integral_boundary_audit.json) explicitly distinguishes fresh counts from the frozen tables. The implementation indexes identical nonzero clusters to avoid repeated full-collar scans; the formula and its equality tests are unchanged.

Separate exploratory screens are not bundled. They comprised 2,362,685 genuine star/window instances and 19,366,833 exact global pattern counts, allowing overlap between batches. The unequal-tail binary-table screens contained $`240^2=57,600`$ stars per geometry. One initially unstable three-direction sample was extended through $`n=40`$ and its finite diagnostic stabilized. Additional checks used colors modulo three and five and changed the collar, common difference period, integer window translate, and determinant-one complementary basis. These historical finite checks are not asymptotic proofs and are not reported as computations rerun by the packaged runner.

## 10. Sources and review status

The CN hypotheses are those of §0.1 of the [identified contest paper](../sources/bibliography.md#cn-target-contest-paper). The repository's [objects](objects.md), [area proof](all-star-extension.md), [eventual quasipolynomial proof](eventual-quasipolynomial.md), [representative coefficients](explicit-coefficients.md), [exact periods](minimal-periods-and-constant-term.md), and [full-support singleton framework](full-support-lines.md) supply the stated preceding definitions and results. The general-tail constancy proof is given in full below; none of those preceding results supplies its conclusion in advance.

This manuscript was developed by GPT-6.1 Sol and established in one fresh-context GPT-6.1 Sol adversarial review, with no fatal, proof-repair or presentation issue. The review record identifies the exact source and report fingerprints. No human review, proof-assistant certification or worldwide-priority claim is made.

## Appendix A. Finite-width defects on strict-vertex windows

**Theorem 3.** Let $`G_i`$ have period $`\kappa_i v_i`$ and be supported on finitely many heights of $`\pi_i`$, with at least one nonzero row. Let $`B`$ be a doubly periodic field and $`\theta=B+\sum_{i=1}^mG_i`$, with pairwise nonparallel primitive $`v_i`$ and $`m\ge2`$. Suppose $`P`$ is a positive-area integral convex polygon and the minimum and maximum of every $`\pi_i`$ on $`P`$ are attained at unique vertices. Equivalently, no edge of $`P`$ is parallel to any $`v_i`$. Then

```math
P_\theta(nP\cap\mathbb Z^2)=\Lambda(P)n^2+\beta n+O(1)
```

with a constant rational $`\beta`$. The following construction is a finite rational formula for it. This theorem permits arbitrary finite transverse defects and multiple parallel colored rows in each directional component. It retains a single global periodic background and the displayed window hypothesis.

### A.1 Strong profiles, weak profiles, and complete collision localization

Write the nonzero rows of component $`i`$ as a finite set $`A_i\subset\mathbb Z`$, with words $`w_{i,a}(s_i(z))`$ on $`\pi_i(z)=a`$. Delete identically zero rows. Put $`a_i^-=\min A_i`$, $`a_i^+=\max A_i`$, and let $`\kappa_i`$ be the least common period of the row words, reduced to its minimal value. Let $`N=\sum_i|A_i|`$ be the total number of physical supporting lines.

The outer band for visible component $`i`$ is

```math
a_i^--nu_i\le\pi_i(u)\le a_i^+-n\ell_i.
```

Choose an integer $`R_i`$ larger than the transverse profile width plus the short-section bound from §2 for $`2N\kappa_i`$ sites. Call an anchor **strong in direction $`i`$** when its height lies in this outer band at distance greater than $`R_i`$ from both endpoints. Every row of a strong component then has at least $`2N\kappa_i`$ consecutive geometric sites. Every other anchor in the band is **weak in direction $`i`$**. Outside the band the component is absent.

Each strong physical row has at least $`2N`$ nonzero sites, and at most $`N-1`$ can be altered by other tracks. Its relative affine line must therefore occur in any competing representation: all noncoincident competing tracks together cannot cover the remaining sites. All $`|A_i|`$ tracks of a strong profile are thus intrinsic. The competing profile has exactly $`|A_i|`$ possible parallel tracks, so the two translated finite height sets must be equal. A nonempty finite subset of the line has no nonzero translational period. Hence the two anchors have the same $`\pi_i(u)`$. In particular the competitor is strong in direction $`i`$ as well. Two strong directions determine the anchor.

For a weak component, all visible tracks have heights within a fixed distance of one extreme of the window. The unique-extreme-vertex hypothesis and §2 imply that its entire visible colored profile is supported on $`nq+C`$, where $`C`$ is one of finitely many bounded colored clusters in $`T_q`$. The finite list is indexed by a collar height and tangential phase. This is an exact eventual list. Integrality of $`q`$ makes the list independent of $`n`$ after translation.

Anchors weak in two directions form a uniformly bounded list, since they prescribe two independent heights from finite collar lists. Also remove anchors where a weak site meets another physical track; they again have the form $`z_*-nq-\delta`$, with $`z_*`$ in the finite physical intersection set. Thus, apart from bounded lists, a nonunique multi-component pattern has exactly one strong component and one nonempty weak colored cluster, disjoint from all tracks of the strong component.

These bounded lists have only boundedly many partners. A strong profile fixes its direction and transverse anchor height. A nonzero site of its off-profile cluster, together with a choice of another physical track through it, fixes the full anchor. There are at most $`N`$ choices. If there is no off-profile site, every non-pure partner belongs to one of the removed lists, and the pure phases are finite. A pattern with no strong direction cannot have a strong partner: the intrinsic full track set would fix the same height and hence the same strong status in the original representation. Its multi-component and pure weak representatives are on bounded lists. This proves the required closure under collision partners.

Unlike the one-track case, a weak cluster may have several nonzero sites and a single other direction may have several possible physical-track assignments. Those possibilities must be retained. Treating each row as an independent star direction would be incorrect.

### A.2 Invisible collar profiles and the raw coefficient

For the unique minimum vertex $`q`$ of direction $`i`$, enumerate collar states

```math
\pi_i(w)=a_i^+-d,\qquad d=0,\ldots,R_i,\qquad
s_i(w)\bmod\kappa_i=r.
```

The actual anchor is $`u=-nq+w`$. On the vertex cone its component pattern is the finite colored map

```math
C_{i,q,d,r}(z)=G_i(w+z),\qquad z\in T_q\cap\mathbb Z^2.
```

The map has finite support because $`q`$ is a strict extreme in direction $`i`$ and only finitely many physical heights are allowed. It depends only on the stated height and phase, not on the particular $`w`$ or on $`n`$. Use $`\pi_i(w)=a_i^-+d`$ at the maximum vertex. Define

```math
\eta_i^- =\frac1{\kappa_i}\sum_{d=0}^{R_i}\sum_{r=0}^{\kappa_i-1}
        \mathbf1[C_{i,q^-,d,r}=0],\qquad
\eta_i^+ =\frac1{\kappa_i}\sum_{d=0}^{R_i}\sum_{r=0}^{\kappa_i-1}
        \mathbf1[C_{i,q^+,d,r}=0].
```

Choosing a still larger collar does not change these sums: all additional states have a fully visible nonzero row. These weights combine geometric holes and phase-dependent invisibility.

Let $`K_i=\mathrm{lcm}(\kappa_i,b_i)`$ as in §3. The full visible profile in the interior determines its transverse height and joint line/background phase; weak pure heights are bounded in number. Hence pure modes contribute $`\sum_iK_iW_i n+O(1)`$. The same integer-support-line and arithmetic-progression argument as §4 gives

```math
\begin{split}
\mathfrak B_{\rm strip}=\sum_iK_iW_i+
\sum_{|I|\ge2}c_I\sum_{i\in I}\big[&
\tau_i^+(I)(a_i^++\tfrac12-\eta_i^-)\\ 
&+\tau_i^-(I)(-a_i^-+\tfrac12-\eta_i^+)\big].
\end{split}
```

For clarity, a component absent on a collar state removes a fraction of anchors on one integer line. With $`u=-nq+w_0+t v_i`$, its absence test has period $`\kappa_i`$ in $`t`$ and no $`n`$-dependence. Intersecting with the other bands gives length $`n\tau+O(1)`$; two collar deletions in distinct directions overlap in only boundedly many anchors. Inclusion-exclusion then gives the displayed raw coefficient and the usual area $`\Lambda`$.

### A.3 Whole-cluster labels: a finite formula for all linear collisions

Let $`\mathcal C_q`$ be the union of all nonempty colored maps $`C_{j,q,d,r}`$ in the finite collar catalogs. Identical maps are retained only once. For $`C\in\mathcal C_q`$, choose the lexicographically first site $`\delta`$ in its nonzero support.

Fix the strong direction $`i`$ and a candidate weak direction $`j\ne i`$. If the marked site belongs to physical row $`a\in A_j`$, the anchor branch is

```math
u=-nq+W_{ij,a,C}(t),\qquad
W_{ij,a,C}(t)=
 \frac{t v_j-(a-\pi_j(\delta))v_i}{d_{ij}}.
```

It satisfies $`\pi_i(W)=t`$ and $`\pi_j(W)=a-\pi_j(\delta)`$. Test integrality, then use the finite profile table to require that $`G_j(W+\cdot)`$ on $`T_q`$ be **exactly the entire colored map $`C`$**, not merely that it contain the marked site. Require the associated state to be in the weak collar catalog. A fixed representation has only one such row $`a`$ at its marked site, so this enumeration does not duplicate it.

Assign the label

```math
\lambda_{ij,a,C}(t)=
 \big(s_i(W_{ij,a,C}(t))\bmod\kappa_i,
        [T^{W_{ij,a,C}(t)}B]\big).
```

All tests and labels have period

```math
T_i=\mathrm{lcm}_{j\ne i}
     |d_{ij}|\mathrm{lcm}(\kappa_i,\kappa_j,Q).
```

Use the same rational intervals $`I_i(q)`$ and $`A_{ij}(q)`$ as §5. For each residue $`r`$ and $`\xi\in I_i(q)`$, collect the branches $`(j,a)`$ that pass the profile and integrality tests at $`t=r`$ and satisfy $`\xi\in A_{ij}(q)`$. Let their number minus their number of distinct labels be $`E_{i,q,C,r}(\xi)`$. Then set

```math
\Gamma_{\rm strip}=
\sum_q\sum_{C\in\mathcal C_q}\sum_i
 \frac1{T_i}\sum_{r=0}^{T_i-1}
       \int_{I_i(q)}E_{i,q,C,r}(\xi)\,d\xi.
```

This is again a finite rational sum of interval lengths times nonnegative integers. For a fixed $`q,C`$, equality of the strong profile forces the same integer $`t`$, and equality of its phase and the background is exactly equality of the labels: the common shift $`-nq`$ cancels. The weak colored cluster has already been matched in full. The short-intersection exclusions remove only boundedly many values of $`t`$. Unassigned directions differ from the leading inactive-band tests at only boundedly many $`t`$, because their relevant slopes $`d_{kj}/d_{ij}`$ are nonzero. Counting each phase class in each interval gives total excess $`n\Gamma_{\rm strip}+O(1)`$. Clusters at distinct vertices are disjoint for large $`n`$, and different strong-profile/cluster families cannot share a pattern.

Together with A.1 and A.2, this proves

```math
\boxed{\quad \beta=\mathfrak B_{\rm strip}-\Gamma_{\rm strip},\qquad p_\beta=1.\quad}
```

The construction terminates on the input tables. The explicit collar and cone thresholds give exact eventual catalogs; every subsequent test is finite integral or rational arithmetic. The zero-tail examples are also covered by the [direct boundary-formula module](../scripts/integral_boundary_coefficients.py); the finite formula is computed before comparison with the count tables.

The argument of A.1–A.3 alone leaves a geometric issue: at a long extreme edge a clipped multi-row profile can remain long while some tracks are absent. The following inward-extreme-row argument closes exactly that issue.

### A.4 Removing the strict-vertex restriction

**Theorem 4 (all equal-tail stars, all integral polygons).** Theorem 3 holds for every positive-area integral convex polygon, without the strict-extreme hypothesis. In the formula, an extreme edge contributes $`\eta_i^\pm=0`$ and no cluster catalog. Thus the same finite construction still computes $`\beta`$.

It suffices to close exactly the clipped-profile issue just identified. Retain the strong/collar classification of A.1. If a collar belongs to an extreme edge parallel to $`v_j`$, every visible physical row in that collar has a geometric section of length growing linearly with $`n`$, uniformly over its finite list of height offsets. Its periodic word is nonzero and therefore contributes a long intrinsic relative line, after any cancellations at other tracks.

At the **minimum** face of $`\pi_j`$ the outer anchor height is

```math
\pi_j(u)=a_j^+-n\ell_j-d,\qquad 0\le d\le R_j.
```

The inward extreme physical row $`a_j^+`$ is always visible: its relative height is $`n\ell_j+d`$. It is the **largest relative $`j`$-height** among all visible tracks of this component. All such visible tracks are long, so their relative lines are intrinsic in a competing representation. At the maximum face the analogous always-visible inward extreme row is $`a_j^-`$, and it is the smallest relative height. It follows at once that no phase can make the component invisible in an edge collar. This proves $`\eta_j^\pm=0`$ there.

Now consider a representation with a strong component $`i`$ and an edge-collar component $`j\ne i`$. Any partner has the same strong status and the same $`\pi_i(u)`$, by the full-track-set argument of A.1. It cannot be strong in direction $`j`$: strong status would then also be forced in the original representation. Nor can its $`j`$-profile lie in a vertex collar, whose entire support is bounded, while the observed relative $`j`$-lines have lengths tending to infinity. A collar at the opposite extreme has relative heights separated by $`nW_j+O(1)`$, and also cannot match. Hence the partner is in the **same edge collar**. Applying intrinsic-line recovery in both representations gives the same set of visible relative $`j`$-lines. Its largest height (at a minimum face), or smallest height (at a maximum face), recovers $`\pi_j(u)`$ from the known inward extreme physical row. The two nonparallel recovered heights determine the same anchor. This proves global uniqueness. Extra components change at most finitely many sites on each recovered line and do not alter the argument.

This argument uses the extremal **inward** row, not an assumption that every physical row remains visible. It is why a clipped long profile introduces no new order-$`n`$ collision family in the zero-tail setting.

If there is no strong direction and at least two visible components, the anchor lies in two finite collar-height lists and is on a bounded list. Its partners cannot acquire a strong direction, by A.1, so they are bounded as well. A pattern with only one visible edge-collar component belongs to the pure modes, with only boundedly many collar heights and phases. Consequently all remaining collision families of linear size have one strong direction and one **vertex** cluster, exactly those counted by A.3. The old bounded exceptions now remove only vertex-cluster intersection sites; edge-collar patterns with a strong partner are already globally unique. This completes the proof for arbitrary integral $`P`$.

Finally, suppose a CN star has equal periodic tails $`L_i=R_i`$ for each component. Put $`G_i=F_i-L_i`$ and $`B=\sum_iL_i`$, also including any initially separate periodic background. Each $`G_i`$ is supported on finitely many heights and has a common tangential period with $`F_i,L_i`$; these periods can always be enlarged. It is nonzero, since $`F_i`$ is not doubly periodic. The finite-tail sum $`B`$ is doubly periodic. Thus Theorem 4 applies. Unequal tails require retaining additional end-field information, as in Appendix B.

## Appendix B. General tails: complete proof and direct finite formula

This appendix proves the full assertion. It retains the original tail information after taking a common-period difference and gives a complete classification of the linear-sized collisions. Numerical stabilization is not a premise.

**Theorem 5 (constant linear coefficient for every CN star).** For every CN star and every positive-area integral convex polygon,

```math
P_\theta(nP\cap\mathbb Z^2)=\Lambda(P)n^2+\beta n+O(1)
```

with a constant rational $`\beta`$. In particular the eventual linear coefficient has minimal period one. The construction below is a finite formula directly in the component/tail tables and the polygon.

### B.1 Finite data and inward support extrema

Choose one positive integer $`q`$ such that every component has period $`qv_i`$ and all periodic tails have period lattice containing $`q\mathbb Z^2`$. This is possible by taking a common multiple of the directional periods and square-lattice tail periods. A separate periodic background is absorbed into one component. All fields are added in the original prime field.

Define the inward deviation extrema

```math
\alpha_i^- =\min\{\pi_i(z):F_i(z)\ne L_i(z)\},\qquad
\alpha_i^+ =\max\{\pi_i(z):F_i(z)\ne R_i(z)\}.
```

Both sets are nonempty, since $`F_i`$ is not doubly periodic. The first is bounded below and the second above. The extrema and their nonzero periodic row words are computable: inspect the finite transition strip and one normal tail period on each side. There is no requirement that $`\alpha_i^-\le\alpha_i^+`$; an empty-transition step can have $`\alpha_i^+=\alpha_i^--1`$.

The outer active band is

```math
V_i(n)=\{u:\alpha_i^--nu_i\le\pi_i(u)\le\alpha_i^+-n\ell_i\}.
```

Below this band the entire component window is exactly $`L_i(u+\cdot)`$; above it the entire window is exactly $`R_i(u+\cdot)`$. All constants in the rest of the proof are fixed from the input, before $`n`$ varies.

### B.2 A characteristic-$`p`$ difference detects strong directions

Choose $`h\in q\mathbb Z^2`$ nonparallel to every $`v_i`$, and put

```math
G_i(z)=F_i(z+h)-F_i(z),\qquad G=\sum_iG_i=\theta(\cdot+h)-\theta.
```

Each $`G_i`$ is supported on finitely many transverse rows and has period $`qv_i`$. It is nonzero: otherwise $`h`$ and $`qv_i`$ would be two independent periods of $`F_i`$. Let $`N`$ be the total number of nonzero supporting rows of all $`G_i`$. These are finite-width profiles of the kind used in Appendix A.

Equality of two $`\theta`$-patterns on $`S_n`$ implies equality of their $`G`$-patterns on

```math
E_n=\{s:s,s+h\in S_n\}.
```

If the facets of $`P`$ are $`A_js\le b_j`$, then $`E_n`$ is defined by $`A_js\le nb_j-\max(0,A_jh)`$. Its facet offsets differ from those of $`nP`$ by fixed constants. The short-section argument of §2 therefore still supplies fixed collars. More explicitly, choose a rational interior point $`c`$ of $`P`$ and a constant $`C`$ with

```math
C(b_j-A_jc)\ge\max(0,A_jh)\quad\text{for all }j.
```

Let $`E_n^{\mathbb R}`$ be the real polygon given by the displayed erosion inequalities, so $`E_n=E_n^{\mathbb R}\cap\mathbb Z^2`$. Then $`(n-C)P+Cc\subseteq E_n^{\mathbb R}`$. Concavity of its section lengths gives an explicit height collar outside which every specified $`G_i`$-row has at least $`2Nq`$ consecutive sites of $`E_n`$ at its relative height. This argument needs no assertion that the erosion's vertices are integral.

Choose collar widths $`R_i`$ so that, at any anchor in $`V_i(n)`$ farther than $`R_i`$ from both endpoints, every nonzero row of $`G_i`$ has at least $`2Nq`$ sites in $`u+E_n`$. Equivalently, its relative height is its physical height minus $`\pi_i(u)`$. Call direction $`i`$ strong at such an anchor; every anchor in the remaining part of $`V_i(n)`$ is in one of its two collars. The widths can also be enlarged to include any fixed finite list of row witnesses used below.

A strong $`G_i`$-row has at least $`2N`$ nonzero sites and loses at most $`N-1`$ to intersections with other gradient tracks. Its relative affine line must occur in any competing representation. All the nonzero supporting rows of $`G_i`$ are consequently intrinsic. The competing finite height set has the same cardinality and must contain this full set, so the two finite translated sets are equal. A finite nonempty height set has no nonzero translation period. Thus equality of $`\theta`$-patterns forces equality of $`\pi_i(u)`$ and strong status. Two strong directions force the same anchor.

This reduction is used for geometric recovery only. The lost periodic information is retained in the full end-profile labels below; equality of gradient patterns alone is never declared equality of original patterns.

### B.3 Vertex collars and edge collars for unequal tails

At a minimum vertex $`b`$ of direction $`j`$, put $`u=-nb+w`$, $`\pi_j(w)=\alpha_j^+-d`$, $`0\le d\le R_j`$. Subtract the bulk right tail. The component correction on the window is, eventually, the colored map

```math
C_{j,b,d,r}(z)=F_j(w+z)-R_j(w+z),\qquad z\in T_b\cap\mathbb Z^2,
\quad s_j(w)\equiv r\pmod q.
```

It is supported on $`0\le\pi_j(z)\le d`$, so is finite at a strict extreme vertex. At a maximum vertex use $`\pi_j(w)=\alpha_j^-+d`$ and subtract $`L_j`$; its support is similarly finite. The tables and exact eventual cone equality are computed as in §2. Every map depends on the fixed height and phase, with no $`n`$ in its color test. A vertex-collar component is called active when this map is nonzero. A strong component is active. Outside $`V_i(n)`$ it is inactive. In an edge collar it is always active for sufficiently large $`n`$.

The last assertion, and the required uniqueness there, have an exact gradient proof. Write $`d_h=\pi_j(h)\ne0`$. The global inward gradient extrema are

```math
\max\pi_j(\mathrm{supp}G_j)
=\alpha_j^++\max(0,-d_h),\qquad
\min\pi_j(\mathrm{supp}G_j)
=\alpha_j^--\max(0,d_h).
```

For example, in $`F_j-R_j`$ the highest nonzero row is $`\alpha_j^+`$. If $`d_h>0`$, only the unshifted term contributes to that highest row of its difference; if $`d_h<0`$, only the shifted term contributes at height $`\alpha_j^+-d_h`$. Thus there is no cancellation at the stated extreme. The lower formula is identical using $`F_j-L_j`$.

On a minimum edge the erosion moves its support line inward by $`\max(0,-d_h)`$ in this coordinate. Consequently the inward gradient extreme row is visible at exactly the same collar depth $`d`$ as the original inward deviation row. Its section has length proportional to $`n`$ and is intrinsic. The analogous statement holds at a maximum edge.

If a pattern has a strong direction $`i`$ and an active edge-collar direction $`j\ne i`$, any partner has the same $`\pi_i(u)`$ and the same strong status. The partner cannot be strong in direction $`j`$, cannot have only a bounded vertex-collar $`j`$-signal, and cannot lie at the opposite extreme edge. In the same edge collar, the set of visible long gradient $`j`$-lines is intrinsic in both representations. Its inward extreme recovers $`\pi_j(u)`$ by the displayed formula. The two independent heights recover the anchor. Hence such patterns are globally unique.

Two collar constraints in different directions prescribe heights from finite lists and determine only boundedly many anchors, whether their corrections are visible or invisible. Apart from bounded lists, every nonunique multi-active pattern therefore has exactly one strong direction $`i`$ and exactly one nonzero **vertex** correction. Edge collars contribute neither invisibility weights nor a linear collision excess.

Define $`\eta_i^-`$ and $`\eta_i^+`$ by averaging the indicator $`[C_{i,b,d,r}=0]`$ over $`r\bmod q`$ and summing $`d=0,\ldots,R_i`$ at the corresponding strict vertex. Put $`\eta=0`$ at an extreme edge. The sum stabilizes if the collar is enlarged, since the inward extreme row is nonzero in every phase on sufficiently long sections.

The same integer-support and deleted-line calculation as A.2 now gives, for anchors active in all directions of $`I`$,

```math
\begin{split}
N_I(n)=\mathrm{area}(K_I)n^2+n\sum_{i\in I}\big[&
\tau_i^+(I)(\alpha_i^++\tfrac12-\eta_i^-)\\ 
&+\tau_i^-(I)(-\alpha_i^-+\tfrac12-\eta_i^+)\big]+O(1).
\end{split}
```

A negative value of $`\alpha_i^+-\alpha_i^-`$ causes no problem: the growing band still has positive width for large $`n`$. All support constants are integral. Intersections of collar deletions in different directions are bounded. Thus inclusion-exclusion still counts every multi-active anchor once, with quadratic coefficient $`\Lambda`$ and a constant linear coefficient.

### B.4 Full end profiles and the pure-mode coefficient

For a direction $`i`$, choose a tail sign $`\sigma_k\in\{L,R\}`$ for every $`k\ne i`$, and define the global one-direction field

```math
H_{i,\sigma}=F_i+\sum_{k\ne i}T_{k,\sigma_k},\qquad
\lambda_{i,\sigma,r}=[T^{r v_i}H_{i,\sigma}],\quad 0\le r\lt q.
```

Brackets mean exact equality of fields on the whole lattice. These labels are finite and decidable from the input: inspect the transition rows of $`F_i`$, one complete normal tail period on each side, and $`q`$ tangential phases. Any two labels have a difference periodic along $`qv_i`$. If the difference is nonzero anywhere, it is nonzero on at least one of this finite list of transverse witness rows, repeated along that row.

At a fixed normal anchor height $`a`$, write $`u=a e_i+t v_i`$. The pure field is $`T^{ae_i}\lambda_{i,\sigma,t\bmod q}`$. The common normal translation cancels in comparisons, so equality of global labels has no dependence on $`a`$ or $`n`$.

It remains to count which labels are realizable. Put $`\xi=a/n`$ and $`\tau=t/n`$. Let $`J_{i,\sigma}`$ consist of the $`\xi\in[-u_i,-\ell_i]`$ for which there is a nonempty open interval of $`\tau`$ satisfying, for every $`k\ne i`$,

```math
\begin{cases}
\pi_k(e_i)\xi+\pi_k(v_i)\tau<-u_k,&\sigma_k=L,\\ 
\pi_k(e_i)\xi+\pi_k(v_i)\tau>-\ell_k,&\sigma_k=R.
\end{cases}
```

This is an interval projection of rational linear inequalities. It is computed by comparing each affine lower bound on $`\tau`$ with each affine upper bound. Two such bounds from different directions have distinct slopes in $`\xi`$: equality of their ratios $`\pi_k(e_i)/\pi_k(v_i)`$ would make the two directions parallel. A leading fiber of zero width can therefore occur at only finitely many $`\xi`$, and cannot persist on an interval of normal heights.

The exact inactive inequalities differ from these leading inequalities by fixed constants. Except at boundedly many integer heights $`a`$ near the finitely many breakpoints, a positive-width fiber contains representatives of **every** residue $`t\bmod q`$, and an empty leading fiber contains none. This follows by counting an arithmetic progression in an interval whose endpoint difference is affine in $`(a,n)`$; away from each root by a fixed number of integer heights its length exceeds $`q`$ plus the fixed endpoint perturbations. Unbounded intervals have every residue automatically.

Define

```math
L_i(\xi)=\#\{\lambda_{i,\sigma,r}:\xi\in J_{i,\sigma},\ 0\le r\lt q\},
\qquad
\boxed{\quad C_{\rm pure}=\sum_i\int_{-u_i}^{-\ell_i}L_i(\xi)\,d\xi.\quad}
```

The integrands are piecewise constant integer functions with rational breakpoints. Thus this is a finite rational expression.

This counts actual pure patterns with only bounded errors, not merely formal tail choices. For a strong normal height away from a fixed larger endpoint collar, the window contains a full period of every witness row for any two distinct end labels. Different labels therefore give different patterns. Strong heights are recovered intrinsically from the gradient; different heights, or different strong directions, cannot give the same pure pattern. The excluded endpoint heights give only boundedly many pure labels.

An anchor with one active component and an invisible vertex-collar component is also covered, apart from bounded exceptions. Translate it along the active direction by a bounded multiple of $`qv_i`$ until the invisible component lies fully on its bulk-tail side. The active component and every tail field are preserved. All other components remain in their previous tail regions unless one is within a fixed distance of a band boundary. Such an exception imposes two independent collar-height conditions and occurs at only boundedly many anchors. Two invisible collars are already on the bounded list. Thus every remaining pure pattern has a representative with all other components fully in their stated inactive tails, exactly the realizability counted above.

Combining pure modes, the bounded number of purely periodic backgrounds, and the multi-active anchors yields raw representations with linear coefficient

```math
\boxed{\begin{split}
\mathfrak B_{\rm all}=C_{\rm pure}+
\sum_{|I|\ge2}c_I\sum_{i\in I}\big[&
\tau_i^+(I)(\alpha_i^++\tfrac12-\eta_i^-)\\ 
&+\tau_i^-(I)(-\alpha_i^-+\tfrac12-\eta_i^+)\big].
\end{split}}
```

### B.5 Whole-cluster collisions with complete tail labels

Let $`\mathcal C_b`$ be the finite set of distinct nonzero vertex corrections from B.3 at vertex $`b`$. Fix $`C\in\mathcal C_b`$ and its lexicographically first nonzero site $`\delta`$. A candidate direction $`j`$ must have $`b`$ as its strict minimum or maximum; its bulk sign $`\sigma_j`$ is respectively $`R`$ or $`L`$. Let $`c`$ run over the finite physical collar heights $`\alpha_j^+-d`$ or $`\alpha_j^-+d`$.

With strong direction $`i\ne j`$ and integer parameter $`t`$, the anchor branch is

```math
u=-nb+W_{ij,c}(t),\qquad
W_{ij,c}(t)=\frac{t v_j-cv_i}{d_{ij}}.
```

Test integrality and require the entire colored map

```math
F_j(W_{ij,c}(t)+z)-T_{j,\sigma_j}(W_{ij,c}(t)+z)
```

on $`T_b`$ to equal $`C`$. This is a lookup in the finite collar catalog. For the remaining components choose signs $`\sigma_k`$. The branch is allowed in the leading interval

```math
\begin{split}
A_{ij,\sigma}(b)=\{\xi\in I_i(b):\ &
\pi_k(b)-\tfrac{d_{kj}}{d_{ij}}\xi>u_k\quad\text{if }\sigma_k=L,\\ 
&\pi_k(b)-\tfrac{d_{kj}}{d_{ij}}\xi<\ell_k\quad\text{if }\sigma_k=R,
\quad k\notin\{i,j\}\}.
\end{split}
```

Assign the full end-profile label

```math
\lambda_{i,\sigma,\,s_i(W_{ij,c}(t))\bmod q}.
```

All integrality, cluster, and label tests are periodic in $`t`$ with period

```math
T_i=q\mathrm{lcm}_{j\ne i}|d_{ij}|.
```

For a residue $`r\bmod T_i`$ and $`\xi\in I_i(b)`$ collect all branches $`(j,c,\sigma)`$ passing these tests and with $`\xi\in A_{ij,\sigma}(b)`$. Let their number minus their number of distinct full end labels be $`E_{i,b,C,r}(\xi)`$. Define

```math
\boxed{\quad
\Gamma_{\rm all}=\sum_b\sum_{C\in\mathcal C_b}\sum_i
 \frac1{T_i}\sum_{r=0}^{T_i-1}
       \int_{I_i(b)}E_{i,b,C,r}(\xi)\,d\xi.\quad}
```

Only rational interval endpoints and finite equality tests occur. Tail choices for a fixed branch are mutually exclusive in the interiors of these intervals. Two different weak directions may name the same anchor only where two collar constraints meet, already on a bounded list.

To prove the formula, first enlarge the list of excluded strong-height endpoints by a **fixed** amount so that every nonzero difference of two end labels is observed at more than twice the maximum vertex-cluster size. The finite row witnesses of B.4 and the section-length bound make this possible. This does not require changing any collar catalog: a newly excluded strong height combined with a weak collar in another direction gives only boundedly many multi-active anchors; for pure modes it excludes only boundedly many normal heights and labels.

Outside this bounded list, compare first in the **original window coordinates**. The common intrinsic normal anchor height is $`a`$. The two end fields are $`T^{ae_i}\lambda_{i,\sigma,s_i(u)\bmod q}`$. A nonzero difference repeats on a witness row too often to be canceled by the two bounded vertex clusters, so these actual global end fields must agree. Subtraction then leaves the two corrections supported on $`nb_1+C_1`$ and $`nb_2+C_2`$. Nonzero corrections at distinct vertices cannot agree for large $`n`$. Hence $`b_1=b_2=b`$, and the entire maps $`C_1,C_2`$ agree. Only after this step do we center both representations at this same vertex. The common shift by $`nb`$ cancels; their global end-field equality is exactly the centered-label test used in the formula. This ordering also rules out a phase condition connecting two different vertices. Conversely equal centered labels and equal $`C`$ give identical original patterns. In the centered coordinates, a candidate pattern is

```math
H_{i,\sigma}(W_{ij,c}(t)+z)+C(z),\qquad z\in n(P-b)\cap\mathbb Z^2.
```

The common normal coordinate is $`t`$. Subtracting its common translation $`te_i`$ proves that the displayed end label is the exact equality test. Translating back by $`nb`$ preserves equality and contributes no $`n`$-phase. This is the missing information that a gradient-only comparison would discard.

Exact strong and inactive inequalities differ from $`t/n\in A_{ij,\sigma}(b)`$ at only boundedly many integers $`t`$. All excluded-direction slopes $`d_{kj}/d_{ij}`$ are nonzero. The possible extra invisible-collar cases also prescribe only boundedly many $`t`$. On each partition interval $`(a,b)`$, a phase $`t\equiv r\pmod{T_i}`$ contributes $`n(b-a)/T_i+O(1)`$ representatives. Grouping by identical end labels therefore gives excess $`n\Gamma_{\rm all}+O(1)`$.

There is no unbounded missing partner set in this argument. Any pattern with a strong direction has that direction and its normal height intrinsic. Every multi-active partner without a second strong direction has another active component at one of finitely many collar heights; choosing its direction and that height fixes the anchor. Thus the number of such partners is uniformly bounded by the sum of the finite collar-list sizes. Its pure partners have only the finite end labels at the fixed strong height. A pattern without a strong direction has no strong partner; its multi-active representations prescribe two collar heights, and its pure weak representatives have finitely many heights and phases. Consequently every discarded bounded representation list has a bounded closure under equality. This justifies all $`O(1)`$ exclusions, including their partners.

Different surviving strong-direction/vertex/correction families cannot share a pattern, by intrinsic strong status, recovered height, and the preceding end-field separation. This proves that the integral counts the **entire** linear representation excess, not merely a collection of possible collisions.

### B.6 Conclusion and exact role of integrality

The complete formula is

```math
\boxed{\quad\beta=\mathfrak B_{\rm all}-\Gamma_{\rm all}.\quad}
```

Every entry is a finite rational sum in the original finite tables and rational polygon geometry. The raw quadratic coefficient is $`\Lambda`$ by pointwise inclusion-exclusion of the direction bands; pure modes and collision families have only linear size. The bounded exceptional sets and their partners were controlled in B.5. Hence the remainder after this constant $`\beta`$ is bounded.

Integrality is used in the collar templates: the physical anchor is $`-nb+w`$ with $`nb\in\mathbb Z^2`$ for every $`n`$. Its residual colors, cluster equality and full end-label tests have no residual $`n`$-phase. All remaining phases are averaged along integer intervals. On a rational polygon the vertex-cone catalogs and integer support displacements must instead be indexed by a dilation residue, and the proof does not assert constant $`\beta`$ there.

The established all-star eventual quasipolynomial theorem now implies that every eventual $`\beta_r`$ equals this $`\beta`$: a nonzero coefficient difference could not have bounded remainder on its infinite residue progression. Thus its minimal period is one, while $`\gamma`$ may retain a nontrivial period. Theorem 2 gives an exact integral-window example with $`p_\gamma=4`$.

Equivalently, for any proved sufficient modulus $`M`$ in the [representative-cell algorithm](explicit-coefficients.md), its exact linear table must satisfy

```math
\mathcal B(X)=\beta(1+X+\cdots+X^{M-1}).
```

Every cyclotomic factor $`\Phi_d`$ with $`d\mid M`$, $`d>1`$, therefore cancels from the linear table. The exact minimal-period algorithm returns $`p_\beta=1`$ and $`p_f=p_\gamma`$. This is a structural cancellation statement about the sum of all cells; individual cells need not have constant linear coefficients.

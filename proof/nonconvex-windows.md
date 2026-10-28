# Directional projections, nonconvex windows, and exact complexity of stars

**Status: established in one fresh-context AI adversarial review.** The review found no fatal or proof-repair issue; its dependency wording clarification is incorporated. Source and review fingerprints are in the [review-status record](../sources/review-status.md). This is AI review, not human certification or proof-assistant verification.

We extend the existing convex-window area theorem for the star configurations of the specified Convex Nivat paper to finite unions of full-dimensional convex bodies. The window may be disconnected and may have holes. Its projections must be retained as unions of intervals: replacing them by their convex hulls can give the wrong quadratic coefficient. We also characterize equality between that coefficient and the area of the window. A separate exact formula for one fixed two-line star on every finite window identifies a genuine obstruction to including zero-width appendages.

The strip-signal and finite-marker constructions are inherited from the contest paper and the [convex area proof](all-star-extension.md), and are proved again below. The new scope is the nonconvex geometry, the equality characterization, and the arbitrary-window collision formula. Section 6.1 additionally resolves the product-region exception for two components supported on periodic affine lines, allowing zero entries in their words. We do not reopen the convex-window work on general linear/constant coefficients, minimal quasiperiods, or Presburger counting.

## 1. Statements

Let

```math
\theta=\sum_{i=1}^{m}F_i:\mathbb Z^2\longrightarrow\mathbb F_p,
\qquad m\ge2,
```

be a star in exactly the sense of CN §0.1. Thus the primitive integer directions $`v_i`$ are pairwise nonparallel, $`\pi_i(z)=\det(v_i,z)`$, and each $`F_i`$:

- has a period $`k_iv_i`$, with $`k_i\ge1`$;
- is not doubly periodic;
- agrees with doubly periodic fields $`L_i,R_i`$ on $`\pi_i<\ell_i`$ and $`\pi_i>r_i`$, respectively, where $`\ell_i\le r_i+1`$.

The strip may be empty, the two tails may differ, tangential periods may exceed one, and a doubly periodic background may be absorbed into a component. Pattern complexity is

```math
P_\theta(S)=\#\{(\theta(u+s))_{s\in S}:u\in\mathbb Z^2\}.
```

Area is normalized by $`\mathrm{area}([0,1]^2)=1`$.

**Definition.** In Theorems 1–2, a region is a set

```math
P=K_1\cup\cdots\cup K_t,
```

where $`t`$ is finite and positive and every $`K_a`$ is a compact convex subset of $`\mathbb R^2`$ with nonempty interior. Neither disjointness nor rationality is required. In particular, every compact polygonal region obtained by taking the closure of its interior, with finitely many polygonal boundary components, is included: triangulate it into nondegenerate closed triangles. Isolated points and zero-width appendages are not included by this definition.

Put

```math
A_i=\pi_i(P),\qquad B_i=\pi_i^{-1}(A_i),\qquad
\mathcal M(P)=\bigcup_{i\lt j}(B_i\cap B_j),\qquad
\Lambda(P)=\mathrm{area}(\mathcal M(P)).
\qquad\text{(1)}
```

Each $`A_i`$ is a finite union of nondegenerate compact intervals. Each $`B_i`$ is a union of parallel strips, and $`\mathcal M(P)`$ is a finite union of bounded parallelograms. The gaps in $`A_i`$ are part of the definition.

**Theorem 1 (nonconvex directional-projection theorem).** For every fixed CN star and every fixed region $`P`$ as above,

```math
\boxed{P_\theta(nP\cap\mathbb Z^2)=\Lambda(P)n^2+O_{\theta,P}(n).}
\qquad\text{(2)}
```

The coefficient depends only on the directions and their actual projection sets $`A_i`$. The result includes both Cases A and B of CN, arbitrary pair determinants, unequal tails and nonzero periodic backgrounds. It is an exact leading coefficient with matching upper and lower bounds.

**Theorem 2 (area equality and strictness).** Under the same geometric hypotheses:

$`\Lambda(P)\ge\mathrm{area}(P)`$ always. If $`m\ge3`$, the inequality is strict. If $`m=2`$, equality holds if and only if

```math
P=\{x:\pi_1(x)\in A_1,\ \pi_2(x)\in A_2\}.
\qquad\text{(3)}
```

Thus the exceptional regions are precisely the images of Cartesian products of two finite interval unions under the inverse of the linear isomorphism $`x\mapsto(\pi_1(x),\pi_2(x))`$.

For connected $`P`$, every $`A_i`$ is an interval and

```math
\Lambda(P)=\Lambda(\mathrm{Conv}(P)).
\qquad\text{(4)}
```

For two directions a connected region satisfies equality in Theorem 2 precisely when it is a parallelogram whose sides are parallel to the two directions. Consequently, except for the two-direction product regions (3),

```math
P_\theta(nP\cap\mathbb Z^2)
=|nP\cap\mathbb Z^2|
+\bigl(\Lambda(P)-\mathrm{area}(P)\bigr)n^2+O(n),
\qquad\text{(5)}
```

with a strictly positive quadratic surplus. In particular $`P_\theta(nP\cap\mathbb Z^2)>|nP\cap\mathbb Z^2|`$ for all sufficiently large $`n`$ in this range.

The statement does not resolve the sign of the lower-order difference in all exceptional product regions. Theorem 4 in §6.1 does resolve it for all two-line stars with zero background, arbitrary nonzero periodic words, and arbitrary pair determinant. Neither result is a periodicity theorem for arbitrary configurations or arbitrary nonconvex finite windows.

**Corollary 2.1 (all fixed connected regions).** For every CN star and every connected region $`P`$ in the sense above,

```math
P_\theta(nP\cap\mathbb Z^2)\ge |nP\cap\mathbb Z^2|+1
\quad\hbox{for all sufficiently large }n.
```

All nonconvex connected regions follow from the strict quadratic surplus. Only the two-direction aligned parallelogram case uses CN's existing Theorem T, since those windows are lattice-convex. Thus the corollary explicitly includes one existing theorem as an input; Theorems 1–2 themselves do not.

In fact, writing $`C=\mathrm{Conv}(P)`$, every connected nonconvex region satisfies

```math
P_\theta(nP\cap\mathbb Z^2)-|nP\cap\mathbb Z^2|
\ge (\mathrm{area}(C)-\mathrm{area}(P))n^2-O(n),
```

with a positive displayed coefficient. This follows from (4), $`\Lambda(C)\ge\mathrm{area}(C)`$, and (10). A closed regular region properly contained in its convex hull has a positive-area deficit by the same regular-closed-set argument used in §6.

For an independent finite-window result, let

```math
\theta_\times(x,y)=\mathbf1[x=0]+\mathbf1[y=0]\pmod2.
\qquad\text{(6)}
```

Given a finite nonempty $`S\subset\mathbb Z^2`$, form its bipartite incidence graph $`G_S`$. The row vertices are the distinct $`y`$-coordinates, the column vertices are the distinct $`x`$-coordinates, and each site $`(x,y)\in S`$ is the edge joining its row and column. Let $`r,c`$ be its numbers of row and column vertices. Let $`e`$ be the number of components with one row and one column, and let $`t`$ be the total number of components with row/column counts $`(1,2)`$, $`(2,1)`$ or $`(2,2)`$.

**Theorem 3 (exact arbitrary-window count for the cross star).**

```math
\boxed{P_{\theta_\times}(S)=(r+1)(c+1)-2e-\binom e2-2t.}
\qquad\text{(7)}
```

This counts all plane translates, including collisions caused by missing sites. Its proof is in §7. The resulting counterexample in §8 shows that positive area, connectedness and an $`O(n)`$ change in the number of sites do not suffice to extend Theorem 1 to arbitrary compact sets.

## 2. Two geometric estimates

We first record the elementary estimates required by (2).

**Lemma 2.1 (long observed rows in an eroded convex piece).** Let $`K`$ be a fixed compact convex body with nonempty interior, let $`v`$ be primitive, and put $`\pi=\det(v,\cdot)`$. For every finite nonempty $`E\subset\mathbb Z^2`$ and integer $`k\ge1`$, there are constants $`C,N`$ such that, for every $`n\ge N`$ and every integer

```math
n\min_K\pi+C\le h\le n\max_K\pi-C,
\qquad\text{(8)}
```

the set $`\{s\in\mathbb Z^2:s+E\subset nK,\ \pi(s)=h\}`$ contains $`k`$ consecutive lattice points in direction $`v`$.

**Proof.** Complete $`v`$ to an integer basis $`(v,w)`$ with determinant one. Chord lengths are measured in $`v`$ units. Choose $`q\in\mathrm{int}K`$, and write $`a=\pi(q)`$, $`l=\min_K\pi`$, $`b=\max_K\pi`$. The chord-length function is concave: convex combinations of two chord segments lie in the intervening section. Since it is positive at $`a`$, there is $`c_0>0`$ such that its value at any $`h\in[l,b]`$ is at least $`c_0\min(h-l,b-h)`$.

Let $`E'=E+\{0,v,\ldots,(k-1)v\}`$. Choose a fixed $`\lambda>0`$ with $`q+e/\lambda\in K`$ for all $`e\in E'`$. Convexity gives

```math
(n-\lambda)K+\lambda q\subseteq\bigcap_{e\in E'}(nK-e)
\qquad(n>\lambda).
```

Its projection endpoints differ from $`nl,nb`$ by constants. The chord estimate shows that, after excluding a constant number of heights at each end, every integer-height chord has length at least one. It therefore contains a lattice point $`s`$. Then $`s+E'\subset nK`$, giving the required run. $`\square`$

**Lemma 2.2 (interval-union counting).** Fix finitely many pairwise nonparallel primitive forms $`\pi_i`$, and for each $`i`$ a finite list of intervals $`[a_{ih},b_{ih}]`$ with $`a_{ih}\lt b_{ih}`$. Allow fixed, possibly different, additive constants at every endpoint of their dilates. Let $`U_n`$ be the set where at least two of the resulting families of strips are present. If $`U`$ is the corresponding unshifted set at scale one, then

```math
|U_n\cap\mathbb Z^2|=\mathrm{area}(U)n^2+O(n).
\qquad\text{(9)}
```

Intervals within one family are counted as the same direction, not as distinct directions.

**Proof.** Express $`U_n`$ as the finite union, over pairs of distinct directions and choices of one interval in each, of their parallelogram intersections. The vertices of each such parallelogram are $`n`$ times the unshifted vertices plus $`O(1)`$, by solving two fixed independent linear equations. Its symmetric difference from the corresponding dilate lies in a constant-width neighborhood of a polygonal boundary of length $`O(n)`$. That neighborhood contains $`O(n)`$ lattice points: cover the segments by $`O(n)`$ fixed-radius disks. A finite union has the same error bound. Finally, lattice count differs from area by $`O(n)`$ for the unshifted union of dilated parallelograms, since only unit squares meeting its finitely many boundary segments can contribute a discrepancy. Overlaps or coincident endpoints do not change this argument. $`\square`$

We will also use

```math
|nP\cap\mathbb Z^2|=\mathrm{area}(P)n^2+O(n)
\qquad\text{(10)}
```

and the fact that a fixed-width neighborhood of $`\partial(nP)`$ contains $`O(n)`$ lattice points. Indeed, $`\partial P\subseteq\bigcup_a\partial K_a`$; each convex body has finite perimeter, and its dilated boundary is covered by $`O(n)`$ fixed-radius disks. The unit-square argument then proves (10). No polygonal or rational boundary is needed here.

## 3. Signals for every star, without Case B

Choose $`H_i=\kappa_i v_i`$, with $`\kappa_i\ge1`$, which is a period of $`F_i`$ and of every tail field $`L_j,R_j`$. Such multiples exist because the tail period lattices have finite index. Fix any injective map $`w:\mathbb F_p\to\mathbb Z`$, and set

```math
\eta=w(\theta),\qquad
Q_i=\prod_{j\ne i}(T^{H_j}-1),\qquad f_i=Q_i\eta,
\quad (T^h\eta)(z)=\eta(z+h).
\qquad\text{(11)}
```

All these operators act in characteristic zero. Let $`E_i`$ be the set of all formal subset sums of the $`H_j`$ with $`j\ne i`$. It contains $`0`$ and contains the support of $`Q_i`$, even if some nonextreme terms cancel. A window on $`s+E_i`$ determines $`f_i`$ at $`s`$.

**Lemma 3.1 (strip signals with common end extrema).** For every $`i`$:

1. $`f_i`$ has support in a fixed bounded-width strip $`\alpha_i\le\pi_i\le\gamma_i`$.
2. Outside a bounded core in that strip, it agrees with one of two nonzero $`H_i`$-periodic end profiles $`g_i^+,g_i^-`$.
3. Both end profiles have the same least and greatest nonzero transverse heights, denoted $`\mu_i^-,\mu_i^+`$.

In particular $`f_i\ne0`$. Outside the core every nonzero value has transverse height in $`[\mu_i^-,\mu_i^+]`$.

**Proof.** Translation by $`nH_i`$ as $`n\to\pm\infty`$ isolates

```math
\Xi_i^\pm=F_i+B_i^\pm,
```

where $`B_i^\pm`$ is a sum of the appropriate other tails. Stabilization holds on each finite set: $`H_i`$ preserves its own component and all tails, and $`\pi_j(H_i)\ne0`$ for $`j\ne i`$. Put $`g_i^\pm=Q_iw(\Xi_i^\pm)`$. These fields are $`H_i`$-periodic.

Let

```math
a_i=\min\pi_i\{z:F_i(z)\ne L_i(z)\},\qquad
b_i=\max\pi_i\{z:F_i(z)\ne R_i(z)\}.
```

The indicated extrema exist, and the sets are nonempty: otherwise $`F_i`$ would equal a doubly periodic tail everywhere. Because addition in $`\mathbb F_p`$ can be canceled and $`w`$ is injective,

```math
\mathrm{supp}\bigl(w(F_i+B)-w(L_i+B)\bigr)=\{F_i\ne L_i\}
```

for either background $`B=B_i^\pm`$, and likewise with $`R_i`$. Each pure encoded tail is killed by $`Q_i`$.

In the formal expansion of $`Q_i`$, the unique greatest projected exponent is obtained by selecting exactly those $`H_j`$ with $`\pi_i(H_j)>0`$; the unique least is obtained from those with $`\pi_i(H_j)<0`$. Write their projected values as $`\beta_i^+,\beta_i^-`$. At the lowest possible output row of the left-tail difference only the greatest term survives; at the highest possible output row of the right-tail difference only the least survives. Hence, for both signs,

```math
\min\pi_i(\mathrm{supp}g_i^\pm)=a_i-\beta_i^+=\mu_i^-,\qquad
\max\pi_i(\mathrm{supp}g_i^\pm)=b_i-\beta_i^-=\mu_i^+.
\qquad\text{(12)}
```

This proves nonvanishing and the common extrema, including when the two end profiles have different values.

To control the global $`f_i`$, widen each original transition strip enough to accommodate all shifts in $`E_i`$. Outside that widened strip, the entire template for the component is in one pure tail. If at most one widened strip is active, and its index is $`j\ne i`$, pair terms of $`Q_i`$ using its factor $`T^{H_j}-1`$. The varying component $`F_j`$ and every other selected tail are preserved by $`H_j`$, so the terms cancel. The same holds with no active component. The region with at least two active widened strips is bounded, since the directions are nonparallel. Thus, outside a bounded set, support of $`f_i`$ lies in its own widened strip. Enlarging this strip accommodates the bounded exception.

Write points of this strip as $`sv_i+tw_i`$ in an integer basis, with $`t`$ in a fixed finite interval. As $`s\to\pm\infty`$, the other components on the entire template agree with exactly the tails defining $`B_i^\pm`$. This is uniform over the finitely many heights. Consequently $`f_i=g_i^\pm`$ at the corresponding ends, with only a bounded core remaining. $`\square`$

This is the [all-star strip argument](all-star-extension.md), with its dependence on the two tails made explicit. In Case A we have not asserted $`Q_i\eta`$ to be globally periodic.

## 4. Globally unique anchors on a nonconvex window

For $`S_n=nP\cap\mathbb Z^2`$, put

```math
E_i(S_n)=\{s\in\mathbb Z^2:s+E_i\subset S_n\},\qquad
Y_i(u)=\{\pi_i(s):s\in E_i(S_n),\ f_i(u+s)\ne0\}.
\qquad\text{(13)}
```

Both sets of observations, including the values, are determined by the pattern at anchor $`u`$.

For piece $`K_a`$ write $`a_{ia}=\min_{K_a}\pi_i`$ and $`b_{ia}=\max_{K_a}\pi_i`$. Apply Lemma 2.1 with template $`E_i`$ and run length $`\kappa_i`$, and use one sufficiently large constant $`C`$ for the finitely many pairs $`(i,a)`$. Define the shrunken projected-anchor sets

```math
J_i^-(n)=\bigcup_a
[\mu_i^+-nb_{ia}+C,\ \mu_i^--na_{ia}-C]
\qquad\text{(14)}
```

and

```math
\mathcal G_n=\{u\in\mathbb Z^2:\pi_i(u)\in J_i^-(n)
\text{ for at least two distinct }i\}.
\qquad\text{(15)}
```

Intervals are nonempty for all sufficiently large $`n`$. If a direction is good in (14), one convex piece supplies complete $`\kappa_i`$-point runs on both local heights $`\mu_i^--\pi_i(u)`$ and $`\mu_i^+-\pi_i(u)`$ inside the erosion. Different good directions may use different pieces.

Lemma 2.2 gives

```math
|\mathcal G_n|=\Lambda(P)n^2+O(n).
\qquad\text{(16)}
```

Reflection of the projection intervals gives $`-\mathcal M(P)`$ in this anchor count, which has the same area. No convexification of a projection gap has occurred.

**Lemma 4.1 (distant anchors are globally unique).** There is a fixed radius $`R`$ such that, if $`u\in\mathcal G_n`$ and $`(u+nP)\cap B_R(0)=\varnothing`$, its $`S_n`$-pattern occurs at no other anchor in the entire plane.

**Proof.** First choose a common core radius $`R_0`$ for Lemma 3.1. Enlarge the cores if necessary so that, at either extreme height, a run of $`\kappa_i`$ consecutive lattice points avoiding the core cannot visit both ends: in the $`v_i`$ coordinate, make the omitted central interval longer than $`\kappa_i`$. Thus each observed run in (14) lies in one periodic end and sees a nonzero value. All other observations outside the core, even if they lie in the opposite end or in a different convex piece, obey the same transverse support bounds (12). It follows that

```math
\min Y_i(u)=\mu_i^--\pi_i(u),\qquad
\max Y_i(u)=\mu_i^+-\pi_i(u)
\qquad\text{(17)}
```

in every good direction.

Suppose $`u'=u+h`$ has the same pattern. A nonzero observation at the first height in (17), together with the global strip support bound for $`f_i`$, implies

```math
\alpha_i-\mu_i^-\le\pi_i(h)\le\gamma_i-\mu_i^-.
\qquad\text{(18)}
```

Two independent good directions bound $`\|h\|`$ by a fixed constant $`K`$. Taking the maximum over all direction pairs makes $`K`$ uniform. This step does not assume the competing anchor to be good or distant.

Choose $`R>R_0+K`$, also accommodating the preceding run condition. Because $`0\in E_i`$, every observation center $`u+s`$ lies in $`u+nP`$; all competing centers $`u'+s`$ then also avoid the cores. Their nonzero heights lie in

```math
[\mu_i^--\pi_i(u'),\ \mu_i^+-\pi_i(u')].
```

Pattern equality makes $`Y_i(u')=Y_i(u)`$. This last interval has exactly the distance between the two extrema in (17), so containing both forces its endpoints to coincide with them. Thus $`\pi_i(h)=0`$ in two independent directions and $`h=0`$.

The entire window need not lie at a single end of a strip. Only the short complete-period runs need this property; the other observations use the common end extrema. This is the point that replaces the convex-window argument. $`\square`$

**Lemma 4.2 (central finite marker).** There is a fixed finite set $`F_0\subset\mathbb Z^2`$, containing $`0`$, such that every anchor satisfying $`F_0\subset u+S_n`$ has a globally unique pattern.

**Proof.** Choose $`i\ne j`$, nonzero points $`z_i,z_j`$ of $`f_i,f_j`$, and $`d=z_j-z_i`$. Then

```math
M(z)=f_i(z)f_j(z+d)
```

is nonzero at $`z_i`$ and has finite support $`F`$: its support lies in the bounded intersection of two nonparallel strips. It is a local function of $`\theta`$ with template $`W=E_i\cup(d+E_j)`$. Set $`F_0=\{0\}\cup(F+W)`$.

If $`F_0\subset u+S_n`$, the pattern observes $`M(z)`$ at every $`z\in F`$. An equal pattern at $`u'=u+h`$ gives $`M(z+h)=M(z)\ne0`$ there, so $`F+h\subseteq F`$. A finite nonempty set cannot contain a nonzero translate of itself: equal cardinality would give equality, and iterating a nonzero translation would give infinitely many points. Hence $`h=0`$. $`\square`$

This finite strip-product marker is the mechanism already used in CN Lemma 4.2 and the [convex area proof](all-star-extension.md); it is not claimed as a new construction.

**Lower bound in Theorem 1.** Enlarge $`R`$ so that $`F_0\subset B_R(0)`$. Let

```math
\mathcal O_n=(-nP)+B_R(0),\qquad
\mathcal I_n=\{u\in\mathbb Z^2:F_0\subset u+S_n\}.
```

Every anchor in $`\mathcal G_n\setminus\mathcal O_n`$ is globally unique by Lemma 4.1, and every anchor in $`\mathcal I_n`$ is globally unique by Lemma 4.2. If a point of $`-nP`$ lies farther than $`R`$ from its boundary, it belongs to $`\mathcal I_n`$. If a point lies in $`\mathcal O_n`$ but outside $`-nP`$, it lies within $`R`$ of that boundary. Therefore

```math
|((\mathcal O_n\cap\mathbb Z^2)\setminus\mathcal I_n)|=O(n).
\qquad\text{(19)}
```

This follows from the boundary estimate after Lemma 2.2, also for nonconvex $`P`$ with holes. The integral nature of $`F_0`$ ensures that containment in $`u+nP`$ here is exactly the required containment in $`u+S_n`$.

All anchors in the union of the two uniqueness sets give different patterns, including between the two sets, since uniqueness is global. From (16) and (19),

```math
P_\theta(S_n)\ge|\mathcal G_n|-O(n)=\Lambda(P)n^2-O(n).
\qquad\text{(20)}
```

## 5. The matching upper bound

For an anchor $`u`$, call direction $`i`$ active when

```math
\pi_i(u)\in J_i^+(n):=\bigcup_a
[\ell_i-nb_{ia},\ r_i-na_{ia}].
\qquad\text{(21)}
```

This definition may overcount actual transitions. If direction $`i`$ is inactive, each individual piece $`u+nK_a`$ lies completely in one of its two tail half-planes. Different pieces can lie in different tails, but there are at most $`2^t`$ choices for this component. The interval definition also works for an empty transition strip; a piece spanning the jump is then included among the active cases.

If at least two directions are active, the number of anchors, and hence the number of their patterns, is at most $`\Lambda(P)n^2+O(n)`$ by Lemma 2.2. The interval lists in (21) have exactly the unshifted projection sets used in (1).

If just one direction $`i`$ is active, fix its transverse height $`\pi_i(u)`$ and the tail choice of every other component on every piece. There are $`O(n)`$ possible active heights and at most $`2^{t(m-1)}`$ such tail assignments. At a fixed height any two integer anchors differ by an integer multiple of primitive $`v_i`$. After also fixing that multiple modulo $`\kappa_i`$, their patterns agree: $`H_i`$ preserves $`F_i`$ and every selected tail. Thus this case contributes only $`O(n)`$ patterns. This argument concerns patterns, not the possibly infinite number of anchors in a one-direction region.

If no direction is active, fix all piecewise tail choices and the anchor modulo the common finite-index tail period lattice. These data determine the pattern, giving $`O(1)`$ possibilities. Overcounting inconsistent or unrealized choices is harmless. The fact that some pieces overlap is likewise harmless: we only use these choices to bound actual patterns from above.

Combining the three cases gives

```math
P_\theta(S_n)\le\Lambda(P)n^2+O(n).
```

Together with (20), this proves Theorem 1. The proof uses neither CN's convex Nivat theorem nor a conjectured nonconvex periodicity statement. $`\square`$

## 6. Equality rigidity and geometric consequences

For every pair $`i\ne j`$,

```math
P\subseteq B_i\cap B_j\subseteq\mathcal M(P).
\qquad\text{(22)}
```

This proves $`\Lambda(P)\ge\mathrm{area}(P)`$.

The sets $`P`$, each $`B_i\cap B_j`$, and $`\mathcal M(P)`$ are closures of their interiors. If two such sets $`P\subseteq Q`$ have equal area and $`P`$ is closed, they are equal: any point of $`Q\setminus P`$ has a neighborhood disjoint from $`P`$ which meets $`\mathrm{int}Q`$, producing a positive-area difference. Thus equality in the area bound forces

```math
P=B_i\cap B_j\quad\hbox{for every pair }i\ne j.
\qquad\text{(23)}
```

For $`m=2`$, this is exactly (3), and conversely (3) gives equality.

Suppose there are three directions and equality holds. Take coordinates

```math
X=\pi_1(z),\quad Y=\pi_2(z),\quad Z=\pi_3(z)=aX+bY.
```

Nonparallelness implies $`a,b\ne0`$. Equip $`P`$ with normalized area measure. By (23) for the pair $`(1,2)`$, its image is $`A_1\times A_2`$ with uniform product measure: the Jacobian is constant. Hence $`X,Y`$ are independent, and $`\mathrm{Var}(X)>0`$. By (23) for $`(1,3)`$, $`X,Z`$ are also independent. However,

```math
0=\mathrm{Cov}(X,Z)
=a\mathrm{Var}(X)+b\mathrm{Cov}(X,Y)
=a\mathrm{Var}(X)\ne0,
```

a contradiction. All variables are bounded, so these elementary second-moment calculations are valid. This rules out equality whenever $`m\ge3`$ and proves Theorem 2.

If $`P`$ is connected, each continuous image $`\pi_i(P)`$ is a compact interval. It has the same extrema as $`\pi_i(\mathrm{Conv}(P))`$, proving (4). In the two-direction equality case both interval unions in (3) are therefore single intervals, giving precisely the asserted parallelogram. Conversely every such parallelogram satisfies (3). Equation (5) follows from (2) and (10). $`\square`$

Here are exact geometric examples with horizontal and vertical component directions:

| Region $`P`$ | Area | $`\Lambda(P)`$ |
|---|---:|---:|
| $`([0,3]\times[0,1])\cup([0,1]\times[0,3])`$ | 5 | 9 |
| $`[0,3]^2\setminus(1,2)^2`$ | 8 | 9 |
| $`[0,1]^2\cup[3,4]^2`$ | 2 | 4 |
| $`([0,1]\cup[3,4])^2`$ | 4 | 4 |

The second example has an area-order hole; its coefficient is unchanged from the filled square. For the third example the convex hull has coefficient $`16`$, whereas the actual region has coefficient $`4`$. More generally, replacing $`3`$ by $`L>1`$ gives the same actual coefficient $`4`$ and convex-hull coefficient $`(L+1)^2`$. Thus convexifying disconnected windows can overestimate the coefficient by an arbitrarily large factor.

Adding any fixed positive-width connecting polygonal corridor to two separated squares makes the region connected and restores the convex-hull projections. The width must be fixed before $`n\to\infty`$. A zero-width corridor is a different problem, as §8 demonstrates.

For rational polygonal regions, the [definable-pattern construction](eventual-quasipolynomial.md#2-finite-descriptions-of-all-star-colors) already implies eventual quasipolynomiality: replace its convex window predicate by the finite Boolean combination defining the region. Combining that existing theorem with (2) identifies its quadratic coefficient as $`\Lambda(P)`$. This is a corollary of the earlier counting result; no new quasipolynomial or period theorem is claimed here.

### 6.1. Closing the product exception for two periodic lines

Assume now that $`m=2`$, that there is no additional background, and that $`F_i`$ is zero off the affine lattice line $`\pi_i=a_i`$, where $`a_i\in\mathbb Z`$. On that line it is a nonzero periodic word, of minimal period $`\kappa_i`$ in primitive $`v_i`$ steps. Individual entries of this word may be zero. The two lines can meet at a lattice point, where their colors are added in $`\mathbb F_p`$ and may cancel.

Let $`P`$ be a product region (3), and put

```math
U_i(n)=nA_i\cap\mathbb Z,\quad N_i(n)=|U_i(n)|,\quad
T(z)=(\pi_1(z),\pi_2(z)),\quad \Gamma=T(\mathbb Z^2),\quad a=(a_1,a_2).
```

**Theorem 4 (exact two-line product count).** For all sufficiently large integers $`n`$,

```math
\boxed{P_\theta(nP\cap\mathbb Z^2)
=C_a(n)+\kappa_1N_1(n)+\kappa_2N_2(n)+1,}
\qquad\text{(24)}
```

where

```math
C_a(n)=|(U_1(n)\times U_2(n))\cap(a+\Gamma)|.
```

Moreover,

```math
C_a(n)=|nP\cap\mathbb Z^2|+O(1),
\qquad\text{(25)}
```

so this entire subclass, including disconnected product regions, has strictly more patterns than sites for all sufficiently large $`n`$.

**Proof.** Set $`d=|\det(v_1,v_2)|`$. A section at any height $`t\in U_i(n)`$ contains lattice points whose other transverse coordinate runs in an arithmetic progression with step $`d`$. Every nondegenerate interval of $`nA_j`$ has length tending to infinity. Consequently, uniformly over all such $`t`$, each section contains an arbitrarily long run of consecutive points in primitive direction $`v_i`$. In particular the actual projection of the lattice window is $`U_i(n)`$ for large $`n`$.

Choose $`n`$ large enough that every section contains at least four nonzero entries of every phase of the corresponding periodic word, and enough consecutive points to distinguish all its $`\kappa_i`$ phases. These are uniform requirements: both words and the interval lists are fixed.

Whenever an affine support line meets the window, at least three of its nonzero sites remain nonzero in the sum, since the other support line can affect at most one of them. The pattern therefore identifies its transverse position. Indeed, a competing pattern with a different parallel support line could cover at most one of these three sites using its other direction. This proves identification even when the intersection color cancels or the two words use the same colors.

If both lines are present, their two transverse positions determine the anchor uniquely. Writing those local positions as $`(t_1,t_2)\in U_1(n)\times U_2(n)`$, the anchor equation is $`T(u)=a-(t_1,t_2)`$. It is solvable in $`\mathbb Z^2`$ precisely when $`(t_1,t_2)\in a+\Gamma`$, and then the anchor is unique. Thus there are exactly $`C_a(n)`$ two-line patterns, all globally unique.

If exactly line $`i`$ is present at height $`t`$, anchors vary by integer multiples of $`v_i`$. Each of the $`\kappa_i`$ phases is realizable while avoiding the other support line: the forbidden translates form a finite set, whereas every residue class of multiples is infinite. The long section distinguishes these phases. Different heights or different visible directions give different patterns by the same support identification. Thus the one-line counts are exactly $`\kappa_iN_i(n)`$. With neither line there is the single zero pattern, which is realizable by avoiding the finitely many anchor lines. This proves (24).

For (25), choose $`w`$ with $`\det(v_1,w)=1`$ and write $`v_2=bv_1+\varepsilon d w`$, where $`\varepsilon\in\{1,-1\}`$ and $`\gcd(b,d)=1`$. In these coordinates

```math
\Gamma=\{(t,s):s\equiv bt\pmod d\}.
```

If $`d=1`$, every coset count is identical. For general $`d`$, let $`n_i(r)`$ count members of $`U_i(n)`$ in residue $`r`$ modulo $`d`$. Since each $`U_i(n)`$ is a fixed finite number of integer intervals,

```math
n_i(r)=N_i(n)/d+\epsilon_i(r),\qquad
\sum_r\epsilon_i(r)=0,\qquad |\epsilon_i(r)|\le C_i,
```

with constants independent of $`n`$. Every coset of $`\Gamma`$ is given by $`s\equiv bt+c\pmod d`$, and its count is

```math
\sum_r n_1(r)n_2(br+c)
=\frac{N_1(n)N_2(n)}d+\sum_r\epsilon_1(r)\epsilon_2(br+c).
```

The linear cross terms vanish because $`b`$ permutes residues. All coset counts therefore differ by $`O(1)`$. The zero-coset count equals the number of sites, proving (25). Finally $`N_i(n)=n\,\mathrm{length}(A_i)+O(1)`$ with positive lengths, so the remaining terms in (24) eventually dominate the bounded difference. $`\square`$

The proof covers sparse periodic words directly. It does not infer their behavior from the earlier full-support convex-line formula, and it does not extend (24) to finite-width transition strips, unequal tails or arbitrary periodic backgrounds.

**Corollary 4.1.** Every two-line star in Theorem 4 satisfies $`P_\theta(nP\cap\mathbb Z^2)>|nP\cap\mathbb Z^2|`$ eventually for every region $`P`$ of §1, whether connected or disconnected. For a nonproduct region use the strict quadratic surplus of Theorem 2; for a product region use (24)–(25). This corollary needs no input from CN Theorem T.

## 7. Proof of the exact cross-star formula

For a translate of (6), the set of ones in the $`S`$-pattern is the symmetric difference of one row section and one column section. Either section may be empty. Every pair of choices is realized independently by an integer anchor: choose the corresponding negative coordinate for an active section and any coordinate outside the finite projection set for an inactive section.

Write $`*`$ for the inactive choice. There are $`(r+1)(c+1)`$ states

```math
(a,b)\in(R\cup\{*\})\times(C\cup\{*\}).
```

Associate to a state the subset of graph vertices consisting of its selected row and selected column. Its pattern is the cut of this subset, read modulo two on the edges of $`G_S`$. Two states have the same pattern if and only if the symmetric difference of their vertex subsets is a union of entire connected components. To verify this, the cut of that difference vanishes precisely when membership is the same at both ends of every edge, which is equivalent to constancy on each connected component.

The difference between two states uses at most two row vertices and at most two column vertices. There are no isolated graph vertices. Consequently a nonempty component union that can be such a difference is exactly one of:

- one component of type $`(1,1)`$, $`(1,2)`$, $`(2,1)`$ or $`(2,2)`$;
- two components, both of type $`(1,1)`$.

This is a complete list, since every component uses at least one vertex of each kind.

Let the isolated-edge components be $`(a_i,b_i)`$, $`1\le i\le e`$. The resulting nontrivial equivalence classes are:

1. The state $`(*,*)`$ and the $`e`$ states $`(a_i,b_i)`$ all give the zero pattern. This reduces the state count by $`e`$.
2. For each $`i`$, $`(a_i,*)`$ and $`(*,b_i)`$ coincide. This reduces it by another $`e`$.
3. For each unordered pair $`i\ne j`$, $`(a_i,b_j)`$ and $`(a_j,b_i)`$ coincide. This reduces it by $`\binom e2`$.
4. A component of type $`(1,2)`$, with row $`a`$ and columns $`b_1,b_2`$, gives the two classes $`\{(*,b_1),(a,b_2)\}`$ and $`\{(*,b_2),(a,b_1)\}`$. The transposed type behaves identically. Each reduces the count by two.
5. A component of type $`(2,2)`$ gives the two classes $`\{(a_1,b_1),(a_2,b_2)\}`$ and $`\{(a_1,b_2),(a_2,b_1)\}`$, again reducing the count by two.

These classes are disjoint. The preceding exhaustive list of possible component differences shows that there are no additional identifications or larger classes. Subtracting the reductions proves (7). $`\square`$

For example, when $`G_S`$ has no component with at most two rows and at most two columns, the state map is injective and $`P_{\theta_\times}(S)=(r+1)(c+1)`$ exactly. For $`e`$ isolated diagonal sites, (7) gives $`1+e+\binom e2`$: the patterns are exactly all subsets of the sites of size zero, one or two.

Applied to the four regions in §6, for every integer $`n\ge2`$ this yields

```math
P_{\theta_\times}(nP\cap\mathbb Z^2)=
\begin{cases}
(3n+2)^2,&P\text{ is the L-region or square frame},\\ 
(2n+3)^2,&P\text{ is either of the two disconnected examples}.
\end{cases}
\qquad\text{(26)}
```

The cross star is only a concrete test case for Theorem 1; (2) covers all stars satisfying the stated hypotheses.

## 8. Zero-width appendages change the leading coefficient

Fix an integer $`L\ge2`$ and let

```math
P_L=[0,1]^2\cup\{(t,t):1\le t\le L\}.
```

This compact connected set has area one, but is not a region in the sense of §1. Its dilated lattice window is

```math
S_n=\{0,\ldots,n\}^2\cup\{(j,j):n\lt j\le Ln\}.
```

For $`n\ge2`$, its incidence graph is one component $`K_{n+1,n+1}`$ and $`e=(L-1)n`$ isolated edges. Thus $`r=c=Ln+1`$ and $`t=0`$. Theorem 3 gives the exact identity

```math
\boxed{P_{\theta_\times}(S_n)
=\frac{L^2+2L-1}{2}n^2+\frac{5L+3}{2}n+4.}
\qquad\text{(27)}
```

The projection sets of $`P_L`$ are full intervals of length $`L`$, so the formula (1), if incorrectly applied outside Theorem 1, would predict $`L^2n^2`$. Its quadratic overestimate is $`(L-1)^2n^2/2`$. The correct coefficient in (27) follows from genuine collisions of anchors, not from a lattice rounding error.

Only $`(L-1)n`$ sites were added to the original square window. Nevertheless its complexity coefficient changed from $`1`$ to $`(L^2+2L-1)/2`$. Thus an $`O(n)`$ change in the number of observed sites does not imply an $`O(n)`$ change in pattern complexity. The boundary estimates in §§2–4 count exceptional anchors and cannot be replaced by such an argument about sites.

This counterexample does not contradict Theorem 1. Its appendage has zero width and supplies no complete-period runs inside an eroded two-dimensional piece. It also does not give a low-complexity nonperiodic counterexample: (27) remains larger than $`|S_n|`$. Its role is to make the geometric hypothesis and the anchor-versus-site distinction precise.

## 9. What remains open here

For two directions and a disconnected product region (3), Theorem 1 proves the exact coefficient but leaves an $`O(n)`$ difference from the number of sites. Theorem 4 proves eventual positivity for every two-line star in its stated scope. For general finite-width components, unequal tails or periodic backgrounds this sign is still unresolved here. The full convex finite-window theorem already handles the connected parallelogram case.

The specific remaining interface is a collision estimate for partially visible transverse transition profiles at the several interval ends. The marker/complete-extrema proof discards a constant-width family of anchor strips, containing $`O(n)`$ anchors; the one-direction modes can also contribute $`O(n)`$ patterns. At an equality region these are the very terms whose difference must be decided. Their cardinalities cannot be compared by the area argument. A separate finite search of 240 axis stars at three sizes found no low-complexity example; that observation supplies no sign theorem.

Periodic sparse masks and arbitrary finite nonconvex windows are outside the theorem. The statements here require precisely the full-dimensional finite-union geometry of §1.

Mixed-dimensional sets such as $`P_L`$ need a collision-sensitive invariant, rather than just their projection sets. Theorem 3 resolves this completely for the cross star and every finite window, but no such general formula is proved here for all stars.

## 10. Evidence, attribution and dependencies

The companion [nonconvex-window diagnostic](../scripts/nonconvex_window_audit.py) uses exact integer/Fraction arithmetic and the Python standard library. It performs:

- direct full-plane pattern enumeration for all 65,535 nonempty subsets of a $`4\times4`$ grid, checking (7) against a separately computed component formula;
- 39 exact instances of (26)–(27);
- exact projection-union area calculations for the four regions, in two and three directions;
- 42 further pattern counts on nonconvex windows for unequal-tail stars, periodic line words, step stars, and three full lines, with the quadratic residuals recorded.

For the cross star, one representative for each active row/column or inactive choice exhausts all anchors. For the explicit axis stars, enumeration of the finite transition envelope and representatives of every tail parity exhausts the plane. For full constant lines, all anchors with two visible lines are enumerated by solving the two transverse equations; zero- and single-line modes are included separately. These are exhaustive finite classifications for the specific examples, not arbitrary sampled boxes. The general asymptotic result is proved in §§2–5 and does not follow from the numerical residuals. The complete deterministic output is [packaged here](../scripts/nonconvex_window_audit.json).

The [two-line product diagnostic](../scripts/two_line_product_audit.py) additionally checks eight exact instances of Theorem 4, with pair determinants $`1`$, $`2`$ and $`3`$, both determinant signs, zero entries in periodic words, lattice and nonlattice line intersections, and cancellation. All eight agree with (24). Their nonzero coset-count differences are $`-4`$, making it necessary to retain $`C_a(n)`$ rather than replace it by the site count in the exact formula. Its [deterministic output](../scripts/two_line_product_audit.json) retains every mathematical case from the reviewed source.

A separate exploratory screen, not bundled, checked 240 generated two-axis stars at three window sizes and checked its normalization against literal pixel patterns in 12 cases. It found no counterexample for the unresolved general product case; that observation does not establish an eventual sign theorem.

CN §0.1 supplies the star definition; Lemmas 4.1–4.2 supply antecedents of the signal and marker mechanisms. The proofs of Theorems 1–4 do not invoke CN Theorem T. Corollary 2.1 invokes it only for the aligned parallelogram case. Exact source versions and printed pages are in the [bibliography](../sources/bibliography.md).

The [convex area theorem](leading-term.md), [all-star extension](all-star-extension.md), [finite-window bound](finite-window-bound.md), [eventual counting theorem](eventual-quasipolynomial.md), and [full-support line formula](full-support-lines.md) delimit the preceding results. The new scope here is the nonconvex projection-union coefficient, the nonconvex equality classification, and the stated exact collision formulas. Convex-window equality, the strip/marker mechanisms, and general definable-pattern counting are not claimed anew. The general [coefficient](explicit-coefficients.md) and [minimal-period](minimal-periods-and-constant-term.md) constructions are separate results.

The manuscript was developed by GPT-6 Astra and established in one fresh-context GPT-6.1 Sol adversarial review. No comprehensive literature-priority search, human certification or proof-assistant verification is claimed. Finite diagnostics do not replace the mathematical proofs.

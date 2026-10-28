# Integral nonconvex windows can have arbitrary linear-coefficient period

**Review status: established in one fresh-context AI adversarial review.** The source manuscript was reviewed on 28 October 2026, with zero fatal, proof-repair or presentation issues. Source and review fingerprints are in the [review record](../sources/review-status.md#integral-nonconvex-windows-established-in-one-fresh-context-ai-adversarial-review).

The linear coefficient need not be constant on an integral nonconvex window. A pair of separated unit squares gives minimal period two. For every integer $`k\ge2`$, a binary two-line star and an integral concave pentagon with connected interior give minimal period $`k`$. The same period is realized by a region with one genuine hole.

The proofs count all plane translates exactly for every $`n\ge2`$. They do not infer stabilization from computation and do not require the general quasipolynomial theorem, coefficient construction, or nonconvex leading-term theorem. The comparison with those results is in Section 8. The AI review is not human certification or proof-assistant verification. No worldwide priority claim is made.

## 1. Definitions and results

For a configuration $`\theta:\mathbb Z^2\to\mathbb F_2`$ and a finite window $`S\subset\mathbb Z^2`$, write

```math
P_\theta(S)=\#\{(\theta(u+z))_{z\in S}:u\in\mathbb Z^2\}.
```

An integral polygon has all vertices in $`\mathbb Z^2`$. All regions below are finite unions of closed, positive-area integral convex polygons. Dilation is about the origin, $`S_n=nP\cap\mathbb Z^2`$, with positive integer $`n`$.

For an exact eventual formula

```math
P_\theta(S_n)=\Lambda n^2+\beta_{n\bmod M}n+\gamma_{n\bmod M},
```

$`p_\beta`$ denotes the least period of the linear-coefficient sequence. The coefficient sequences are intrinsic: two eventual quasipolynomial presentations agree as polynomials on every infinite progression in their common refinement. In particular, a varying linear coefficient cannot be absorbed into a bounded constant term.

For integers $`k\ge2`$ and $`s`$, use the single explicit family

```math
\theta_{k,s}(x,y)
=\mathbf1[x=0]+\mathbf1[y=sx,\ x\equiv0\pmod k]
\quad\text{in }\mathbb F_2.
\qquad\text{(1)}
```

Thus the origin has value zero, because its two contributions cancel. Every proof below keeps this cancellation.

**Theorem 1 (two separated squares).** Let

```math
P_{\mathrm{sq}}=[0,1]^2\ \cup\ ([0,1]\times[2,3]),
\qquad \theta=\theta_{2,1}.
```

For every integer $`n\ge2`$,

```math
\boxed{
P_\theta(nP_{\mathrm{sq}}\cap\mathbb Z^2)
=\begin{cases}
4n^2+\dfrac{25}{2}n+2,&n\equiv0\pmod2,\\ 
4n^2+13n+2,&n\equiv1\pmod2.
\end{cases}}
\qquad\text{(2)}
```

Consequently $`p_\beta=p_f=2`$ and $`p_\gamma=1`$. Thus constant $`\beta`$ fails for finite unions of positive-area integral convex polygons.

**Theorem 2 (every period, connected interior, no hole).** Fix any integer $`k\ge2`$ and put $`L=k+1`$. Let

```math
\begin{split}
P_k={}&([0,L]\times[-1,0])\\ 
&\cup\mathrm{conv}\{(0,0),(1,0),(0,1)\}\\ 
&\cup\mathrm{conv}\{(1,0),(L,0),(L,1)\}.
\end{split}
\qquad\text{(3)}
```

This is the closed concave pentagon with boundary vertices, in order,

```math
(0,-1),\ (L,-1),\ (L,1),\ (1,0),\ (0,1).
```

It has connected interior and no hole. Define

```math
\rho=n\bmod k\in\{0,\ldots,k-1\},\qquad
A_\rho=\begin{cases}1,&\rho=0,\\2,&\rho\ne0.\end{cases}
```

For $`\theta=\theta_{k,0}`$ and every integer $`n\ge2`$,

```math
\boxed{
P_\theta(nP_k\cap\mathbb Z^2)
=2Ln^2+3Ln+2+A_\rho\left(\left\lfloor\frac{Ln}{k}\right\rfloor+2\right).
}
\qquad\text{(4)}
```

Equivalently,

```math
\boxed{
\Lambda=2L,\qquad
\beta_\rho=3L+\frac{L}{k}A_\rho,\qquad
\gamma_\rho=2+A_\rho\left(2-\frac{\rho}{k}\right).
}
\qquad\text{(5)}
```

All three minimal periods satisfy $`p_\beta=p_\gamma=p_f=k`$.

In particular the fixed binary star

```math
\theta_{2,0}(x,y)=\mathbf1[x=0]+\mathbf1[y=0,\ x\equiv0\pmod2]
```

on the pentagon with vertices $`(0,-1),(3,-1),(3,1),(1,0),(0,1)`$ gives

```math
P_{\theta_{2,0}}(nP_2\cap\mathbb Z^2)
=\begin{cases}
6n^2+\dfrac{21}{2}n+4,&n\text{ even},\\ 
6n^2+12n+5,&n\text{ odd},
\end{cases}
\qquad n\ge2.
\qquad\text{(6)}
```

**Theorem 3 (every period, one hole).** Retain the two triangles in (3), and replace its base rectangle by

```math
\begin{split}
B_k={}&([0,1]\times[-3,0])\cup([2,L]\times[-3,0])\\ 
&\cup([1,2]\times[-3,-2])\cup([1,2]\times[-1,0]).
\end{split}
```

Write $`Q_k=B_k\cup T_{\mathrm{left}}\cup T_{\mathrm{right}}`$. This is the closed polygonal region bounded externally by

```math
(0,-3),\ (L,-3),\ (L,1),\ (1,0),\ (0,1),
```

with the open rectangle $`(1,2)\times(-2,-1)`$ removed. The hole lies strictly inside the outer polygon. The region has connected interior, and its displayed six pieces are positive-area integral convex polygons. For $`\theta=\theta_{k,0}`$ and every $`n\ge2`$,

```math
\boxed{
P_\theta(nQ_k\cap\mathbb Z^2)
=4Ln^2+5Ln+2+A_\rho\left(\left\lfloor\frac{Ln}{k}\right\rfloor+2\right).
}
\qquad\text{(7)}
```

Here $`\beta_\rho=5L+(L/k)A_\rho`$, the constant table is the same as in (5), and $`p_\beta=p_\gamma=p_f=k`$.

## 2. The configurations are CN stars

The two summands in (1) have primitive directions

```math
v_{\mathrm{vert}}=(0,1),\qquad v_{\mathrm{sparse}}=(1,s).
```

They are nonparallel and their pair determinant has absolute value one. The vertical summand has period $`(0,1)`$. The sparse summand has period $`k(1,s)`$, and its minimal positive period in primitive tangential steps is exactly $`k`$. Each summand is zero off its own supporting line, so it has zero doubly periodic tails on the two sides of a width-zero transition strip. Neither summand is doubly periodic: a period must preserve its supporting line, and hence must be parallel to that line. They therefore satisfy the star hypotheses in CN §0.1, with two components and prime field $`\mathbb F_2`$.

The complete configuration also has no nonzero period. Indeed, its support has a cofinite subset of the vertical line $`x=0`$, whereas on any other vertical line it has at most one point. A period $`(a,b)`$ must therefore have $`a=0`$. Its infinitely many off-vertical support points lie on $`y=sx`$; invariance under $`(0,b)`$ then forces $`b=0`$.

Thus the examples do not rely on a periodic configuration, dependent directions, nonintegral vertices, a higher-index pair lattice, unequal tails, a nonzero background, or a zero-width geometric appendage. The zero-width supports of the *components* are allowed CN stars; the *windows* are full-dimensional polygonal regions.

## 3. A complete finite-window identity

The following elementary lemma is the counting tool. Its proof is also the completeness proof for the translation representatives used by the verification scripts.

**Lemma 4 (a full column and a periodically sparse line).** Fix $`k\ge2`$, an integer $`s`$, and a nonempty finite $`S\subset\mathbb Z^2`$. Let

```math
X=\{x:(x,y)\in S\text{ for some }y\},\qquad c=|X|,
```

and suppose **every nonempty column of $`S`$ contains at least three sites**. For $`r\in\mathbb Z/k\mathbb Z`$, put

```math
\begin{split}
c_r&=\#\{x\in X:x\equiv r\pmod k\},\\ 
H_r&=\#\{h\in\mathbb Z:\text{there is }(x,y)\in S,
\ y-sx=h,\ x\equiv r\pmod k\}.
\end{split}
```

Then

```math
\boxed{P_{\theta_{k,s}}(S)=c+1+\sum_{r=0}^{k-1}(c_r+1)H_r.}
\qquad\text{(8)}
```

**Proof.** Given any anchor $`u=(U,V)\in\mathbb Z^2`$, define

```math
a=-U,\qquad h=sU-V.
```

These are independent integer parameters: any $`(a,h)\in\mathbb Z^2`$ is realized by $`u=(-a,-sa-h)`$. Identify a binary pattern with its set of ones. Set

```math
C_a=\{(x,y)\in S:x=a\},\qquad
E_{h,r}=\{(x,y)\in S:y-sx=h,\ x\equiv r\pmod k\}.
```

The pattern of this anchor is exactly

```math
C_a\mathbin{\triangle}E_{h,a\bmod k},
\qquad\text{(9)}
```

where $`\triangle`$ is symmetric difference. This formula accounts for the cancellation at an intersection.

Suppose first $`a\in X`$. At most one site of $`C_a`$ lies on $`y-sx=h`$. Consequently the pattern has at least two ones on column $`a`$. Every other column has at most one one. Thus the pattern intrinsically identifies $`a`$; it cannot coincide with a pattern with a different active column or with no active column.

For this fixed $`a`$, put $`r=a\bmod k`$. The nonempty sets $`E_{h,r}`$, as $`h`$ varies, are mutually disjoint and therefore distinct. Symmetric difference with the fixed $`C_a`$ is injective on sets. There are $`H_r`$ nonempty choices and one empty choice, all realizable by the independent parameter $`h`$. Hence this column contributes $`H_r+1`$ distinct patterns. The total with an active vertical line is

```math
\sum_r c_r(H_r+1)=c+\sum_r c_rH_r.
```

If $`a\notin X`$, the vertical contribution vanishes. For each residue $`r`$, there are infinitely many such $`a\equiv r\pmod k`$; thus every $`E_{h,r}`$ is realizable. The nonempty sets indexed by different pairs $`(h,r)`$ are disjoint, so all are distinct, and the empty set is one additional pattern. These patterns have at most one one per column and cannot coincide with any of the preceding patterns. Their number is $`1+\sum_r H_r`$. Adding the two counts proves (8). Every integer anchor belongs to one of these cases. $`\square`$

No lower bound on the number of sparse-line sites is needed. In particular, an empty phase, a singleton, two separated sites on one line, and a canceled intersection are all included in the same exact count.

## 4. Proof of Theorem 1

The dilated two-square window is

```math
S_n=\{0,\ldots,n\}\times
\bigl(\{0,\ldots,n\}\cup\{2n,\ldots,3n\}\bigr).
```

Each column has $`2n+2\ge3`$ sites. Apply Lemma 4 with $`k=2`$, $`s=1`$, and $`c=n+1`$.

For parity $`r\in\{0,1\}`$, let $`a_r`$ and $`b_r`$ be the smallest and largest integers of that parity in $`[0,n]`$. The possible heights $`h=y-x`$ in the lower square, restricted to this parity of $`x`$, form exactly the integer interval

```math
I_r^-=[-b_r,n-a_r]\cap\mathbb Z.
```

Indeed, the intervals $`[-x,n-x]\cap\mathbb Z`$ from consecutive allowed columns are shifted by two and have no integer gaps for $`n\ge2`$. The upper square similarly supplies

```math
I_r^+=[2n-b_r,3n-a_r]\cap\mathbb Z.
```

Thus $`H_r=|I_r^-\cup I_r^+|`$.

If $`n`$ is even, then

```math
\begin{array}{c|cccc}
r&a_r&b_r&c_r&H_r\\ \hline
0&0&n&n/2+1&4n+1\\ 
1&1&n-1&n/2&4n-2
\end{array}
```

For $`r=0`$, the two intervals meet at the single height $`n`$; each has length $`2n+1`$ in number of sites. For $`r=1`$, each has $`2n-1`$ sites and they are disjoint, with the height $`n`$ missing. Substitution in (8) gives

```math
\begin{split}
P_\theta(S_n)
&=n+2+\left(\frac n2+2\right)(4n+1)
       +\left(\frac n2+1\right)(4n-2)\\ 
&=4n^2+\frac{25}{2}n+2.
\end{split}
```

If $`n`$ is odd, both parity classes have $`c_r=(n+1)/2`$. For each parity the lower and upper intervals are adjacent, disjoint intervals, each of cardinality $`2n`$. Hence $`H_0=H_1=4n`$, and

```math
P_\theta(S_n)=n+2+2\left(\frac{n+3}{2}\right)4n
=4n^2+13n+2.
```

This proves (2) at every stated scale. Its two different linear coefficients prove minimal period two. A period-one polynomial would have quadratic coefficient four and would require both $`25/2`$ and $`13`$ as its linear coefficient, by polynomial identity on the even and odd infinite progressions. The constant table is identically two. $`\square`$

## 5. Proof of Theorem 2

The geometric claims follow directly from (3). For $`0\le x\le L`$, the vertical fiber is the interval $`[-1,g(x)]`$, where

```math
g(x)=\begin{cases}
1-x,&0\le x\le1,\\ 
(x-1)/k,&1\le x\le L.
\end{cases}
```

It is everywhere of positive length. Every interior point can be connected down its vertical fiber to the interior base rectangle. There is no hole. The top segment joining $`(0,1)`$ to $`(L,1)`$ is not contained in the region, so the polygon is nonconvex. All three pieces in (3) have positive area and integral vertices.

Now fix $`n\ge2`$. Its observed columns are every integer $`x\in[0,Ln]`$, and every column contains $`\{-n,\ldots,0\}`$, so Lemma 4 applies. We have

```math
c=Ln+1,\qquad \sum_r c_r=c.
```

Every row $`-n\le y\le0`$ contains all columns and hence all $`k`$ residues. For $`1\le y\le n-1`$, the lattice section is

```math
\bigl([0,n-y]\cup[n+ky,Ln]\bigr)\cap\mathbb Z.
\qquad\text{(10)}
```

Its right interval has length $`k(n-y)\ge k`$ and thus contains all residues modulo $`k`$. At the top row $`y=n`$, however, the section consists of exactly two sites,

```math
(0,n),\qquad (Ln,n).
\qquad\text{(11)}
```

Their $`x`$-residues are $`0`$ and $`\rho`$, since $`L=k+1`$. With $`s=0`$ in the star, the line height is just $`h=y`$. The $`2n`$ rows from $`-n`$ through $`n-1`$ are seen in every phase, whereas the last row is seen only in the residue set $`\{0,\rho\}`$. Therefore

```math
H_r=2n+\mathbf1[r\in\{0,\rho\}].
\qquad\text{(12)}
```

Put $`q=\lfloor Ln/k\rfloor`$. Because $`Ln=kq+\rho`$, the integer interval $`[0,Ln]`$ contains $`q+1`$ integers of each residue $`0,\ldots,\rho`$, and $`q`$ of each remaining residue. In particular

```math
c_0=c_\rho=q+1
```

when $`\rho\ne0`$; when $`\rho=0`$, there is just the one selected residue, again with count $`q+1`$. In either case

```math
\sum_{r\in\{0,\rho\}}(c_r+1)=A_\rho(q+2).
\qquad\text{(13)}
```

Applying (8), (12), and (13) now gives

```math
\begin{split}
P_\theta(S_n)
&=c+1+2n\sum_r(c_r+1)+A_\rho(q+2)\\ 
&=Ln+2+2n(Ln+1+k)+A_\rho(q+2)\\ 
&=2Ln^2+3Ln+2+A_\rho(q+2).
\end{split}
```

This is (4). Finally, $`q=Ln/k-\rho/k`$, which proves (5).

The $`k`$-entry linear table has its unique smaller value at residue zero:

```math
\beta_0=3L+\frac Lk,\qquad
\beta_\rho=3L+\frac{2L}{k}\quad(1\le\rho\lt k).
```

Thus a positive shift preserving the table must take residue zero to itself; it must be divisible by $`k`$. The table already has period $`k`$, proving $`p_\beta=k`$. Likewise $`\gamma_0=4`$, whereas

```math
\gamma_\rho=6-\frac{2\rho}{k}>4\quad(1\le\rho\lt k),
```

so $`p_\gamma=k`$. Formula (5) supplies a quasipolynomial period $`k`$, and any eventual quasipolynomial period must be a period of its intrinsic linear table, by the common-refinement polynomial identity. Hence $`p_f=k`$. This proves Theorem 2 without fitting or invoking a general period theorem. $`\square`$

## 6. Proof of Theorem 3

The four rectangles making up $`B_k`$ fill $`[0,L]\times[-3,0]`$ except for the open interior rectangle $`(1,2)\times(-2,-1)`$, with all four hole edges retained. The latter rectangle lies strictly inside the base. Its boundary is an actual hole boundary, with positive surrounding widths; the two upper triangles do not meet the hole. This proves the stated geometry and finite integral convex decomposition.

At every scale $`n\ge2`$, the observed columns remain $`0,\ldots,Ln`$. A column outside the open horizontal span of the hole contains all $`3n+1`$ base sites. A column inside it contains at least the two intervals of rows $`[-3n,-2n]`$ and $`[-n,0]`$, of total cardinality $`2n+2`$. Thus every column has at least three sites.

For any row $`-3n\le y\le0`$, if the hole removes some sites, the right interval of columns $`[2n,Ln]\cap\mathbb Z`$ is still present. Its number of sites is $`(k-1)n+1\ge k`$, so it contains every residue modulo $`k`$. In all other base rows the full interval is present. The upper rows have exactly the sections (10) and (11). Consequently

```math
H_r=4n+\mathbf1[r\in\{0,\rho\}].
```

The column counts and (13) have not changed. Substitution in (8) gives

```math
\begin{split}
P_\theta(nQ_k\cap\mathbb Z^2)
&=Ln+2+4n(Ln+1+k)+A_\rho(q+2)\\ 
&=4Ln^2+5Ln+2+A_\rho(q+2),
\end{split}
```

which is (7). Substituting $`q=Ln/k-\rho/k`$ gives its coefficients. The same unique-minimum argument proves all the asserted minimal periods. $`\square`$

## 7. Why integrality no longer removes the phase

The precise obstruction is simultaneous sampling at *different* macroscopic points of a short section, not lattice-point rounding of one vertex.

In a convex polygon a directionally extreme section of bounded size is centered at a single strict extreme vertex. Translating its finite local cone profile by that integral vertex times $`n`$ removes the dilation phase; averaging the remaining tangential phases gives constant linear weights. This is the relevant mechanism in the [integral convex-window theorem](integral-windows.md).

For the connected pentagon, the maximum horizontal section instead consists of the two separated vertices in (11), whose displacement is $`(Ln,0)`$. Centering at one tip leaves the other at tangential distance $`Ln\equiv n\pmod k`$. Thus no common recentering removes their relative phase. There are exactly

```math
|\{0,n\bmod k\}|=A_\rho
```

phases that see a nonzero sparse-line correction on that section. The resulting contribution to complexity is exactly

```math
A_\rho\left(\left\lfloor\frac{Ln}{k}\right\rfloor+2\right).
```

It contains the varying linear term $`(L/k)A_\rho n`$. Each visible phase is paired with order $`n/k`$ distinguishable vertical-column positions, as well as its one no-column mode. If a sparse phase sees neither tip, its top-row representation collapses to the corresponding pure-column pattern; Lemma 4 counts that collapse exactly. Thus a change on a fixed number of endpoint sites can change the *number of patterns* by order $`n`$.

In the two-square example the critical height is $`h=y-x=n`$. Its two sites are $`(0,n)`$ and $`(n,2n)`$. Their tangential separation has parity $`n`$: even scales put them in the same sparse phase, whereas odd scales put them in different phases. The interval-union count in Section 4 is an exact accounting of this effect, including all the other rows.

Both projection sets of $`P_k`$, in the two star directions, are single intervals: the coordinate projections are $`[0,L]`$ and $`[-1,1]`$. Therefore connectedness of the region or connectedness of all directional projections does not prevent this obstruction. The one-hole version shows that a hole does not remove it either. The original example has only two directions, determinant one, and equal zero tails, so none of those simplifications suffices to extend the convex conclusion.

This identifies the obstruction to extending the convex conclusion to all integral nonconvex windows: the convex proof explicitly uses a single strict extreme vertex for a bounded endpoint cluster. Our nonconvex sections have two distant clusters at the same extreme height. No all-star classification of which nonconvex regions avoid such interactions is claimed here.

## 8. Relation to the other results

The exact counts agree with the [nonconvex projection-union leading coefficient](nonconvex-windows.md). There are two directions and a unimodular coordinate pair, so $`\Lambda`$ is the product of the lengths of the two actual transverse projections:

| Region and star | Transverse coordinates | Projection lengths | Exact quadratic coefficient |
| --- | --- | --- | --- |
| $`P_{\mathrm{sq}},\theta_{2,1}`$ | $`x,y-x`$ | $`1,4`$ | $`4`$ |
| $`P_k,\theta_{k,0}`$ | $`x,y`$ | $`L,2`$ | $`2L`$ |
| $`Q_k,\theta_{k,0}`$ | $`x,y`$ | $`L,4`$ | $`4L`$ |

That leading theorem allows an $`O(n)`$ remainder, so the varying coefficients do not contradict it. The disconnected two-line *product-region* formula in [the nonconvex-window theorem](nonconvex-windows.md#61-closing-the-product-exception-for-two-periodic-lines) does not apply to these regions: none is the complete inverse image of the Cartesian product of its two transverse projections. For example the two squares in $`(x,y-x)`$ coordinates do not fill $`[0,1]\times[-1,3]`$. Our use of a long full vertical column does not assert that the sparse-line sections are uniformly long.

The [representative-cell](explicit-coefficients.md) and [exact minimal-period](minimal-periods-and-constant-term.md) constructions can also be applied to rational nonconvex windows: replace the window's facet conjunction in the definable representative formula by the finite disjunction for its convex pieces, and use the nonconvex leading coefficient. This does not change the finite elimination or coefficient-table reduction. Its general allowance for a periodic linear coefficient is essential here. The present proof does not implement the full compiler; it obtains the same relevant information directly from Lemma 4. It is also independent of the convex constant-$`\beta`$ result, which is used only to explain the convexity boundary.

The [integral convex-window proof](integral-windows.md), especially Appendix B, uses recentering at a single strict extreme vertex for bounded endpoint clusters. The [nonconvex-window proof](nonconvex-windows.md) gives the actual projection-union leading term, its equality classification, and the two-line product formula under its product hypothesis. The [object definitions](objects.md), [representative coefficients](explicit-coefficients.md), and [minimal-period construction](minimal-periods-and-constant-term.md) specify the shared star scope and intrinsic coefficient conventions.

The additional results here are the negative answer on integral nonconvex windows, the exact separated-square formula (2), arbitrary prescribed periods on connected pentagons and one-hole regions, and the complete two-line finite-window identity (8). No new general Presburger or Ehrhart principle, comprehensive literature search or worldwide priority claim is asserted.

## 9. Exact diagnostics and their scope

The standard-library programs below are self-contained: they read no external sources, write only adjacent deterministic JSON, and contain no timings or timestamps. Run them with `python3 -B scripts/verify.py`; the runner compares each complete output with its frozen baseline.

- [`nonconvex_beta_search_audit.py`](../scripts/nonconvex_beta_search_audit.py) retains all 810 assertions: complete-plane pattern deduplication, the two-square formula, connected examples for $`2\le k\le12`$, one-hole examples for $`2\le k\le6`$, and 23 literal crosschecks. Its five shape screens include axis-aligned separated squares, a square frame and an L-region; finite constant-linear diagnostics for these controls are not general constancy theorems.
- [`nonconvex_beta_audit.py`](../scripts/nonconvex_beta_audit.py) uses independent integer-inequality geometry and literal pixel evaluation at a proved complete list of anchors. It checks 29 square inputs, 300 connected inputs, 300 one-hole inputs, 180 irregular long-column windows, and exact coefficient tables for $`2\le k\le128`$. The literal path covers 217 inputs.
- [`nonconvex_beta_independent_audit.py`](../scripts/nonconvex_beta_independent_audit.py) imports neither preceding program. Integer winding on ordered boundaries agrees with convex-piece membership on 390 inputs. It directly deduplicates full patterns on 357 family inputs, including $`n\ll k`$; checks 1,860 exhaustive and 400 fixed-seed irregular long-column inputs; recovers coefficients on 110 progressions by exact interpolation with a fourth-point check; and checks 1,023 exact period tables through $`k=1024`$. Its two counting paths compare entire pattern sets on 336 inputs, using 271,522 actual anchors and 68,995,483 pixel evaluations. A negative control with two sites per column has true complexity 7 versus 11 from the incorrectly extended identity, confirming the importance of Lemma 4's threshold. These groups overlap and are not summed as independent samples.

Whole-plane completeness for the first two programs follows from (9): take every observed column parameter $`a`$, one absent $`a`$ in each residue modulo $`k`$, every observed height $`h=y-sx`$, and one absent height. For the third program, the event counter uses the same complete support states. Its literal counter instead uses the entire finite anchor strip

```math
-\max X-k\le U\le-\min X+k,\qquad
sU-\max H-1\le V\le sU-\min H+1,
```

where $`H=\{y-sx:(x,y)\in S\}`$. Every active column and an absent column of every phase are present; each possible intersection height and an empty sparse-line state are present for every $`U`$. An outside anchor can be moved into this strip while preserving its column state, phase and intersection height, or its empty sparse state. Thus both paths exhaust all plane translates, even without the three-site-column assumption.

The exact formulas and all stated minimal periods follow from Sections 2–6, independently of these finite checks. A characterization of all integral nonconvex windows with constant $`\beta`$, or a direct all-star nonconvex boundary formula, remains outside this result. The [review record](../sources/review-status.md) identifies the reviewed source and report; no human expert or proof assistant has certified the work.

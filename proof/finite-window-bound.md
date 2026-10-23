# An explicit finite-window bound for every star

Prerequisites: [objects and statements](objects.md), [extremal recovery](signals-and-anchors.md), and [the all-star strip and marker arguments](all-star-extension.md).

Theorem 2 gives a finite, computable core correction valid for every star. It does **not** assert that the uncorrected expression $`|\mathcal A(S)|`$ from the Case B anchor lemma works for every Case A star.

All norms in this section are the infinity norm, except where a subscript 1 is written. Let $`\mathcal B_r=[-r,r]^2`$. Keep an arbitrary integer injection, $`\eta =w(\theta )`$, and $`f_i=Q_i\eta`$, without an annihilator assumption. Choose a lattice basis $`(v_i,e_i)`$ and write
```math
z=s_i(z)v_i+\pi_i(z)e_i,\qquad s_i(z)=\det(z,e_i).
```
The end fields $`g_i^\pm`$ of [the common-extrema argument](all-star-extension.md#1-common-transverse-extrema) have common nonzero transverse extrema $`\mu _i^-\le \mu _i^+`$.

## 1. Certificate data and Theorem 2

Choose integer certificates $`\alpha _i\le \gamma _i`$ and nonnegative integers $`T_i`$ with the following properties:

1. $`\mathrm{supp} f_i`$ and both end supports lie in $`\Pi_i=\{z:\alpha_i\le\pi_i(z)\le\gamma_i\}`$.
2. On $`\Pi _i`$, $`f_i(z)=g_i^+(z)`$ when $`s_i(z)>T_i`$, and $`f_i(z)=g_i^-(z)`$ when $`s_i(z)<-T_i`$.

The explicit construction in [the certificate construction](finite-window-bound.md#3-automatic-certificate-construction) always supplies these certificates; smaller directly verified ones are allowed. Put
```math
R_0=\max_i\max_{\substack{s\in\{-T_i,T_i\}\\t\in\{\alpha_i,\gamma_i\}}}
\|sv_i+te_i\|_\infty,
\qquad
c_i=\max\{|\alpha_i-\mu_i^-|,|\gamma_i-\mu_i^-|\},
```
```math
K=\left\lceil\max_{i\lt j}
\frac{\|v_j\|_\infty c_i+\|v_i\|_\infty c_j}
{|\det(v_i,v_j)|}\right\rceil,
\qquad R=R_0+K+1.
\qquad\text{(F1)}
```
All these numbers are fixed integers depending only on the supplied star decomposition and encoding.

Choose two indices $`i\ne j`$ and nonzero points $`z_i,z_j`$ of $`f_i,f_j`$. Set
```math
d=z_j-z_i,\quad M(z)=f_i(z)f_j(z+d),\quad
W=\mathrm{supp}Q_i\cup(d+\mathrm{supp}Q_j),
```
```math
F_0=\{0\}\cup(\mathrm{supp}M+W).
\qquad\text{(F2)}
```
As proved in [the complete-marker argument](all-star-extension.md#4-a-complete-finite-marker), $`M`$ has nonempty finite support and is a local function with window $`W`$. An explicitly computed finite superset of $`\mathrm{supp} M`$ may be used in (F2).

For a nonempty finite window $`S`$, use the erosion $`E_i(S)`$ and the complete-period row set $`L_i(S)`$ from the Case B anchor lemma, but use the **end** extrema even if they differ from the global extrema of $`f_i`$. Define
```math
C_i^{\rm end}(S)=(\mu_i^- -L_i(S))\cap(\mu_i^+-L_i(S)),
```
```math
\mathcal A^{\rm end}(S)=\{u\in\mathbb Z^2:
\pi_i(u)\in C_i^{\rm end}(S)\text{ for at least two }i\},
```
```math
\mathcal J(S)=\bigcap_{f\in F_0}(f-S),\qquad
\mathcal O_R(S)=(\mathcal B_R-\mathrm{Conv}S)\cap\mathbb Z^2.
\qquad\text{(F3)}
```
These are finite explicitly enumerable sets. In particular $`0\in F_0`$ implies $`\mathcal J(S)\subseteq -S`$. The occurrence of $`S`$, rather than its convex hull, in $`\mathcal J(S)`$ is essential for nonconvex windows.

**Theorem 2 (all-star finite bound).** For every CN star and every nonempty finite $`S`$,
```math
\boxed{
P_\theta(S)\ge
\left|\mathcal J(S)\ \cup\ 
\left(\mathcal A^{\rm end}(S)\setminus\mathcal O_R(S)\right)\right|.
}
\qquad\text{(F4)}
```
Every anchor counted on the right has a globally unique $`S`$-pattern. In particular,
```math
P_\theta(S)\ge |\mathcal A^{\rm end}(S)|
-|\mathcal O_R(S)\setminus\mathcal J(S)|.
\qquad\text{(F5)}
```
No Case A/B test occurs in these statements. On a lattice-convex window one may also take the maximum with CN's $`|S|+1`$ bound. The union expression (F4) remains nonnegative even when the weaker difference in (F5) is negative.

For a lattice-convex $`S=C\cap \mathbb Z^2`$, this also supplies a purely geometric finite correction. Let $`L=\max_{f\in F_0}\|f\|_\infty`$, and write $`D\ominus\mathcal B_L=\{x:x+\mathcal B_L\subseteq D\}`$. Then
```math
P_\theta(S)\ge |\mathcal A^{\rm end}(S)|-
\#\left[\left(((-C)+\mathcal B_R)\setminus((-C)\ominus \mathcal B_L)\right)\cap\mathbb Z^2\right].
\qquad\text{(F10)}
```
Indeed $`((-C)\ominus \mathcal B_L)\cap \mathbb Z^2\subseteq \mathcal J(S)`$. This is an explicit shell count, with fixed, computable $`R,L`$, not an unspecified boundary error. It is permitted to be conservative on small windows.


**Proof.** Newton-polygon additivity gives $`0\in\mathrm{conv}(\mathrm{supp}Q_i)=\sum_{j\ne i}[0,H_j]`$. Therefore $`E_i(S)\subseteq\mathrm{conv}S`$, even for a window with holes. The origin coefficient of $`Q_i`$ may cancel, as shown in [objects](objects.md#3-observation-centers-and-newton-polygons). This is why we use the convex hull rather than assert that observation centers belong to $`S`$.

Suppose first that $`u\in \mathcal A^{\mathrm{end}}(S)\setminus \mathcal O_R(S)`$. The convex region $`u+\mathrm{conv}(S)`$ avoids $`\mathcal B_R`$, hence avoids the central rectangles enclosed by $`\mathcal B_{R_0}`$. On its intersection with $`\Pi _i`$, one common end profile describes all observations of $`f_i`$: a segment between the two ends would cross the excluded rectangle. The complete-period runs on the two extreme rows therefore give the extrema of equation (5) in [the extremal-row proof](signals-and-anchors.md#2-extremal-row-recovery) in each good direction.

If $`u'`$ has the same pattern, a nonzero observation on the lowest row gives
$`\alpha _i-\mu _i^-\le \pi _i(u'-u)\le \gamma _i-\mu _i^-`$, hence $`|\pi _i(u'-u)|\le c_i`$. For two good directions the inverse coordinate formula is
```math
h=\frac{\pi_i(h)v_j-\pi_j(h)v_i}{\det(v_i,v_j)}.
```
Consequently $`\|u'-u\|_\infty \le K`$ by (F1). The second convex region also avoids $`\mathcal B_{R_0}`$, since the first avoids $`\mathcal B_{R_0+K+1}`$. Its observed nonzero heights are contained in the interval with the same end extrema $`\mu _i^\pm`$. The equal-length argument of [extremal recovery](signals-and-anchors.md#2-extremal-row-recovery) now gives equality of both transverse coordinates, and $`u'=u`$.

If $`u\in \mathcal J(S)`$, all local windows needed to observe the full translated nonzero support of $`M`$ fit inside the actual $`S`$. The complete-marker argument of [the complete-marker argument](all-star-extension.md#4-a-complete-finite-marker) gives $`u'=u`$ for any equal pattern. These two globally unique families may therefore be united, proving (F4); subtracting at most $`|\mathcal O_R(S)\setminus\mathcal J(S)|`$ points from $`\mathcal A^{\mathrm{end}}`$ proves (F5). ∎

## 2. A Case A numerical example

Take
```math
\theta(x,y)=\mathbf1[y\ge0]+\mathbf1[x\ge0]\pmod2,
\qquad v_1=e_1,\ v_2=e_2,\ H_i=v_i,
```
and the natural integer injection. This is Case A, since $`D\eta (-1,-1)=-2`$. Directly,
```math
f_1=(1-2\mathbf1[x\ge0])\mathbf1[y=-1],\qquad
f_2=(1-2\mathbf1[y\ge0])\mathbf1[x=-1].
```
Use bases $`e_1^{\mathrm{basis}}=e_2`$, $`e_2^{\mathrm{basis}}=-e_1`$, and certificates
```math
(\alpha_1,\gamma_1)=(-1,-1),\quad
(\alpha_2,\gamma_2)=(1,1),\quad T_1=T_2=1.
```
Formula (F1) gives $`R_0=1`$, $`c_1=c_2=0`$, $`K=0`$, $`R=2`$. Both fields are nonzero at $`(-1,-1)`$, so choose $`d=0`$. Then
```math
M=\mathbf1[\{(-1,-1)\}],\quad W=\{0,e_1,e_2\},\quad
F_0=\{0,(-1,-1),(-1,0),(0,-1)\}.
```
For the standard triangle $`S_n`$, $`n\ge 4`$, the good anchors form $`[-n,-1]^2`$. The central count is $`n(n-1)/2`$. Its region is disjoint from the distant anchors, since it contains the origin in the translated window. The latter count is $`(n-4)(n-3)/2`$: avoidance of $`\mathcal B_2`$ within the good square is equivalent to
$`u_x,u_y\le -3`$ and $`u_x+u_y\le -n-5`$.
Thus (F4) gives the explicit Case A inequality
```math
\boxed{P_\theta(S_n)\ge n^2-4n+6\qquad(n\ge4).}
\qquad\text{(F6)}
```
At $`n=30`$ this is 786, strictly above CN's 497 on the same 496-site triangle. No Case B hypothesis is used.

## 3. Automatic certificate construction

This makes explicitness independent of guessing optimized certificates. The input consists of the chosen periods, bases, tail lattices/tables and the finite transition data of a known star decomposition.

For each $`i`$, form the finite formal exponent set
$`E_i^{\mathrm{formal}}=\{\sum_{j\in C}H_j:C\subseteq[m]\setminus\{i\}\}`$ and, for each component $`k`$, put
```math
m_{ik}=1+\max_{q\in E_i^{\rm formal}}|\pi_k(q)|,
\quad
b_{ik}=\max\{|\ell_k-m_{ik}|,|r_k+m_{ik}|\},
```
```math
A_i=\left\lceil\max_{j\lt k}
\frac{\|v_k\|_\infty b_{ij}+\|v_j\|_\infty b_{ik}}
{|\det(v_j,v_k)|}\right\rceil.
\qquad\text{(F7)}
```
[The widened-template argument](all-star-extension.md#2-bounded-cores-and-end-agreement) puts $`\mathrm{supp} f_i`$ inside its own widened strip or a pairwise strip intersection. The latter lies in $`\mathcal B_{A_i}`$, by the inverse coordinate formula. Therefore valid global symmetric bounds are
```math
B_i=\max\{b_{ii},\ \|v_i\|_1A_i,
|\mu_i^-|,|\mu_i^+|\},\qquad
\alpha_i=-B_i,\quad\gamma_i=B_i.
\qquad\text{(F8)}
```
Here $`|\pi _i(z)|\le \|v_i\|_1\|z\|_\infty`$ explains the middle term.

The end extrema themselves require only finite evaluations. Put
```math
\beta_i^- =\min_{q\in E_i^{\rm formal}}\pi_i(q),\qquad
\beta_i^+=\max_{q\in E_i^{\rm formal}}\pi_i(q).
```
Both $`g_i^\pm`$ have support only on integer rows in
$`[\ell _i-\beta _i^+,r_i-\beta _i^-]`$. Evaluate $`g_i^+`$ at $`sv_i+te_i`$ for those finitely many rows and $`0\le s<\kappa _i`$. Its nonzero rows give $`\mu _i^\pm`$; their equality for the other end was proved in equation (24) in [the common-extrema proof](all-star-extension.md#1-common-transverse-extrema). Nonvanishing guarantees that this finite computation finds at least one row.

Set
```math
T_i=1+\max_{j\ne i}\left\lceil
\frac{B_i|\pi_j(e_i)|+
\max_{q\in E_i^{\rm formal}}|\pi_j(q)|+
\max(|\ell_j|,|r_j|)+1}
{|\pi_j(v_i)|}
\right\rceil.
\qquad\text{(F9)}
```
For $`z=sv_i+te_i`$, $`|t|\le B_i`$ and $`|s|>T_i`$, every $`z+q`$ in the template is in the tail of component $`j\ne i`$ prescribed by $`\mathrm{sign}(s\pi _j(v_i))`$. This follows from the triangle inequality in the numerator of (F9). Hence the required end agreement is certified. Equations (F1) now compute $`R`$ using only integers and rational ceilings.

Finally select a nonzero sampled row/phase of $`g_i^+`$ and translate it by a sufficiently large nonnegative multiple of $`H_i`$ so its tangential coordinate exceeds $`T_i`$. This gives an explicit nonzero point of $`f_i`$. Do this in two directions to construct $`M`$. Its support is in the bounded integer parallelogram $`\Pi _i\cap (\Pi _j-d)`$, which can be enumerated and evaluated; a finite enclosing rectangle also suffices. This computes $`F_0`$. Thus (F4) has no unknown asymptotic constants or hidden pattern-complexity oracle.


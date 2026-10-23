# Signals and globally unique anchors

Prerequisites: [objects and statements](objects.md). Throughout this file assume Case B.

## 1. Nonzero strip signals in Case B

We reprove exactly the part of CN Lemma 4.1 used here (p. 13). The proof uses an arbitrary injective integer encoding.

Fix $`i`$. As $`n\to +\infty`$, $`T^{nH_i}\theta`$ eventually stabilizes on each finite set to
```math
\Xi_i=F_i+B_i^+,
```
where $`B_i^+`$ is the sum of the tails of the other components selected by the signs of $`\pi _j(H_i)`$. Indeed, these projections are nonzero, and $`H_i`$ preserves its component and all the tails. Write
```math
G_i^R=R_i+B_i^+,\quad G_i^L=L_i+B_i^+,
\quad\delta_i^{R/L}=w(\Xi_i)-w(G_i^{R/L}).
```
Since $`(T^{H_i}-1)f_i=D\eta =0`$, $`f_i`$ is $`H_i`$-periodic. Taking the isolating limit in $`f_i=T^{nH_i}f_i`$, a finite-window identity, gives
```math
f_i=Q_iw(\Xi_i)=Q_i\delta_i^R=Q_i\delta_i^L.
\qquad\text{(3)}
```
Here every factor of $`Q_i`$ preserves each encoded $`G_i^{R/L}`$; $`m\ge 2`$ ensures that $`Q_i`$ has at least one factor. The right difference vanishes on $`\pi _i>r_i`$, and the left difference vanishes on $`\pi _i<\ell _i`$. Finitely many shifts in (3) confine $`\mathrm{supp} f_i`$ to a strip of finite width in direction $`v_i`$.

Moreover, $`\delta _i^R\ne 0`$. Otherwise injectivity of $`w`$ would give $`F_i=R_i`$, contradicting non-double-periodicity. Its nonzero transverse heights have a maximum $`t_*`$. Choose $`z_*`$ there with $`\delta _i^R(z_*)\ne 0`$. In the formal expansion
```math
Q_i=\sum_{C\subseteq[m]\setminus\{i\}}
(-1)^{m-1-|C|}T^{H_C},\qquad H_C=\sum_{j\in C}H_j,
```
the unique minimum of $`\pi _i(H_C)`$ occurs at
```math
C_* =\{j\ne i:\pi_i(H_j)<0\}.
```
Changing any membership strictly increases the projection. At $`z_*-H_{C_*}`$, only that term of $`Q_i\delta _i^R`$ can be nonzero, and it equals $`\pm \delta _i^R(z_*)`$. Thus $`f_i\ne 0`$. This uses subsets in the formal expansion; coincidences between other subset sums do not affect the unique extreme term.

Define the global transverse support extrema by
```math
\mu_i^-=\min\pi_i(\mathrm{supp}f_i),\qquad
\mu_i^+=\max\pi_i(\mathrm{supp}f_i).
```
Consequently the integers $`\mu _i^-\le \mu _i^+`$ in the Case B anchor lemma exist. We have used only the star hypotheses, an integer injection and $`D\eta =0`$; in particular, no KS ideal-structure theorem has entered.

## 2. Extremal-row recovery

For an anchor $`u`$, the $`S`$-pattern determines all values
```math
f_i(u+s)=(Q_i\eta)(u+s),\qquad s\in E_i(S).
```
Define the finite set of observed nonzero heights
```math
Y_i(u,S)=\{\pi_i(s):s\in E_i(S),\ f_i(u+s)\ne0\}.
```
Regardless of visibility,
```math
Y_i(u,S)\subseteq
[\mu_i^- -\pi_i(u),\mu_i^+-\pi_i(u)].
\qquad\text{(4)}
```

Suppose $`\pi _i(u)\in C_i(S)`$. At each of the two local heights $`\mu _i^--\pi _i(u)`$ and $`\mu _i^+-\pi _i(u)`$, the erosion contains $`\kappa _i`$ consecutive lattice points. After translation by $`u`$, these are on an extreme nonzero row of $`f_i`$. That row is periodic with period $`\kappa _i v_i`$ and is not identically zero. Any $`\kappa _i`$ consecutive points visit every residue class of that period, so at least one of the observations is nonzero. It follows that
```math
\min Y_i(u,S)=\mu_i^- -\pi_i(u),\qquad
\max Y_i(u,S)=\mu_i^+ -\pi_i(u).
\qquad\text{(5)}
```
This also covers $`\mu _i^-=\mu _i^+`$.

If another anchor $`u'`$ has the same $`S`$-pattern, its $`Y_i`$ set is identical. The interval in (4) for $`u'`$ has length $`\mu _i^+-\mu _i^-`$, exactly the distance between the two extrema in (5). An interval of that length containing both extrema must have the same endpoints. Thus
```math
\pi_i(u')=\pi_i(u).
\qquad\text{(6)}
```
The second anchor was not assumed to be good, or even to have a second visible direction.

For $`u\in \mathcal A(S)`$, apply (6) in two distinct directions. Their linear forms are independent, hence $`u'=u`$. This proves the Case B anchor lemma. The map $`(\pi _i,\pi _j)`$ need not be unimodular: independence suffices for uniqueness. Integrality restrictions must be retained when counting possible anchors.


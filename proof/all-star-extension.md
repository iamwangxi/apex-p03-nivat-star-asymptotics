# All-star extension of the exact leading term

Prerequisites: [objects](objects.md), [signals and anchors](signals-and-anchors.md), and [the geometric estimates and upper bound](leading-term.md).

**Theorem 1 (all stars).** For every CN star, without imposing (B), every fixed compact convex $`P`$ with nonempty interior satisfies the formula in [Theorem 1](objects.md#2-main-results), with exactly the same coefficient $`\Lambda (P)`$ defined from its component directions. The uncorrected Case B anchor bound is not asserted for general Case A windows.

[The upper bound](leading-term.md#3-the-matching-upper-bound) already applies. We now give a lower-bound argument for every star without a Case A/B assumption. This argument is necessary because $`Q_i\eta`$ need not be $`H_i`$-periodic globally. The central finite marker below uses the strip-product mechanism at the start of CN Lemma 4.2, after establishing its strip hypotheses for all stars.

## 1. Common transverse extrema

Keep any integer injection and write $`f_i=Q_i\eta`$, without assuming $`D\eta =0`$. For a sign $`\sigma\in\{+,-\}`$, let $`B_i^\sigma`$ be the sum of the other components' tails, choosing $`R_j`$ when $`\sigma\pi_j(v_i)>0`$ and $`L_j`$ when $`\sigma\pi_j(v_i)<0`$. These are the limits obtained by translation along $`\sigma H_i`$. For both isolating limits $`\Xi _i^\sigma =F_i+B_i^\sigma`$, $`\sigma \in \{+,-\}`$, define
```math
g_i^\sigma=Q_iw(\Xi_i^\sigma).
```
Each $`g_i^\sigma`$ is $`H_i`$-periodic. It is nonzero and supported in a bounded-width strip, by applying the two difference identities and unique-extreme argument of [the strip-signal argument](signals-and-anchors.md#1-nonzero-strip-signals-in-case-b) directly to $`\Xi _i^\sigma`$; that argument does not require $`D\eta =0`$ when applied to the isolated configuration itself.

We need the following stronger detail. Define
```math
a_i=\min\pi_i\{z:F_i(z)\ne L_i(z)\},\qquad
b_i=\max\pi_i\{z:F_i(z)\ne R_i(z)\}.
```
Both sets are nonempty, since otherwise $`F_i`$ would be a doubly periodic tail. Their indicated extrema exist by the two half-plane agreements. In the formal support of $`Q_i`$ put
```math
\beta_i^- =\sum_{j\ne i,\ \pi_i(H_j)<0}\pi_i(H_j),\qquad
\beta_i^+=\sum_{j\ne i,\ \pi_i(H_j)>0}\pi_i(H_j).
```
Then, for both signs,
```math
\min\pi_i(\mathrm{supp}g_i^\sigma)=a_i-\beta_i^+ =:\mu_i^-,
\quad
\max\pi_i(\mathrm{supp}g_i^\sigma)=b_i-\beta_i^- =:\mu_i^+.
\qquad\text{(24)}
```
To verify this, subtract the encoded left or right pure tail using equation (3) in [the strip-signal proof](signals-and-anchors.md#1-nonzero-strip-signals-in-case-b). Because addition in $`\mathbb F_p`$ can be cancelled and $`w`$ is injective, the support of
$`w(F_i+B_i^\sigma )-w(L_i+B_i^\sigma )`$ is exactly $`\{F_i\ne L_i\}`$, independent of $`\sigma`$; the corresponding right support is exactly $`\{F_i\ne R_i\}`$. The unique largest projection in the formal $`Q_i`$ expansion leaves one nonzero term at the lowest possible output row $`a_i-\beta _i^+`$. All rows below it vanish. The unique smallest projection similarly gives the highest row $`b_i-\beta _i^-`$ using the right difference. This proves (24), including possible coincidences of nonextreme formal exponents. Thus the two end profiles may have different values, but they have the same nonzero transverse extrema.

## 2. Bounded cores and end agreement

For each $`i`$ there is a bounded strip
```math
\Pi_i=\{z:\alpha_i\le\pi_i(z)\le\gamma_i\}
```
containing $`\mathrm{supp} f_i`$ and both profile supports. There is also a fixed disk containing a bounded core, outside which $`f_i`$ agrees with the appropriate positive or negative end profile on $`\Pi _i`$ and is zero off $`\Pi _i`$.

Here are details without an annihilator assumption. Let $`E_i^{\mathrm{formal}}=\{H_C:C\subseteq[m]\setminus\{i\}\}`$. For each component $`k`$ take $`M_k=1+\max_{q\in E_i^{\mathrm{formal}}}|\pi _k(q)|`$ and the widened strip $`\Sigma_k=\{z:\ell_k-M_k\le\pi_k(z)\le r_k+M_k\}`$. Outside $`\Sigma _k`$, the entire formal template is in a single pure tail of component $`k`$; this remains valid when the original transition strip is empty. If at most one widened strip is active and its index is $`j\ne i`$, pair the formal terms using the factor $`T^{H_j}-1`$ in $`Q_i`$: the sole varying component and every other pure tail are preserved by $`H_j`$, so the output is zero. The same conclusion holds with no active component. The region where at least two widened strips are active is a finite union of bounded pairwise intersections. Therefore, apart from a bounded set, $`\mathrm{supp} f_i`$ lies in the widened $`i`$-th strip. Enlarge that strip to include the bounded exception and the profile supports, obtaining $`\Pi _i`$.

Complete $`v_i`$ to a lattice basis $`(v_i,e_i)`$ and write $`z=s v_i+t e_i`$. For bounded $`t`$, as $`s\to \sigma \infty`$ all other components, on the whole finite $`Q_i`$ template, agree with exactly the tails defining $`B_i^\sigma`$. Hence $`f_i(z)=g_i^\sigma (z)`$ for all sufficiently large $`\sigma s`$, uniformly over the finitely many integer heights in $`\Pi _i`$. There are only finitely many remaining points of $`\Pi _i`$ with bounded $`s`$. This proves the core assertion. A sufficiently large disk also contains the real rectangle in these coordinates between the two ends of $`\Pi _i`$.

If a convex region $`C`$ avoids that disk, $`C\cap \Pi _i`$, when nonempty, lies entirely at one end: otherwise its segment joining a positive-end point to a negative-end point would cross the excluded central rectangle. Thus one common $`g_i^\sigma`$ describes every observation of $`f_i`$ in $`C\cap \Pi _i`$. In particular, outside the core all nonzero values of $`f_i`$ have transverse height between the same $`\mu _i^-`$ and $`\mu _i^+`$ from (24).

## 3. Distant anchors

Use [the fixed-row erosion lemma](leading-term.md#1-uniform-interior-rows) with the $`\mu _i^\pm`$ from (24), and define $`I_i^-(n)`$ by equation (10) in [the geometric proof](leading-term.md#1-uniform-interior-rows). Let
```math
\mathcal G_n=\{u\in\mathbb Z^2:\pi_i(u)\in I_i^-(n)
\text{ for at least two distinct }i\}.
```
By [the lattice estimate](leading-term.md#2-counting-band-intersections),
```math
|\mathcal G_n|=\Lambda(P)n^2+O(n).
\qquad\text{(25)}
```
There is a fixed radius $`R`$ such that every $`u\in \mathcal G_n`$ with $`(u+nP)\cap B_R(0)=\varnothing`$ has a globally unique $`S_n`$-pattern.

Newton-polygon additivity over the Laurent polynomial domain gives
```math
0\in\mathrm{conv}(\mathrm{supp}Q_i)
=\sum_{j\ne i}[0,H_j].
```
Taking convex hulls of $`s+\mathrm{supp}Q_i\subseteq S`$ gives $`s\in\mathrm{conv}S`$, hence $`E_i(S)\subseteq\mathrm{conv}S`$ by taking convex hulls. In particular all the difference centers used below lie in the translated convex window; the origin need not be an actual exponent of the simplified polynomial. See the cancellation example in [objects](objects.md#3-observation-centers-and-newton-polygons).

To prove the distant-anchor assertion, first choose a radius containing all the cores from [the bounded-core argument](all-star-extension.md#2-bounded-cores-and-end-agreement). For such a distant window and a good direction $`i`$, its observed set $`Y_i(u,S_n)`$ has exactly the extrema of equation (5) in [the extremal-row proof](signals-and-anchors.md#2-extremal-row-recovery): all observations are in a single periodic end profile, and the erosion supplies a complete $`\kappa _i`$-point run on each extreme row.

Suppose an arbitrary anchor $`u'`$ has the same pattern, and put $`h=u'-u`$. An observed nonzero value on the row $`\mu _i^--\pi _i(u)`$ remains nonzero at $`u'`$. The global strip bound for $`f_i`$ therefore gives
```math
\alpha_i-\mu_i^-\le\pi_i(h)\le\gamma_i-\mu_i^-.
\qquad\text{(26)}
```
For the two good directions these are two independent linear inequalities of bounded width. Consequently $`\|h\|\le K`$ for a fixed constant $`K`$, taking the maximum over the finitely many pairs. This conclusion does not assume that $`u'`$ is good or distant.

Now choose $`R`$ larger than the common core radius plus $`K`$. The translated region $`u'+nP=(u+nP)+h`$ then also avoids every core. Its observed nonzero heights lie in
```math
[\mu_i^- -\pi_i(u'),\mu_i^+-\pi_i(u')],
```
regardless of which end profile it sees, by (24). This interval has the same length as the two observed extrema at $`u`$. Exactly the argument of [extremal recovery](signals-and-anchors.md#2-extremal-row-recovery) gives $`\pi _i(u')=\pi _i(u)`$ in both good directions, hence $`u'=u`$. This proves the distant-anchor assertion. The bounded-displacement step (26) is essential; periodicity of $`f_i`$ was not assumed globally.

## 4. A complete finite marker

Choose two distinct indices $`i,j`$. Both global fields $`f_i,f_j`$ are nonzero: [the bounded-core argument](all-star-extension.md#2-bounded-cores-and-end-agreement) makes them agree far along their strips with the nonzero periodic end profiles from [the common-extrema argument](all-star-extension.md#1-common-transverse-extrema). Pick $`z_i,z_j`$ where they are nonzero and put $`d=z_j-z_i`$. Define
```math
M(z)=f_i(z)f_j(z+d),\qquad
W=\mathrm{supp}Q_i\ \cup\ (d+\mathrm{supp}Q_j).
```
Then $`M(z_i)\ne 0`$. Its support is contained in $`\Pi _i\cap (\Pi _j-d)`$, a bounded intersection of nonparallel strips; hence $`M`$ has finite nonempty support. It is a local function of $`\theta`$ with window $`W`$. This is the finite-strip product construction used in the first part of CN Lemma 4.2 (pp. 13–14); we do not need to apply $`D`$ to that product or localize a quadratic displacement.

Let
```math
F_0=\{0\}\cup(\mathrm{supp}M+W),\qquad
\mathcal I_n=\{u\in\mathbb Z^2:F_0\subseteq u+nP\}.
```
All points in $`F_0`$ are integral. For every $`u\in \mathcal I_n`$, the $`S_n`$-pattern therefore determines the complete translated nonzero support $`\mathrm{supp} M-u`$ among its $`M`$-observations: every required local window lies in $`S_n`$. If another anchor has the same pattern, it has exactly these $`|\mathrm{supp} M|`$ nonzero observed positions. Since the global $`M`$ has only that many nonzero positions, it too must show its entire translated support there. Hence $`\mathrm{supp} M-u'=\mathrm{supp} M-u`$. A finite nonempty subset of $`\mathbb Z^2`$ cannot be invariant under a nonzero translation; thus $`u'=u`$. Central-marker patterns are globally unique as well.

Increase $`R`$ further so that $`F_0\subseteq B_R(0)`$, and put
```math
\mathcal O_n=(-nP)+B_R(0).
```
The condition $`u∉\mathcal O_n`$ is exactly that $`u+nP`$ avoids the disk. Since $`0\in F_0`$, $`\mathcal I_n\subseteq -nP`$. Every point of $`-nP`$ at distance more than $`R`$ from its boundary belongs to $`\mathcal I_n`$. Therefore
```math
|(\mathcal O_n\cap\mathbb Z^2)\setminus\mathcal I_n|=O(n).
\qquad\text{(27)}
```
Indeed the difference lies in a fixed-width neighborhood of $`\partial (-nP)`$. A fixed planar compact convex body has finite perimeter; after dilation that perimeter is $`O(n)`$. Cover its boundary by $`O(n)`$ unit-length arc pieces. A fixed-width neighborhood is covered by $`O(n)`$ fixed-radius disks and contains only $`O(n)`$ lattice points. This justification also applies to nonpolygonal $`P`$.

Every anchor in
```math
\mathcal I_n\ \cup\ \bigl(\mathcal G_n\setminus\mathcal O_n\bigr)
```
has a globally unique pattern, by the two preceding arguments. Equations (25)–(27) give
```math
P_\theta(S_n)\ge |\mathcal G_n|
-|(\mathcal O_n\cap\mathbb Z^2)\setminus\mathcal I_n|
\ge\Lambda(P)n^2-O(n).
```
Together with [the general upper bound](leading-term.md#3-the-matching-upper-bound), this proves Theorem 1. ∎


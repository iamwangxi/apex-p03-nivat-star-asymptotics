# The exact linear coefficient in a three-line example


Use the binary configuration (with the displayed sum also equal to its integer encoding)
```math
c(x,y)=\sum_{i=0}^2\mathbf1[y-2ix=i],\qquad v_i=(1,2i),
```
and the fixed triangle $`P=\mathrm{conv}\{(0,0),(1,0),(0,1)\}`$, $`S_n=nP\cap\mathbb Z^2`$. The lines meet at real first coordinate $`-1/2`$, so their lattice supports are disjoint. Each component is singly periodic with zero tails, $`H_i=v_i`$, and all color indicators are killed by $`D`$. Thus this is a CN star in Case B.

The projection bands are
```math
0\le y\le1,\qquad -2\le y-2x\le1,\qquad -4\le y-4x\le1.
```
Their pairwise intersection areas are $`3/2,5/4,15/2`$. The first and third intersection is contained in the middle band, so the triple intersection has area $`5/4`$. Consequently
```math
\Lambda(P)=\frac32+\frac54+\frac{15}2-2\frac54=\frac{31}4.
```
For the general area formula and inclusion-exclusion identity, see [equations (16)–(17) in the geometric proof](leading-term.md#4-evaluating-the-area-coefficient). Here $`W_0,W_1,W_2`$ are $`1,3,5`$, while the absolute pair determinants are $`2,4,2`$, giving the three pair areas above. Since $`\pi_1=(\pi_0+\pi_2)/2`$, the first and third bands force $`-2\le\pi_1\le1`$. Their intersection is therefore contained in the middle band, proving that the triple area is $`5/4`$. This computes the coefficient in [Theorem 1](objects.md#2-main-results). We now determine the exact linear coefficient for this particular example.

**Theorem 3.** For the three-line configuration above and $`S_n=\{(x,y)\in\mathbb Z_{\ge0}^2:x+y\le n\}`$,
```math
\boxed{P_c(S_n)=\frac{31}{4}n^2+\frac{11}{4}n+O(1).}
\qquad\text{(S1)}
```
In particular the coefficient $`\beta (c,P)`$ exists and equals **11/4**. We do not infer an exact constant remainder or a general second-order formula for all stars from this example.

## 1. Sections and representations

Write
```math
Q_i(t)=\{(x,y)\in S_n:y-2ix=t\},\qquad
\mathcal Q_i=\{t:Q_i(t)\ne\varnothing\}.
```
For an anchor $`(X,Y)`$, the three section parameters are
```math
t_i=i-Y+2iX=t_0+i d,\qquad d=1+2X\text{ odd}.
\qquad\text{(S2)}
```
The real supporting lines meet at a half-integer first coordinate, so their lattice sections are disjoint. A pattern is their union, with no mod-2 cancellation.

For $`n\ge 3`$, direct section enumeration gives
```math
\mathcal Q_0=[0,n]\cap\mathbb Z,
```
```math
\mathcal Q_1=([-2n,n]\cap\mathbb Z)\setminus\{-2n+1\},
```
```math
\mathcal Q_2=([-4n,n]\cap\mathbb Z)
\setminus\{-4n+1,-4n+2,-4n+3,-4n+6,-4n+7,-4n+11\}.
\qquad\text{(S3)}
```
For verification, at $`x=n-r`$, the possible heights in $`\mathcal Q_i`$ are the integer interval
$`[-2in+2ir,-2in+(2i+1)r]`$. For $`i=1`$ all successive intervals join after $`r=1`$; for $`i=2`$ they join after the first four intervals. The gaps are exactly those in (S3). Thus the cardinalities are
```math
|\mathcal Q_0|=n+1,\quad |\mathcal Q_1|=3n,\quad |\mathcal Q_2|=5n-5.
\qquad\text{(S4)}
```

Let $`N_{ij}`$ count anchors for which sections $`i,j`$ are both nonempty. The congruences from (S2) are respectively
```math
t_1-t_0\text{ odd},\qquad
t_2-t_0\equiv2\pmod4,\qquad
t_2-t_1\text{ odd}.
```
If finite integer intervals with at most a fixed number of holes have sizes $`A,B`$, the number of pairs in two prescribed matching residue classes modulo fixed $`q`$, summed over those classes, is $`AB/q+O(1)`$. To see the bounded, rather than linear, error, write the residue counts as $`A/q+a_r`$ and $`B/q+b_r`$. The deviations are bounded and sum to zero; both cross terms cancel, leaving only $`\sum_r a_rb_{r+s}=O(1)`$.

Applying this with (S4),
```math
N_{01}=\tfrac32n^2+\tfrac32n+O(1),\qquad
N_{02}=\tfrac54n^2+O(1),
```
```math
N_{12}=\tfrac{15}{2}n^2-\tfrac{15}{2}n+O(1).
\qquad\text{(S5)}
```
For any compatible $`t_0\in \mathcal Q_0,t_2\in\mathcal Q_2`$, their midpoint is in the full interval $`[-2n,n]`$. It fails to belong to $`\mathcal Q_1`$ only if $`t_0+t_2=-4n+2`$. Because $`t_0\ge 0`$ and $`t_2\ge -4n`$, this allows at most three pairs. Consequently the triple-hit count is $`N_{012}=N_{02}+O(1)`$.

Let $`A_n`$ be the number of anchors with at least two hits. Inclusion-exclusion for these three pair sets gives
```math
A_n=N_{01}+N_{02}+N_{12}-2N_{012}
=\frac{31}{4}n^2-6n+O(1).
\qquad\text{(S6)}
```
Every pure section $`Q_i(t)`$ occurs as an actual pattern: move an anchor along $`v_i`$ and avoid the finitely many other-line intersections with the finite window. A zero pattern also occurs. Thus a complete list of representations consists of all $`A_n`$ multiple-hit anchors, one copy of every pure section, and zero. Its length is
```math
R_n=A_n+\sum_i|\mathcal Q_i|+1
=\frac{31}{4}n^2+3n+O(1).
\qquad\text{(S7)}
```
Different entries of this list may be the same pattern; this is the entire remaining issue.

## 2. Bounded exceptional collisions

Call a section **long** if it contains at least three points. A long section in direction $`i`$ is intrinsic to the pattern, including its offset: if another representation used a different parallel $`i`$-line, its other two lines could cover at most two points of this section. It therefore must use the same $`i`$-line. This also shows that the set of long directions is intrinsic.

Two long directions determine the anchor uniquely. A pure long section cannot coincide with a multiple-hit pattern: after its intrinsic line is fixed, every additional nonempty section contributes a point off that line.

Every nonempty section of size at most two lies in a fixed-size neighborhood of one of the vertices
```math
V=(0,n),\qquad U=(n,0).
```
For example the union of their infinity-radius-20 neighborhoods suffices, and has at most $`B=2\cdot 41^2`$ lattice points, independently of $`n`$. Here is a direct bound: for slope $`k=2`$ or $`4`$, if $`t\ge 0`$, a section of size at most two starts at $`(0,t)`$ with $`n-t<2(k+1)`$. If $`t<0`$, its first point has $`0\le y_0\lt k`$; the absence of a third point gives $`n-x_0\lt y_0+2(k+1)\le 3k+1`$. Both its possible points are consequently in the stated vertex neighborhoods. Horizontal short sections are in the top two rows next to $`V`$.

Patterns with no long direction have only $`O(1)`$ representations in the complete list: every short parameter is in a bounded-size set; a multiple-hit anchor is determined by two such parameters. Pure short sections also number $`O(1)`$.

Now fix a pattern with exactly one long direction $`i`$. That line and section are fixed in every representation. Let $`E`$ be the set of its other points. If two distinct representations assign an external point to the same remaining direction, that point fixes a second line and hence the anchor, a contradiction. Thus **every external point must swap between the two remaining directions**. If $`|E|\ge 3`$, two points in one old direction would have to lie on one line in a different direction, which is impossible. If $`|E|=2`$, one point belongs to each remaining direction and both swap; these two external points, lying in the fixed-size vertex neighborhoods, determine both anchors. There are only $`O(1)`$ such representations. These bounds are uniform in $`n`$: no-long representations contribute at most $`3B^2+3B+1`$, and the two-external-point cases contribute at most $`6B^2`$ after choosing their points, direction and assignment. A very loose bound $`20B^2`$ for all this exceptional excess is therefore sufficient.

It remains to count collisions with exactly one external point.

## 3. The unique linear-size collision family

Let that point be $`p=(a,b)`$, and let the intrinsic long parameter be $`t_i=t`$. It must be a singleton section in each of the two alternative directions. The following parity check rules out long directions 0 and 2.

- For long direction 0, assigning $`p`$ to direction 1 requires $`b-t`$ odd. Assigning it to direction 2 requires $`b-t\equiv2\pmod4`$. These cannot both hold.
- For long direction 2, assigning $`p`$ to direction 0 requires $`t-b\equiv2\pmod4`$. Assigning it to direction 1 requires $`t-b`$ odd. Again they cannot both hold.
- For long direction 1, the two assignments both require $`t-b`$ odd.

In the last case one alternative direction is horizontal. Its only singleton section is $`Q_0(n)=\{V\}`$, so $`p=V=(0,n)`$. The other singleton is $`Q_2(n)=\{V\}`$. The two possible parameter triples are exactly
```math
(t_0,t_1,t_2)=(n,t,2t-n),\qquad (2t-n,t,n).
\qquad\text{(S8)}
```
They give two distinct integer anchors precisely when $`t-n`$ is odd; distinctness also follows from that parity. They produce the same pattern $`Q_1(t)\cup\{V\}`$ precisely when their unused section is empty. Accordingly the number of these two-element collision classes is
```math
F_n=\#\{t\in\mathcal Q_1:
|Q_1(t)|\ge3,\ t-n\text{ odd},\ 
2t-n\notin\mathcal Q_0\cup\mathcal Q_2\}.
\qquad\text{(S9)}
```
No third representation exists, since the external point has only these two assignments. Different $`t`$ give different intrinsic long lines, so different classes.

As $`t\le n`$, absence from $`\mathcal Q_0`$ forces $`2t-n<0`$. By (S3), absence from $`\mathcal Q_2`$ then means either
```math
t<-3n/2,
```
or one of only six exceptional values $`2t-n=-4n+\delta`$, $`\delta\in\{1,2,3,6,7,11\}`$. In the main interval $`[-2n,-3n/2)`$, the long-section condition removes only a bounded number of parameters near $`-2n`$: for slope 2, all sections are long once their parameter is at least 15 above that endpoint. The parity condition retains one half of the remaining integers. Hence
```math
F_n=\frac n4+O(1).
\qquad\text{(S10)}
```
For instance the direct comparison just given yields the loose uniform estimate $`|F_n-n/4|\le 20`$ for $`n\ge 3`$.

All collisions outside (S8) have bounded total excess by Section 2 of this file. Since each class in (S8) has two representations, it subtracts exactly one from $`R_n`$. Therefore
```math
P_c(S_n)=R_n-F_n-O(1)
=\frac{31}{4}n^2+\frac{11}{4}n+O(1),
```
proving Theorem 3. ∎

The linear coefficient has the following decomposition for this example:
```math
\beta(c,P)=\underbrace{-6}_{\text{multiple-hit anchors}}
+\underbrace{9}_{\text{pure line sections}}
-\underbrace{\tfrac14}_{\text{top-vertex ambiguity}}
=\tfrac{11}{4}.
\qquad\text{(S11)}
```
One can make the bounded remainder entirely explicit, albeit very conservatively. The residue deviations in (S5) are bounded by 4, 28 and 28; the midpoint exclusion has at most three pairs. These give an error below 72 in (S7). With $`|F_n-n/4|\le 20`$ and the exceptional bound in Section 2 of this file, a valid bound for $`n\ge 3`$ is
```math
\left|P_c(S_n)-\tfrac{31}{4}n^2-\tfrac{11}{4}n\right|
\le20(2\cdot41^2)^2+100.
\qquad\text{(S12)}
```
No numerical fit or unproved periodicity of the constant term is needed for this bound.



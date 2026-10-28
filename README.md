# Exact complexity asymptotics and eventual quasipolynomials for Nivat stars

Every star in the specified Convex Nivat paper has an exact area coefficient on arbitrary fixed convex bodies and on finite unions of convex bodies with nonempty interior, including disconnected windows and holes. Actual directional projection sets determine the coefficient and its area-equality cases. On rational polygons complexity is eventually a quadratic quasipolynomial. On every integral convex polygon its linear coefficient is constant for **every** star, while the constant term can have minimal period four. Finite constructions give exact coefficients and decide their minimal periods.

## Results and review status

The target is Apex Intelligence, *The Convex Nivat Conjecture: A Complexity Lower Bound for Star Configurations, and a Reduction from Low Convex Complexity to Star Configurations*, contest version of 12 September 2026. Its PDF SHA-256 is `7fd67831155f4226c010fd6af32de558b76643771e6be12dea4eebafa21745a8`. Citations use printed pages. The package concerns quantitative results for stars, not a stronger periodicity threshold for arbitrary configurations.

For primitive component directions let $`\pi_i(z)=\det(v_i,z)`$, $`A_i=\pi_i(P)`$, $`B_i=\pi_i^{-1}(A_i)`$, and $`\Lambda(P)=\mathrm{area}(\bigcup_{i\lt j}(B_i\cap B_j))`$.

- **Convex area and finite-window results:** $`P_\theta(nP\cap\mathbb Z^2)=\Lambda(P)n^2+O(n)`$ for every star and fixed convex body with interior; an explicit finite-window lower bound counts globally unique marker and two-direction anchors.
- **Nonconvex windows, established in one fresh-context AI adversarial review:** [directional projections and equality classification](proof/nonconvex-windows.md). The same exact leading term holds for every finite union of full-dimensional compact convex bodies. Projection gaps must be retained. Always $`\Lambda(P)\ge\mathrm{area}(P)`$; the inequality is strict for at least three directions. For two directions equality holds exactly for $`P=\{x:\pi_1(x)\in A_1,\ \pi_2(x)\in A_2\}`$, the inverse images of products of finite interval unions. Connected regions share their convex hull's coefficient. The page also proves an arbitrary-finite-window cross-star formula and an exact sparse two-line product count. Theorems 1–4 do not use CN Theorem T; Corollary 2.1 uses it only for aligned parallelograms.
- **Eventual exact form:** on every positive-area rational convex polygon, $`P_\theta(S_n)=\Lambda(P)n^2+\beta_{n\bmod M}n+\gamma_{n\bmod M}`$ eventually. This applies Bogart–Goodrick–Woods's existing counting theorem to finite star descriptions and definable representatives. The window-predicate construction also covers finite rational polygonal unions.
- **Integral convex polygons, established in one fresh-context AI adversarial review:** [constant linear coefficients for every star](proof/integral-windows.md). For all stars, including unequal tails, sparse components, periodic backgrounds and cancellation, $`P_\theta(S_n)=\Lambda(P)n^2+\beta n+\gamma_{n\bmod M}`$ eventually, with $`\beta\in\mathbb Q`$ and $`p_\beta=1`$. Appendix B gives a direct finite formula in the original component and tail tables. An exact binary diagonal example has $`\beta=8`$ and $`(\gamma_0,\gamma_1,\gamma_2,\gamma_3)=(2,3,2,-1)`$, with $`p_\gamma=p_f=4`$.
- **Direct full-support line formula:** for nonzero periodic words on full affine lattice lines, arbitrary rational convex polygons have the explicit coefficient $`\mathfrak B_r-\Gamma_r`$, including intersections and cancellation. A vertex-denominator modulus suffices for this linear coefficient; the full remainder may require a different period.
- **Representative coefficients, established in one fresh-context AI adversarial review:** [explicit coefficients and sufficient periods](proof/explicit-coefficients.md), through finite elimination and geometric sums, including thin and lower-dimensional cells. This is an algorithmic formula and can be large.
- **Exact periods, established in one fresh-context AI adversarial review:** [constant terms and minimal periods](proof/minimal-periods-and-constant-term.md). A sufficient period $`M=qE`$ comes from vertex-velocity denominators; finite cone sums give $`\gamma`$; cyclotomic reduction decides the exact minimal periods of $`\beta`$, $`\gamma`$ and the whole quasipolynomial. A three-line colored family with $`\beta=3-\tfrac14\mathbf1[c_0=c_2]`$ rules out formulas using only uncolored geometry or independent component terms.

The [review-status record](sources/review-status.md) identifies reviewed manuscripts, review scopes and SHA-256 fingerprints. The two added manuscripts each passed one fresh-context GPT-6.1 Sol AI adversarial review on 28 October 2026. The nonconvex review found no fatal or proof-repair issue and one dependency wording issue, now incorporated; the integral-window review found no issues in any category. AI adversarial review is not human certification or proof-assistant verification.

## Proof guide

[Objects and original quantitative statements](proof/objects.md) defines stars and conventions. The convex area and finite-window proof comprises [signals and anchors](proof/signals-and-anchors.md), [leading-term geometry](proof/leading-term.md), [the all-star extension](proof/all-star-extension.md), and [the finite-window bound](proof/finite-window-bound.md). The [nonconvex page](proof/nonconvex-windows.md) gives its full extended proof, equality rigidity and exact collision examples.

[Eventual quasipolynomiality](proof/eventual-quasipolynomial.md) gives the finite tables, the norm-first representative, the counting-theorem interface and the irrational-window obstruction. [Representative coefficients](proof/explicit-coefficients.md) and [constant terms and periods](proof/minimal-periods-and-constant-term.md) make the resulting tables effective. These counting constructions do not assume the integral-window constancy theorem.

[Full-support periodic lines](proof/full-support-lines.md) proves cone certificates, projection holes, facet terms, singleton collisions and bounded intersection corrections. [Integral windows](proof/integral-windows.md) extends the direct boundary analysis to sparse lines, finite-width equal-tail components and finally all stars; its Appendix B is the complete general proof. Theorem numbers on each page are local.

[The three-line example](proof/three-line-second-order.md) independently proves the coefficient $`11/4`$. [Comparisons and scope](proof/remarks-and-scope.md) supplies the 4059 bound, complete enumeration coverage, classical sharpness and the remaining nonconvex product exception. The [bibliography](sources/bibliography.md) identifies external sources and checked locations.

## Reproducing the computations

Python 3.9 or later and its standard library suffice. No downloads, credentials or external source files are needed. From the repository root run:

~~~sh
python3 -B scripts/verify.py
~~~

The runner executes twelve diagnostics and compares their complete deterministic outputs with eleven JSON baselines and one text baseline. Programs that save output write only an adjacent file; the new programs also print deterministic JSON. Shared helpers import only repository scripts. Timings, timestamps and source-directory paths are excluded from the new outputs.

- `area_and_anchor_audit.py`: exact rational band areas, the 4059 count and global line patterns.
- `end_profile_audit.py`: unequal periodic tails, nonlinear encodings, end extrema and origin-exponent cancellation.
- `three_line_enumeration.py`: all $`3\le n\le164`$, including exact collision multiplicities.
- `finite_window_audit.py`: all 511 nonempty subsets of a three-by-three window, 120 additional finite windows and general certificate examples.
- `full_support_lines_audit.py`: cone-derived coefficients before pattern enumeration; arbitrary-window and intersecting-line collision checks.
- `presburger_examples_audit.py`: 2,401 finite-table comparisons and complete pattern representations for the specified examples.
- `minimal_period_audit.py`: cone coefficients against independent vertical sections, including thin, empty and lower-dimensional cells; cyclic versus cyclotomic period reduction; all eight color assignments in the three-line family.
- `beta_geometry_audit.py`: finite elimination, exact affine-polygon coefficients, thin cells and lower-dimensional cases.
- `nonconvex_window_audit.py`: all 65,535 nonempty four-by-four windows for the exact cross formula; four projection-union geometries, 39 exact cross/appendage examples, and 42 further star counts.
- `two_line_product_audit.py`: eight exact product-region counts with sparse periodic words, determinants of absolute value one, two and three, both signs, intersections and cancellation; retains the necessary affine-coset correction.
- `integral_sparse_audit.py`: finite sparse-line coefficients before 126 exact counts, 18 independent support enumerations, seven rectangle counts, exact period-four cone tables and both minimal-period reductions.
- `integral_boundary_audit.py`: the direct general-tail formula against 157 frozen exact-count tables, including every auxiliary collar/catalog value; fresh exhaustive counts on a fixed coverage set. The fixture is `integral_boundary_cases.json`; helpers are `integral_boundary_coefficients.py` and `star_pattern_tools.py`. Finite-difference comparisons are diagnostics, not certified eventual coefficients.

The large exploratory screens are not bundled or rerun. Their reported scales and limits are described on the two new proof pages. In particular, the integral-window source reports 2,362,685 genuine star/window instances and 19,366,833 exact global counts, with overlap between batches; those are historical exploratory totals.

Pattern equality uses exact masks or color strings. The [proved scan rectangle](proof/remarks-and-scope.md#11-complete-line-enumeration), [quadrant threshold cells](proof/remarks-and-scope.md#12-complete-quadrant-enumeration), [full-line representation lists](proof/full-support-lines.md), and new pages' complete classifications justify coverage. No fitted polynomial or observed period is used as a proof.

`MANIFEST.sha256` fingerprints every package file except itself and `.DS_Store`:

~~~sh
shasum -a 256 -c MANIFEST.sha256
~~~

## Numerical landmarks and limits

On the three-line triangle at $`n=30`$, the window has 496 sites, the finite lower bound is 4059, CN gives 497, and exact complexity is 7037. Its eventual form is $`31n^2/4+11n/4+\gamma_{n\bmod M}`$. The distinct quadrant star has finite bound 786 and exact count 1022 on the same window.

For horizontal and vertical directions, $`P=[0,1]^2\cup[3,4]^2`$ has area two and coefficient four; filling the projection gaps gives sixteen. A square with a diagonal zero-width appendage lies outside the full-dimensional-union hypothesis and changes the quadratic coefficient through genuine pattern collisions. For general two-direction stars on disconnected product regions, the leading coefficient equals window area and the sign of the lower-order difference remains unresolved here. The sparse two-line subclass is resolved by the exact product formula.

The full cross on a half-width rational rectangle has linear coefficients $`3,5/2`$; the sparse cross has $`9/2,4`$. An irrational width can prevent any periodic linear coefficient with bounded remainder. Integral convex polygons have constant linear coefficient for all stars, but need not have constant remainder.

CN supplies the finite-difference and marker mechanisms. Quadratic growth alone was already present in earlier work, including Kari–Szabados. The contributions here are the exact area coefficient and nonconvex equality classification, the star application of existing counting theory, effective coefficients and periods, and the stated direct boundary and collision formulas. No stronger general Nivat threshold, formal certification, practical complexity bound for elimination or comprehensive priority certification is claimed.

## AI disclosure

GPT-6.1 Sol, in OpenAI Codex under human direction, developed the initial proofs and computations and the integral-window manuscript. GPT-6 Astra developed the constant-term/minimal-period supplement and the nonconvex-window manuscript, and advised on research direction and scope. Separate GPT-6.1 Sol sessions in fresh contexts performed the recorded adversarial reviews. Claude planned and edited earlier source stages and checked key source interfaces. GPT-6.1 Sol prepared this integration and the portable diagnostics. No human expert has certified the work.

## License

Original prose and selection/arrangement: CC BY 4.0. Code in `scripts/`: MIT. External publications are cited, not bundled. See [LICENSE](LICENSE).

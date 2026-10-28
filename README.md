# Exact complexity asymptotics and eventual quasipolynomials for Nivat stars

Every star in the specified Convex Nivat paper has an exact area coefficient on arbitrary fixed convex bodies with nonempty interior. On rational polygons its pattern complexity is eventually a quadratic quasipolynomial. Full-support periodically colored affine lines additionally have a direct geometric formula for the linear coefficient. For general stars, the linear and constant coefficients are explicit finite sums over representative cells, and the exact minimal period is decidable.

## Results and review status

The target is Apex Intelligence, *The Convex Nivat Conjecture: A Complexity Lower Bound for Star Configurations, and a Reduction from Low Convex Complexity to Star Configurations*, contest version of 12 September 2026. Its PDF SHA-256 is `7fd67831155f4226c010fd6af32de558b76643771e6be12dea4eebafa21745a8`. Citations use printed pages. The package concerns a quantitative supplement for stars, not a stronger periodicity threshold for arbitrary configurations.

Let $`B_i(P)`$ be the supporting band parallel to direction $`v_i`$, and $`\Lambda(P)=\mathrm{area}(\bigcup_{i\lt j}(B_i(P)\cap B_j(P)))`$.

- **Reviewed area and finite-window results:** $`P_\theta(nP\cap\mathbb Z^2)=\Lambda(P)n^2+O(n)`$ for every star and fixed convex body with interior; an explicit finite-window lower bound counts globally unique marker and two-direction anchors.
- **Reviewed eventual exact form:** on every positive-area rational polygon, $`P_\theta(S_n)=\Lambda(P)n^2+\beta_{n\bmod M}n+\gamma_{n\bmod M}`$ eventually. This is an application of Bogart–Goodrick–Woods's existing counting theorem to finite star descriptions and definable representatives.
- **Reviewed direct line formula:** for nonzero periodic words on full affine lattice lines, arbitrary rational polygons have the explicit coefficient $`\mathfrak B_r-\Gamma_r`$, including intersections and cancellation. A vertex-denominator modulus suffices for this linear coefficient; the full remainder may require a different period. Integral polygons have constant linear coefficient in this class.
- **Separate supplement, established in one fresh-context AI adversarial review:** [representative-polygon coefficients and sufficient periods](proof/explicit-coefficients.md). This covers general stars through finite elimination and a geometric sum, rather than a compact formula on the original direction bands. The other results do not use this supplement as a premise.
- **Second supplement, established in one fresh-context AI adversarial review:** [exact constant terms and minimal periods](proof/minimal-periods-and-constant-term.md). A smaller sufficient period $`M=qE`$ from vertex-velocity denominators; an explicit constant term $`\gamma`$; the exact minimal periods of $`\beta`$, of $`\gamma`$ and of the whole quasipolynomial by cyclotomic reduction; and a three-line colored family with $`\beta=3-\tfrac14\mathbf 1[c_0=c_2]`$, which rules out formulas using only uncolored geometry or independent component terms.

The [review-status record](sources/review-status.md) identifies the source and review fingerprints, including the fresh-context GPT-6.1 Sol reviews of this consolidation and the first supplement (25 October 2026) and of the second supplement (26 October 2026). No fatal or proof-repair issue was found. This is an AI adversarial review, not human certification or proof-assistant verification.

## Proof guide

[Objects and original quantitative statements](proof/objects.md) defines stars and all conventions. The complete area and finite-window proof comprises [signals and anchors](proof/signals-and-anchors.md), [leading-term geometry](proof/leading-term.md), [the all-star extension](proof/all-star-extension.md), and [the finite-window bound](proof/finite-window-bound.md).

[Eventual quasipolynomiality](proof/eventual-quasipolynomial.md) gives the full finite-table description, the unique norm-first representative, the exact counting-theorem interface, the common quadratic coefficient, and the irrational-window obstruction.

[Full-support periodic lines](proof/full-support-lines.md) includes complete proofs of the cone certificates, projection holes, two-line formula, facet contribution, singleton collision integral, and bounded intersection correction. Theorem numbers on that page are local to that module.

[The three-line example](proof/three-line-second-order.md) independently proves the coefficient $`11/4`$. [Comparisons and scope](proof/remarks-and-scope.md) supplies the 4059 bound, the complete enumeration-rectangle and quadrant-cell coverage proofs, and classical sharpness. [The bibliography](sources/bibliography.md) records exact versions and checked source locations.

The supplementary general coefficient proof is confined to its separately marked page. Its review is documented in the review-status record; its own finite diagnostics do not confer independent-review status.

## Reproducing the computations

Python 3.9 or later and its standard library suffice. No downloads, credentials or external source files are needed for the diagnostics. From the repository root run:

~~~sh
python3 -B scripts/verify.py
~~~

The runner executes eight programs and compares their complete deterministic outputs to the packaged baselines. Each program writes only its adjacent output.

- `area_and_anchor_audit.py`: exact rational band areas, the 4059 count and global line patterns.
- `end_profile_audit.py`: unequal periodic tails, nonlinear encodings, end extrema and origin-exponent cancellation.
- `three_line_enumeration.py`: all $`3\le n\le164`$, including exact collision multiplicities.
- `finite_window_audit.py`: all 511 nonempty subsets of a three-by-three window, 120 additional finite windows and general certificate examples.
- `full_support_lines_audit.py`: cone-derived coefficients before pattern enumeration; arbitrary-window and intersecting-line collision checks. The deterministic mathematical payload is preserved from the reviewed source output.
- `presburger_examples_audit.py`: 2,401 finite-table comparisons and complete pattern representations for the specified examples.
- `minimal_period_audit.py`: exact coefficients from the cone formulas compared with independently counted vertical sections, including nonprimitive normals, thin and empty cells and extrapolated vertices; cyclotomic period reduction against cyclic-table equality; the colored three-line family under all eight color assignments.
- `beta_geometry_audit.py`: the reviewed supplement's finite elimination step, exact affine-polygon coefficients, thin cells and lower-dimensional cases. These are finite cross-checks; the separate-context review is recorded in `sources/review-status.md`.

Pattern equality uses exact masks or exact color strings. The [proved scan rectangle](proof/remarks-and-scope.md#11-complete-line-enumeration), [all quadrant threshold cells](proof/remarks-and-scope.md#12-complete-quadrant-enumeration), and [full-line representation lists](proof/full-support-lines.md) justify global coverage. No fit of a polynomial or period is used as a proof.

`MANIFEST.sha256` fingerprints every package file except itself:

~~~sh
shasum -a 256 -c MANIFEST.sha256
~~~

## Numerical landmarks and limits

On the three-line triangle at $`n=30`$, the window has 496 sites, the explicit lower bound is 4059, CN gives 497, and the exact complexity is 7037. Its eventual form is $`31n^2/4+11n/4+\gamma_{n\bmod M}`$, with no fitted constant table claimed. The distinct Case A quadrant star has finite bound 786 and exact count 1022 on that same window.

The full cross on a half-width rational rectangle has linear coefficients $`3,5/2`$; the sparse cross has $`9/2,4`$. For an irrational width, even the full cross admits no periodic linear coefficient with bounded remainder.

CN supplies the finite-difference and marker mechanisms. Quadratic growth alone was already present in earlier work, including Kari–Szabados; it is not a novelty claim. The precise contribution here is the exact area coefficient, the star application of known counting theory, and the stated direct line formula. The general coefficient supplement applies standard finite elimination and lattice geometry and has passed one fresh-context AI adversarial review. No smallest-period result, all-star constant coefficient on integral polygons, formal certification, comprehensive priority search or stronger general Nivat threshold is claimed.

## AI disclosure

The following disclosure records the source-work roles. Fresh-context GPT-6.1 Sol reviews covered this consolidation and the first supplement (25 October 2026) and the second supplement (26 October 2026); the second supplement was developed by GPT-6 Astra. Claude completed the final assembly and editorial review.

GPT-6.1 Sol, in OpenAI Codex under human direction, developed the proof, ran the computations, and drafted this text; separate GPT-6.1 Sol sessions in fresh contexts performed adversarial reviews. Claude planned the work, checked key steps against the sources, and reviewed and edited the final text. GPT-6 Astra, also in OpenAI Codex, advised on research direction and submission scope. No human expert has certified the work.

## License

Original prose and selection/arrangement: CC BY 4.0. Code in `scripts/`: MIT. External publications are cited, not bundled. See [LICENSE](LICENSE).

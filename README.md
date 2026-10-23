# Exact complexity asymptotics and finite-window bounds for Nivat star configurations

Every star configuration in the specified Convex Nivat contest paper has an exact quadratic complexity coefficient on every fixed compact convex shape with nonempty interior. The coefficient is the area of the union of pairwise intersections of its direction bands; a matching upper bound is part of the proof. The result covers Case A and Case B.

## Results and target

The target is Apex Intelligence, *The Convex Nivat Conjecture: A Complexity Lower Bound for Star Configurations, and a Reduction from Low Convex Complexity to Star Configurations*, contest version dated 12 September 2026. Its PDF SHA-256 is `7fd67831155f4226c010fd6af32de558b76643771e6be12dea4eebafa21745a8`. Citations use printed pages and numbered results. This package supports a breakthrough contribution (submission kind `supplement`) concerning stars.

For a fixed convex body $`P`$, let $`B_i(P)`$ be the full band between its two supporting lines parallel to component direction $`v_i`$, and put

```math
\mathcal M(P)=\bigcup_{i\lt j}(B_i(P)\cap B_j(P)).
```

The three main results are:

- **Theorem 1:** for every CN star and every fixed compact convex $`P`$ with nonempty interior, $`P_\theta(nP\cap\mathbb Z^2)=\mathrm{area}(\mathcal M(P))n^2+O_{\theta,P}(n)`$.
- **Theorem 2:** an explicit finite-window bound for every star and every nonempty finite window, counting globally unique central-marker and distant two-direction anchors. Automatic certificate construction is included. For some small windows the bound is weaker than CN's $`|S|+1`$.
- **Theorem 3:** the binary three-line example on the standard triangle satisfies $`P_c(S_n)=31n^2/4+11n/4+O(1)`$. This is an example-specific linear coefficient.

On the same 496-site triangle $`S_{30}`$, the three-line finite bound is 4059, CN's Theorem T gives 497, and the exact complexity is 7037. For the Case A quadrant star the respective numbers are 786, 497 and 1022.

## Proof guide

Read [objects and statements](proof/objects.md) first. The remaining modules are:

1. [Signals and anchors](proof/signals-and-anchors.md): nonzero Case B strip signals and recovery from two extreme rows; includes the sharper uncorrected bound available in Case B.
2. [Leading-term geometry](proof/leading-term.md): fixed erosion losses for arbitrary convex bodies, lattice counting, the Case B lower bound, the matching upper bound for every star, and exact area evaluation by finite inclusion-exclusion.
3. [All-star extension](proof/all-star-extension.md): common end extrema, bounded cores, distant-anchor uniqueness and complete finite-marker observation; completes Theorem 1.
4. [Finite-window bound](proof/finite-window-bound.md): Theorem 2, explicit certificate construction and the Case A numerical example.
5. [Three-line second order](proof/three-line-second-order.md): Theorem 3 and its collision classification.
6. [Remarks and scope](proof/remarks-and-scope.md): the 4059 calculation, complete line-rectangle and quadrant-cell coverage proofs, classical sharpness, attribution and limitations.

The [bibliography](sources/bibliography.md) identifies the exact versions and checked source locations. The convex-hull observation-center issue, including an example where the origin exponent cancels, is addressed explicitly in the objects, all-star and finite-window modules.

## Reproducing the computations

Python 3.9 or later and its standard library suffice. No downloads, account, environment file or third-party packages are required. From the repository root run:

```sh
python3 -B scripts/verify.py
```

The runner executes all four evidence programs and compares the generated JSON with the checked-in deterministic outputs. Each evidence program writes only its adjacent JSON file. It can also be run separately with `python3 -B scripts/NAME.py`:

- `area_and_anchor_audit.py`: exact rational integration of band areas, the 4059 count and global line-pattern enumeration. Output: `area_and_anchor_audit.json`.
- `end_profile_audit.py`: 24 direction cases across six four-direction stars, including unequal periodic tails, nonlinear encoding and cancellation of the origin exponent. Output: `end_profile_audit.json`.
- `three_line_enumeration.py`: exact global enumeration for every integer $`3\le n\le164`$, with representation multiplicities and the linear-size collision family. Output: `three_line_enumeration.json`.
- `finite_window_audit.py`: all 511 nonempty subsets of a three-by-three window, 120 seeded additional windows, the quadrant examples and 16 general certificate trials. Output: `finite_window_audit.json`.

Pattern equality uses exact integer masks. The line enumeration scans [the proved sufficient rectangle](proof/remarks-and-scope.md#11-complete-line-enumeration) for every multiple-hit anchor, then separately realizes every pure section and zero. The quadrant check covers [all threshold cells, including unbounded ones](proof/remarks-and-scope.md#12-complete-quadrant-enumeration). Finite experiments are evidence about instances, not a replacement for the proofs. Runtime depends on the machine; the larger line enumeration is the dominant computation.

`MANIFEST.sha256` fingerprints all package files except itself. For example, it can be verified using `shasum -a 256 -c MANIFEST.sha256`.

## Attribution and limits

CN supplies the star structure, finite-difference geometry, nonzero-strip argument and shifted-strip product marker (Lemma 4.1, p. 13; Lemma 4.2, pp. 13–14). The proof identifies how extreme rows and complete marker observations turn these into globally unique patterns and an exact area coefficient. It does not use the Kari–Szabados ideal-structure theorem, a spectral polynomial or a scalar ideal bridge.

Quadratic growth alone is not claimed as new. Kari–Szabados already prove rectangular area-order lower bounds and obtain quadratic increments in parts of that proof. The quantitative contribution here is the exact coefficient with a matching upper bound for all stars on arbitrary fixed convex shapes. No comprehensive literature search has been carried out. The single-defect threshold example is presented only as a classical sharpness remark, with attribution through Szabados's introduction and its reference to Cassaigne.

No stronger periodicity threshold for arbitrary configurations, general second-order formula, exact eventual quasipolynomial, or formal proof-assistant certification is claimed.

## AI disclosure

GPT-6.1 Sol, in OpenAI Codex under human direction, developed the proof, ran the computations, and drafted this text; separate GPT-6.1 Sol sessions in fresh contexts performed adversarial reviews. Claude planned the work, checked key steps against the sources, and reviewed and edited the final text. GPT-6 Astra, also in OpenAI Codex, advised on research direction and submission scope. No human expert has certified the work.

## License

Original prose and selection/arrangement: CC BY 4.0. Code in `scripts/`: MIT. External publications retain their own rights and are cited, not bundled. See [LICENSE](LICENSE).

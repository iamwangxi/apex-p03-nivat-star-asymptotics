# Bibliography and source checks

All source checks for this package used archived PDF/text copies. No online search or comprehensive priority review was performed. Page numbers below are printed pages of the identified versions, not extraction line numbers. External publications are not bundled.

## CN: target contest paper

Apex Intelligence. *The Convex Nivat Conjecture: A Complexity Lower Bound for Star Configurations, and a Reduction from Low Convex Complexity to Star Configurations*. Contest PDF, 12 September 2026, 37 pages.

PDF SHA-256: `7fd67831155f4226c010fd6af32de558b76643771e6be12dea4eebafa21745a8`.

Public source: [competition-distributed PDF](https://math.apexin.net/papers/convex-nivat.pdf), dated 12 September 2026. The SHA-256 above identifies the checked artifact; no later version is substituted.

Locations checked:

- §0.1–0.2, p. 4: star hypotheses, complexity and Theorem T.
- §0.3, pp. 4–5: characteristic-zero differences, common tail lattice, periods and formal exponent sets.
- Lemma 1.1, pp. 5–6: finite-difference geometry for local functions; §1, p. 6: the Case A/B split.
- §2.1, p. 7: the two isolating limits and pure-tail fields.
- §2.2, equation (2.3), p. 8: Newton-polygon additivity.
- Lemma 4.1, p. 13: nonzero strip signals in Case B.
- Lemma 4.2, pp. 13–14: product of shifted strip signals with finite nonempty support. The product construction itself is inherited, not claimed as new here.
- Theorem 7.3, p. 18 (Theorem T): $`P_\theta(S)\ge|S|+1`$ for a nonempty finite lattice-convex window. This is the numerical comparison baseline. The [nonconvex connected-region corollary](../proof/nonconvex-windows.md) also uses Theorem T, only for the remaining aligned-parallelogram case; its Theorems 1–4 do not use it.

## Kari–Szabados

Jarkko Kari and Michal Szabados. *An Algebraic Geometric Approach to Nivat's Conjecture*. [arXiv:1605.05929v1](https://arxiv.org/abs/1605.05929v1), deposited 19 May 2016; manuscript dated 20 May 2016. Later journal citation: *Information and Computation* 271 (2020), 104481. The theorem numbers below refer specifically to the checked arXiv v1, not an unchecked journal renumbering.

Acquisition: an archived arXiv v1 PDF and its text extraction.

Locations checked: Lemmas 38–40, pp. 25–27; Lemma 42, pp. 28–29; Theorem 44 and Corollary 45, p. 30. Theorem 44 says a nonperiodic configuration has $`P_c(M,N)>MN`$ outside finitely many rectangle sizes; parts of the supporting proof already exhibit quadratic increments. These facts are used only to delimit novelty, not as premises of the proofs in this package.

## Szabados and the classical sharpness attribution

Michal Szabados. *Nivat's conjecture holds for sums of two periodic configurations*. [arXiv:1710.05360v1](https://arxiv.org/abs/1710.05360v1), 15 October 2017. Conference version: *SOFSEM 2018*, LNCS 10706, pp. 539–551. The source checked here is the arXiv v1.

Acquisition: an archived arXiv v1 PDF and its text extraction.

Locations checked: Introduction, p. 2, recalling nonperiodic configurations with rectangular complexity $`MN+1`$; reference Cas99, p. 13. These support the historical placement of the single-defect observation as classical sharpness. No theorem of this paper is used as a proof premise.

Julien Cassaigne. *Double sequences with complexity mn+1*. *Journal of Automata, Languages and Combinatorics* 4(3), 153–170 (1999). Bibliographic details and the classification attribution were checked through Szabados, p. 2 and reference Cas99. Cassaigne's original article was not available in the source set and was not directly checked; no claim about an individual classification entry is made.

## Bogart–Goodrick–Woods: counting theorem

Tristram Bogart, John Goodrick and Kevin Woods. *Parametric Presburger Arithmetic: Logic, Combinatorics, and Quasi-polynomial Behavior*. *Discrete Analysis* 2017:4, 34 pages. [arXiv:1608.08520v2](https://arxiv.org/abs/1608.08520v2), 16 January 2017; [DOI 10.19086/da.1254](https://doi.org/10.19086/da.1254).

Checked in the archived source: Definitions 1.1 and 1.3, p. 2; Property 2, p. 5; Theorem 1.15, p. 6; Definitions 3.1 and 3.3, p. 14. The theorem allows integer-coordinate fibers, quantifiers and Boolean combinations; the nonnegative-coordinate restriction belongs only to its separate generating-function property. The present counting application independently proves that every representative fiber is finite.

PDF SHA-256: `e516f97579717ba0ec11551e0b0b06ea79b827165d5b6cebc2ffc1fc6b586f39`.
Text SHA-256: `e8badbabe2aafff01681c3d363e70142c09f525a53479d60cf17d5578157c3eb`.

## Woods: earlier fixed-normal context

Kevin Woods. *The unreasonable ubiquitousness of quasi-polynomials*. [arXiv:1308.4694v2](https://arxiv.org/abs/1308.4694v2), 2 March 2014. Theorem 1.25, p. 9, and Theorem 3.5(b), p. 15, cover relevant fixed-normal Presburger counting. These older results already suffice for rational polygon dilations after a definable representative choice; the later Bogart–Goodrick–Woods theorem is used for a uniform counting interface. No claim is made that the general quantified counting principle is new.

## Classical ingredients in the supplementary construction

Pick's theorem and finite residue-and-bound quantifier elimination are classical. The supplementary page proves the precise elimination step, cell decomposition, thin-cell estimate and vertical-slice period bound it uses. It does not appeal to an unchecked numerical output of a quantifier-elimination library. That application passed one fresh-context GPT-6.1 Sol adversarial review on 25 October 2026; see [review status](review-status.md) for its scope and the distinction from human certification or proof-assistant verification.

## Ingredients of the nonconvex and integral-window extensions

The [nonconvex extension](../proof/nonconvex-windows.md) reuses CN’s star definition and strip-product marker and proves its interval-union geometry, equality rigidity and graph collision formula in full. The [integral-window extension](../proof/integral-windows.md) uses the established area/quasipolynomial results and classical Pick and tangent-cone arguments; its original-table boundary classification is proved on that page. No additional external theorem or unchecked new publication is invoked. The exact period-four example uses the package’s [cone and period constructions](../proof/minimal-periods-and-constant-term.md).

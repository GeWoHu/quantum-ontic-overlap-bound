# REVIEW REQUEST — we ask experts to try to break this

**PREPRINT FOR EXPERT REVIEW — NOT PEER REVIEWED.**
We invite hostile, specific criticism. Please open an issue with the
front number. Four fronts, in decreasing order of how much sleep this
costs us:

## Front 1 — Stratification completeness
The classification routes every endpoint pattern through exactly one of:
generic strata, eight coincidence lines/edges, midpoint, corners. Frozen
zero-sets are chosen by exact rational probe samples (Remark 6.1 of the
manuscript). Questions:
- Is the containment argument airtight — at every (p, q) ∈ [0,1]², is the
  true zero-set really contained in the frozen support of exactly one
  enumerated stratum, so that support-restricted infeasibility excludes all
  models outside Case A?
- Could a model exist whose true zero-sets are strictly larger than every
  frozen support considered, without lying on an enumerated locus?
- Is the probe grid (fixed in `certificates/export_certificates.py`)
  sufficient, or can you construct a pattern whose generic zero-structure
  is missed?

## Front 2 — Algebraic certificates
Each of 844 excluded strata is certified by 1 ∈ I for the ideal I of the
support-restricted Born + simplex polynomial system over ℚ (Hilbert's
Nullstellensatz ⇒ no complex solution). Questions:
- Recompute from `certificates/CERT_LEDGER.json` (run
  `python3 certificates/verify_certificates.py`). Any mismatch?
- Does the elimination order (response variables, then parameters, then
  the dummy y) correctly encode "no solution in any variable"?
- Independent recomputation with Macaulay2/Singular/Magma on a sample: any
  disagreement with sympy's reduced Gröbner bases?
- Is the relaxation direction (larger frozen support ⇒ weaker exclusion,
  never stronger) correctly applied in every stratum kind?

## Front 3 — Born-rule consistency
- Verify the four models in Appendix B against all 16 Born equations + 4
  simplex normalizations with exact rationals (minutes by hand).
- The Endpoint Lemma derives elementwise complementarity
  ξ_x(λ) π_x(λ) = 0 from the zero diagonal of B. Does the derivation use
  only the frozen convention (outcome index = preparation index), and is
  the label-covariance claim under μ₁ ↔ μ₂ correct?
- Marginal consistency: check the four exhibited couplings have the stated
  marginals and endpoint parameters (Appendix B lists all t_x values).
- Is the mixed-side unconditional inequality (proof of Theorem 1) valid at
  the region boundaries p = 1/2 or q = 1/2?

## Front 4 — Novelty
Our claim is deliberately narrow (Section 7 + `novelty` checks below).
- Does any published result imply ω² ≤ ½ Σ_x δ_x for k = 2 with
  **unrestricted couplings** (no preparation independence)? Closest known:
  Branciard PRL 113, 020409 (arXiv:1407.3005); Kunjwal–Spekkens PRL 115,
  110403 (arXiv:1506.04150); Mansfield PRA 94, 042124 (arXiv:1412.0669);
  Aaronson et al. PRA 88, 032111 (arXiv:1303.2834); Abramsky et al. PRL
  119, 050504 (arXiv:1705.07918). All IDs arXiv-API title-matched.
- Is the exhaustive four-model rigidity classification for this exclusion
  table at k = 2 known in the literature (e.g., as a special case of
  finite-ontic-size analyses)?
- Is the constant 1/2, or the overlap-vs-disturbance form, anticipated in
  the error–disturbance or preparation-contextuality literature?

## What would change our verdict
A single exact counterexample model (rational ξ, π, p, q satisfying all 16
Born equations + simplex + marginals with ω² > ½Σδ) refutes the headline.
A missed stratum with a feasible model refutes the classification. Both
would be acknowledged prominently and the preprint revised.

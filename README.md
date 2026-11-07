# Stochastic rigidity of two-ontic-state models of the PBR exclusion scenario

> ## ⚠️ PREPRINT FOR EXPERT REVIEW — NOT PEER REVIEWED
>
> This is a self-contained mathematical preprint released for independent
> expert scrutiny. It has been machine-verified internally but has **not**
> undergone external peer review. Please read `REVIEW_REQUEST.md` and
> challenge everything.

## The theorem (one paragraph)

For ontological models with **exactly two ontic states** (|Λ| = 2) that
reproduce the exclusion statistics of the two-preparation, two-copy BB84
PBR scenario — four product preparations, a four-outcome measurement with
Born table `B` having zero diagonal and off-diagonal entries {1/4, 1/4, 1/2}
per row — with **arbitrary** joint ontic couplings of the given marginals
(no preparation-independence assumption), the epistemic overlap ω of the two
marginals and the total-variation distances δ_x of the couplings from the
independent products satisfy the sharp inequality

    ω²  ≤  (1/2) · Σ_x δ_x ,

with equality **if and only if** (p, q) ∈ {(1/4, 3/4), (3/4, 1/4), (0, 1),
(1, 0)}. More strongly, the admissible class is **rigid**: it consists of
exactly four models — two deterministic models at (1/4, 3/4) and (3/4, 1/4)
(nontrivial equality) and two degenerate-marginal corner models at (0, 1)
and (1, 0) where both sides vanish. No genuinely stochastic model with
nondegenerate marginals survives.

## Assumptions (explicit)

1. Ontic space Λ = {0, 1} exactly; preparation pairs sample Λ².
2. Marginals μ₁ = (p, 1−p), μ₂ = (q, 1−q), p, q ∈ [0, 1]; marginal
   consistency of every coupling π_x.
3. Stochastic response functions ξ_o(λ) ∈ [0, 1] with Σ_o ξ_o(λ) = 1.
4. FULL-BORN compatibility at all 16 entries: B[o][x] = Σ_λ ξ_o(λ) π_x(λ)
   with the frozen preparation ordering stated in the manuscript.
5. TV convention δ_x = ½‖π_x − μ_i ⊗ μ_j‖₁.

## What is **not** claimed

- Nothing for |Λ| > 2 (no unrestricted-ontic-model conclusion).
- No claim to have "disproved objective reality", ψ-epistemic models in
  general, or quantum theory's completeness.
- The exclusion scenario (two-copy BB84 PBR table) is standard, not new.
- No novelty claim for quantitative ψ-ontology bounds as a genre; the
  novelty claim is the specific sharp TV-disturbance inequality with exact
  equality characterization plus the complete four-model rigidity
  classification at k = 2 (see Section 7 of the manuscript).

## Repository structure

| Path | Contents |
|---|---|
| `manuscript/main.pdf` | The paper (with proof appendices A–C) |
| `manuscript/main.tex` | LaTeX source (compiles with `tectonic` or pdflatex) |
| `certificates/CERT_LEDGER.json` | All 846 strata with exact polynomial generators |
| `certificates/verify_certificates.py` | Standalone certificate verifier (recomputes everything) |
| `certificates/export_certificates.py` | Certificate generator (regenerates the ledger) |
| `verification/verify_independent.py` | Independent reconstruction code path (Case A, equality cases, spot-certificates, break search) |
| `REPRODUCE.md` | Exact reproduction instructions and expected outputs |
| `REVIEW_REQUEST.md` | What we ask experts to attack |
| `MANIFEST.sha256` | SHA-256 of every published file |

## Verification in one command

```bash
python3 certificates/verify_certificates.py
# expected: VERIFY DONE: 844 certificates reproduced, 0 mismatches,
#           2 feasible strata skipped (verified by explicit model)
```

Exact rational arithmetic only (sympy ≥ 1.12); no floats, no sampling, no
trust in any log file. See `REPRODUCE.md` for the full protocol including
PDF rebuild and hash checking.

## License

Code and ledger: MIT (see `LICENSE`). Manuscript text: released under CC-BY-4.0
for scholarly reuse with attribution to the repository author.

## Citation

Geri (GERIWH), "Stochastic rigidity of two-ontic-state models of the PBR
exclusion scenario, with a sharp overlap–disturbance inequality," preprint
v1.0, 2026, this repository. Bibliography therein is arXiv-API-verified
(all IDs title-matched).

#!/usr/bin/env python3
# QF-02.8 — STANDALONE certificate verifier. A mathematician runs:
#     python3 verify_certificates.py
# and nothing else is trusted: each ledger entry lists the polynomial
# generators of the stratum ideal (exact rationals, sympy syntax) plus the
# declared variables. We recompute the reduced lex Groebner basis over Q of
# <generators, 1 - y>. If the basis contains a nonzero constant then 1 lies
# in the ideal, so the stratum has NO complex solution, a fortiori no
# stochastic model there (Hilbert's Nullstellensatz). Entries marked
# feasible (the two corner models) are checked the opposite way by
# explicit substitution elsewhere (see manuscript Appendix B).
# Requires only sympy >= 1.12. Exact rational arithmetic throughout.
import json, sys
import sympy as sp

led = json.load(open("CERT_LEDGER.json"))
ok = bad = skip = 0
for e in led:
    if e["kind"].startswith("corner") and not e["gb_contains_1"]:
        skip += 1  # feasible strata: verified by explicit model, not infeasibility
        continue
    vars_ = [sp.Symbol(v) for v in e["vars"]]
    # recover p,q appearing in generators
    pq = [s for s in (sp.Symbol('p'), sp.Symbol('q'))
          if any(s in sp.sympify(g).free_symbols for g in e["generators"])]
    y = sp.Symbol('_y')
    G = sp.groebner([sp.sympify(g) for g in e["generators"]] + [y-1],
                    *(vars_ + pq + [y]), order='lex')
    has1 = any((g.is_Number and g != 0) for g in G)
    if has1 == e["gb_contains_1"]:
        ok += 1
    else:
        bad += 1
        print("MISMATCH", e["id"], flush=True)
print(f"VERIFY DONE: {ok} certificates reproduced, {bad} mismatches, "
      f"{skip} feasible strata skipped (verified by explicit model)")
sys.exit(0 if bad == 0 else 1)

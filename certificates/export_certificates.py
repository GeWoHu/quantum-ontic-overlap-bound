#!/usr/bin/env python3
# QF-02.8 — export exact Gröbner certificates for EVERY excluded stratum of the
# k=2 FULL-BORN classification, as a machine-readable ledger:
#   {id, kind, pattern, condition, support, vars, generators, gb_constant}
# A verifier recomputes the reduced lex Groebner basis of `generators` over Q
# and checks it contains the constant 1 (=> empty variety => no model there).
# Exact rational arithmetic only; sympy >= 1.12.
import sympy as sp, itertools, json, sys

p,q = sp.symbols('p q')
B = [[0,sp.Rational(1,4),sp.Rational(1,4),sp.Rational(1,2)],
     [sp.Rational(1,4),0,sp.Rational(1,2),sp.Rational(1,4)],
     [sp.Rational(1,4),sp.Rational(1,2),0,sp.Rational(1,4)],
     [sp.Rational(1,2),sp.Rational(1,4),sp.Rational(1,4),0]]
preps = [(0,1),(0,0),(1,1),(1,0)]
SGN = [1,-1,-1,1]

def prods(mu):
    out=[]
    for (i,j) in preps:
        out.append([mu[i][0]*mu[j][0], mu[i][0]*mu[j][1],
                    mu[i][1]*mu[j][0], mu[i][1]*mu[j][1]])
    return out

def endpoint(pr, ch):
    return [-pr[0], pr[1], pr[2], -pr[3]][ch]

def pi_of(pr, ch):
    t = endpoint(pr, ch)
    return [sp.expand(pr[k] + SGN[k]*t) for k in range(4)]

SAMPLES = [(sp.Rational(1,5),sp.Rational(2,5)),(sp.Rational(3,10),sp.Rational(7,10)),
           (sp.Rational(7,10),sp.Rational(3,10)),(sp.Rational(3,5),sp.Rational(1,5)),
           (sp.Rational(1,10),sp.Rational(9,10)),(sp.Rational(9,10),sp.Rational(1,10)),
           (sp.Rational(2,5),sp.Rational(2,5)),(sp.Rational(1,4),sp.Rational(3,4))]

def stratum_generators(P, Zs, unknown_params):
    """P: list of 4 pi-vectors (polys); Zs: support per prep; unknown_params: [p,q] subset."""
    unk=[]; idx={}
    for o in range(4):
        for lam in Zs[o]:
            idx[(o,lam)] = len(unk); unk.append(sp.Symbol(f'x{o}_{lam}'))
    eqs=[]
    for x in range(4):
        for o in range(4):
            eqs.append(sp.expand(sum(unk[idx[(o,lam)]]*P[x][lam] for lam in Zs[o]) - B[o][x]))
    for lam in range(4):
        s = sum(unk[idx[(o,lam)]] for o in range(4) if (o,lam) in idx)
        if s != 0: eqs.append(sp.expand(s-1))
    return unk, eqs

def certify(unk, eqs, params):
    varorder = list(unk) + list(params)
    y = sp.Symbol('_y')
    G = sp.groebner([sp.expand(e) for e in eqs] + [y-1], *(varorder + [y]), order='lex')
    has1 = any((g.is_Number and g != 0) for g in G)
    return has1

ledger = []

def probe_and_build(pattern, subs, cond, kind, params_after):
    """Probe pattern under substitution; return ledger entry or None (dropped/permutation)."""
    ok = False
    for (pv,qv) in SAMPLES:
        good = True; P=[]; Zs=[]
        pr_all = prods([[pv,1-pv],[qv,1-qv]])
        for x in range(4):
            pis = pi_of(pr_all[x], pattern[x])
            if any(v < 0 for v in pis): good=False; break
            Z = [k for k in range(4) if pis[k]==0]
            if not Z: good=False; break
            P.append(pis); Zs.append(Z)
        if good: ok=True; break
    if not ok: return None
    if set().union(*Zs) != {0,1,2,3}: return None
    if all(len(Z)==1 for Z in Zs): return None  # Case A handled analytically
    # symbolic pi under substitution
    mu = [[p,1-p],[q,1-q]]
    Ps=[]
    for x in range(4):
        pr = prods(mu)[x]
        Ps.append([sp.expand(z.subs(subs)) for z in pi_of(pr, pattern[x])])
    unk, eqs = stratum_generators(Ps, Zs, params_after)
    has1 = certify(unk, eqs, params_after)
    return {"id": f"{kind}_{pattern}", "kind": kind,
            "pattern": list(pattern), "condition": cond,
            "support": Zs, "vars": [str(u) for u in unk],
            "param_condition": str(subs),
            "generators": [str(sp.expand(e)) for e in eqs],
            "gb_contains_1": bool(has1)}

# --- generic strata (no substitution) ---
for pattern in itertools.product(range(4), repeat=4):
    e = probe_and_build(pattern, {}, "generic", "generic", [p,q])
    if e: ledger.append(e); print("generic", pattern, e["gb_contains_1"], flush=True)

# --- coincidence lines and degenerate edges ---
LINES = [({q:p}, "p=q", [p]), ({q:1-p}, "p+q=1", [p]),
         ({p:sp.Rational(1,2)}, "p=1/2", [q]), ({q:sp.Rational(1,2)}, "q=1/2", [p]),
         ({p:0}, "p=0", [q]), ({p:1}, "p=1", [q]),
         ({q:0}, "q=0", [p]), ({q:1}, "q=1", [p])]
for subs, cond, pa in LINES:
    for pattern in itertools.product(range(4), repeat=4):
        e = probe_and_build(pattern, subs, cond, "line_"+cond.replace("=","").replace("/",""), pa)
        if e: ledger.append(e); print(cond, pattern, e["gb_contains_1"], flush=True)

# --- midpoint (1/2,1/2): all 16 collapsed endpoint couplings, linear systems ---
mu = [[sp.Rational(1,2),sp.Rational(1,2)],[sp.Rational(1,2),sp.Rational(1,2)]]
pr0 = prods(mu)
seen = set()
for combo in itertools.product([0,1], repeat=4):  # t=-1/4 (ch0) or +1/4 (ch1)
    P = [pi_of(pr0[x], combo[x]) for x in range(4)]
    key = tuple(tuple(pi) for pi in P)
    if key in seen: continue
    seen.add(key)
    Zs = [[k for k in range(4) if P[x][k]==0] for x in range(4)]
    if any(not Z for Z in Zs): continue
    if set().union(*Zs) != {0,1,2,3}: continue
    unk, eqs = stratum_generators(P, Zs, [])
    has1 = certify(unk, eqs, [])
    ledger.append({"id": f"midpoint_{combo}", "kind": "midpoint",
                   "pattern": list(combo), "condition": "p=q=1/2",
                   "support": Zs, "vars": [str(u) for u in unk],
                   "param_condition": "{p:1/2,q:1/2}",
                   "generators": [str(sp.expand(e)) for e in eqs],
                   "gb_contains_1": bool(has1)})
    print("midpoint", combo, has1, flush=True)

# --- corners: full (unrestricted-support) linear systems ---
for (pv,qv) in [(0,0),(1,1),(0,1),(1,0)]:
    mu = [[sp.Integer(pv),1-sp.Integer(pv)],[sp.Integer(qv),1-sp.Integer(qv)]]
    P = prods(mu)  # point-mass marginals force pi=prod
    Zs = [[0,1,2,3] for _ in range(4)]  # full support: no endpoint restriction
    unk, eqs = stratum_generators(P, Zs, [])
    has1 = certify(unk, eqs, [])
    ledger.append({"id": f"corner_{pv}{qv}", "kind": "corner",
                   "pattern": None, "condition": f"(p,q)=({pv},{qv})",
                   "support": Zs, "vars": [str(u) for u in unk],
                   "param_condition": f"{{p:{pv},q:{qv}}}",
                   "generators": [str(sp.expand(e)) for e in eqs],
                   "gb_contains_1": bool(has1)})
    print("corner", (pv,qv), has1, flush=True)

json.dump(ledger, open("CERT_LEDGER.json","w"))
n = len(ledger); c = sum(1 for e in ledger if e["gb_contains_1"])
print(f"=== EXPORT DONE: {n} strata, {c} certified 1-in-ideal", flush=True)

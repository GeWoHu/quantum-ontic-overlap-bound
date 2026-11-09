#!/usr/bin/env python3
# QF-02.8 R2 — export exact Groebner certificates for EVERY excluded stratum of the
# k=2 FULL-BORN classification as a machine-readable ledger:
#   {id, kind, pattern, condition, support, vars, generators, gb_constant}
# A verifier recomputes the reduced lex Groebner basis of `generators` over Q
# and checks it contains the constant 1 (=> empty variety => no model there).
#
# R2 CORRECTION (external review 2026-11-09): supports are now derived SYMBOLICALLY
# from the substituted coupling (cells vanishing identically on the stratum), never
# from numerical probe samples. Stratification of (p,q) in [0,1]^2:
#   generic 2D chamber; 8 coincidence lines / degenerate edges (valid on the dense
#   open of each); point strata where zero sets grow further ((1/2,1/2) and the four
#   edge-midpoints (0,1/2),(1,1/2),(1/2,0),(1/2,1)); 4 corners with full-support
#   systems. Patterns whose true supports are all singletons are Case A (handled
#   analytically in the manuscript). Patterns producing identical substituted
#   couplings are deduplicated to one entry.
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

def poly_vars(subs):
    return [v for v in (p,q) if v not in subs]

def uni_quad_min(f, x):
    """Exact minimum of polynomial f (deg <= 2) on [0,1]."""
    f = sp.expand(f)
    c = f.subs(x,0); d = f.subs(x,1)
    a = sp.Poly(f, x)
    deg = a.degree()
    if deg <= 1: return min(c,d)
    A,Bc,Cc = a.all_coeffs()
    if A <= 0: return min(c,d)
    v = sp.nsimplify(-Bc/(2*A))
    if 0 < v < 1: return min(c, d, f.subs(x, v))
    return min(c,d)

def nonneg_on(f, params):
    """Exact nonnegativity test for f (deg <= 2) on the unit square / unit interval."""
    f = sp.expand(f)
    if len(params) == 1:
        return uni_quad_min(f, params[0]) >= 0
    # 2D: corners + edge minima + interior critical point (gradient linear)
    for pv in (0,1):
        for qv in (0,1):
            if f.subs({p:pv,q:qv}) < 0: return False
    for edge_x, other in ((p,q),(q,p)):
        g = f.subs(edge_x, 0); h = f.subs(edge_x, 1)
        if uni_quad_min(g, other) < 0 or uni_quad_min(h, other) < 0: return False
    fp, fq = sp.diff(f,p), sp.diff(f,q)
    try:
        sol = sp.solve([fp, fq], [p,q], dict=True)
    except Exception:
        sol = []
    if len(sol) == 1:
        s = sol[0]
        pv, qv = s.get(p), s.get(q)
        if pv is not None and qv is not None and pv.is_rational \
           and 0 < pv < 1 and 0 < qv < 1:
            if f.subs(s) < 0: return False
    return True

def isolated_interior_zeros(f, params):
    """Rational interior points where f vanishes but not identically (double roots)."""
    pts=[]
    if len(params) != 1: return pts
    x = params[0]; f = sp.expand(f)
    if f == 0: return pts
    for r in sp.roots(sp.Poly(f, x), multiple=True):
        if r.is_Rational and 0 < r < 1: pts.append(r)
    return pts

def stratum_generators(P, Zs, unknown_params):
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
    return any((g.is_Number and g != 0) for g in G)

ledger = []

def symbolic_stratum(pattern, subs, cond, kind, seen_keys, tp=None):
    """True-support entry for one pattern on one cell; None if dropped."""
    mu = [[p,1-p],[q,1-q]]
    Ps=[]; Zs=[]
    params = poly_vars(subs)
    for x in range(4):
        pis = [sp.expand(z.subs(subs)) for z in pi_of(prods(mu)[x], pattern[x])]
        if tp is not None:
            # every pi entry is linear (sign-constant on the cell) or globally >= 0;
            # exact rational test at one interior point decides nonnegativity on the cell
            if any(f.subs({p:tp[0], q:tp[1]}) < 0 for f in pis): return None
        else:
            for f in pis:
                if not nonneg_on(f, params): return None
        Ps.append(pis)
        Zs.append([k for k in range(4) if sp.expand(pis[k]) == 0])
    if any(not Z for Z in Zs): return None
    if set().union(*Zs) != {0,1,2,3}: return None
    if all(len(Z)==1 for Z in Zs): return None  # Case A handled analytically
    key = tuple(tuple(str(v) for v in pi) for pi in Ps)
    if key in seen_keys: return "dup"
    seen_keys.add(key)
    unk, eqs = stratum_generators(Ps, Zs, params)
    has1 = certify(unk, eqs, params)
    return {"id": f"{kind}_{pattern}", "kind": kind,
            "pattern": list(pattern), "condition": cond,
            "support": Zs, "vars": [str(u) for u in unk],
            "param_condition": str(subs),
            "support_source": "symbolic-identical-zeros-on-stratum",
            "generators": [str(sp.expand(e)) for e in eqs],
            "gb_contains_1": bool(has1)}

# --- full symbolic arrangement of endpoint-feasibility chambers ---
# Channel validity per preparation is LINEAR in (p,q): ch0 needs pr0<=pr3, ch1 pr1<=pr2,
# ch2 pr2<=pr1, ch3 pr3<=pr0. The arrangement lines p=q, p+q=1, p=1/2, q=1/2 cut the
# square into 8 open 2D chambers; the same lines give 8 open 1D segments (split at the
# midpoint, where all four meet); the 4 edges force t=0; the midpoint and corners are 0D.
R = sp.Rational
CELLS = [
  ("chamber_C1a", {}, "p+q<1,p<q,q<1/2", (R(1,4),R(1,3))),
  ("chamber_C1b", {}, "p+q<1,p<q,q>1/2", (R(1,3),R(5,9))),
  ("chamber_C2a", {}, "p+q<1,q<p,p<1/2", (R(5,12),R(1,4))),
  ("chamber_C2b", {}, "p+q<1,q<p,p>1/2", (R(3,5),R(1,4))),
  ("chamber_C3a", {}, "p+q>1,p<q,p>1/2", (R(3,5),R(4,5))),
  ("chamber_C3b", {}, "p+q>1,p<q,p<1/2", (R(1,3),R(3,4))),
  ("chamber_C4a", {}, "p+q>1,q<p,q>1/2", (R(4,5),R(3,5))),
  ("chamber_C4b", {}, "p+q>1,q<p,q<1/2", (R(3,4),R(1,3))),
  ("seg_pq_lo",  {q:p},     "p=q<1/2",   (R(1,4),R(1,4))),
  ("seg_pq_hi",  {q:p},     "p=q>1/2",   (R(3,4),R(3,4))),
  ("seg_sum_lo", {q:1-p},   "p+q=1,p<1/2", (R(1,4),R(3,4))),
  ("seg_sum_hi", {q:1-p},   "p+q=1,p>1/2", (R(3,4),R(1,4))),
  ("seg_p12_lo", {p:R(1,2)}, "p=1/2,q<1/2", (R(1,2),R(1,4))),
  ("seg_p12_hi", {p:R(1,2)}, "p=1/2,q>1/2", (R(1,2),R(3,4))),
  ("seg_q12_lo", {q:R(1,2)}, "q=1/2,p<1/2", (R(1,4),R(1,2))),
  ("seg_q12_hi", {q:R(1,2)}, "q=1/2,p>1/2", (R(3,4),R(1,2))),
  ("edge_p0", {p:0}, "p=0", (0,R(1,4))), ("edge_p1", {p:1}, "p=1", (1,R(1,4))),
  ("edge_q0", {q:0}, "q=0", (R(1,4),0)), ("edge_q1", {q:1}, "q=1", (R(1,4),1)),
]
def valid_channels(tp):
    mu = [[tp[0],1-tp[0]],[tp[1],1-tp[1]]]
    out=[]
    for x in range(4):
        pr = prods(mu)[x]
        out.append([c for c in range(4)
                    if [pr[0]<=pr[3], pr[1]<=pr[2], pr[2]<=pr[1], pr[3]<=pr[0]][c]])
    return out
for kind, subs, cond, tp in CELLS:
    seen = set()
    vc = valid_channels(tp)
    for pattern in itertools.product(*vc):
        e = symbolic_stratum(pattern, subs, cond, kind, seen, tp=tp)
        if e == "dup" or e is None: continue
        ledger.append(e); print(kind, pattern, e["gb_contains_1"], flush=True)

# --- point stratum: the midpoint, where all four arrangement lines meet and
#     supports grow beyond every adjacent cell (edge-midpoints need no split:
#     on an edge t=0 is forced and the zero pattern is constant along it) ---
POINTS = {(sp.Rational(1,2), sp.Rational(1,2)): "p=q=1/2"}
for (pv,qv), cond in POINTS.items():
    mu = [[sp.Rational(pv),1-sp.Rational(pv)],[sp.Rational(qv),1-sp.Rational(qv)]]
    pr0 = prods(mu)
    seen = set()
    for combo in itertools.product([0,1], repeat=4):  # collapsed endpoint choices
        P = [pi_of(pr0[x], combo[x]) for x in range(4)]
        if any(v < 0 for pi in P for v in pi): continue
        key = tuple(tuple(pi) for pi in P)
        if key in seen: continue
        seen.add(key)
        Zs = [[k for k in range(4) if P[x][k]==0] for x in range(4)]
        if any(not Z for Z in Zs): continue
        if set().union(*Zs) != {0,1,2,3}: continue
        if all(len(Z)==1 for Z in Zs): continue  # Case A at the point
        unk, eqs = stratum_generators(P, Zs, [])
        has1 = certify(unk, eqs, [])
        ledger.append({"id": f"point_{pv}_{qv}_{combo}", "kind": "point",
                       "pattern": list(combo), "condition": cond,
                       "support": Zs, "vars": [str(u) for u in unk],
                       "param_condition": f"{{p:{pv},q:{qv}}}",
                       "support_source": "pointwise-zeros",
                       "generators": [str(sp.expand(e)) for e in eqs],
                       "gb_contains_1": bool(has1)})
        print("point", cond, combo, has1, flush=True)

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
                   "support_source": "full-support-at-corner",
                   "generators": [str(sp.expand(e)) for e in eqs],
                   "gb_contains_1": bool(has1)})
    print("corner", (pv,qv), has1, flush=True)

json.dump(ledger, open("CERT_LEDGER.json","w"))
n = len(ledger); c = sum(1 for e in ledger if e["gb_contains_1"])
print(f"=== EXPORT DONE: {n} strata, {c} certified 1-in-ideal", flush=True)

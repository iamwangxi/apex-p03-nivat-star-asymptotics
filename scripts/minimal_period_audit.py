#!/usr/bin/env python3
"""Exact finite cross-checks; no fitting, network, third-party code or compiler claim.

The coefficient routines use rational vertex cones. The comparison routine uses
independent vertical slices of the original inequalities. The star routine lists
all multi-line anchors and all realizable pure patterns, not a guessed anchor box.
Run with python3 -B scripts/minimal_period_audit.py from the repository root;
JSON is printed to stdout.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, lcm
import json
import random


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def scale(a, x):
    return (a[0] * x, a[1] * x)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def floor(x):
    x = F(x)
    return x.numerator // x.denominator


def ceil(x):
    return -floor(-x)


def primitive(a):
    d = lcm(a[0].denominator, a[1].denominator)
    x, y = int(d * a[0]), int(d * a[1])
    g = gcd(abs(x), abs(y))
    assert g
    return (x // g, y // g)


def intersection(c1, c2):
    a, b, s, t = c1
    c, d, u, v = c2
    D = a * d - b * c
    if not D:
        return None
    return ((F(s*d-b*u, D), F(a*u-s*c, D)),
            (F(t*d-b*v, D), F(a*v-t*c, D)))


def at(v, k):
    return add(scale(v[0], k), v[1])


def root_bound(coeffs):
    """An integer strictly above the moduli of all roots (Cauchy's bound)."""
    coeffs = list(map(F, coeffs))
    while coeffs and not coeffs[-1]:
        coeffs.pop()
    if len(coeffs) <= 1:
        return 1
    return floor(1 + max(abs(x / coeffs[-1]) for x in coeffs[:-1])) + 1


def hull(points):
    ps = sorted(set(points))
    if len(ps) <= 1:
        return ps
    def turn(a, b, c):
        return det(sub(b, a), sub(c, a))
    lo = []
    for p in ps:
        while len(lo) >= 2 and turn(lo[-2], lo[-1], p) <= 0:
            lo.pop()
        lo.append(p)
    hi = []
    for p in reversed(ps):
        while len(hi) >= 2 and turn(hi[-2], hi[-1], p) <= 0:
            hi.pop()
        hi.append(p)
    return lo[:-1] + hi[:-1]


def stable_vertices(cs):
    """All test inputs are bounded; stability is certified by affine signs."""
    candidates = list(dict.fromkeys(v for x, y in combinations(cs, 2)
                                    if (v := intersection(x, y)) is not None))
    accepted = []
    N = 1
    for v in candidates:
        ok = True
        for a, b, s, t in cs:
            slope = a*v[0][0] + b*v[0][1] - s
            const = a*v[1][0] + b*v[1][1] - t
            N = max(N, root_bound([const, slope]))
            if slope > 0 or (not slope and const > 0):
                ok = False
        if ok:
            accepted.append(v)
    for v, w in combinations(accepted, 2):
        for j in (0, 1):
            N = max(N, root_bound([v[1][j]-w[1][j], v[0][j]-w[0][j]]))
    for v, w, z in combinations(accepted, 3):
        a, b = sub(w[0], v[0]), sub(w[1], v[1])
        c, d = sub(z[0], v[0]), sub(z[1], v[1])
        N = max(N, root_bound([det(b, d), det(a, d)+det(b, c), det(a, c)]))
    by_point = {at(v, N): v for v in accepted}
    vs = [by_point[p] for p in hull(by_point)]
    E = lcm(*(x.denominator for v in vs for x in v[0])) if vs else 1
    return vs, N, E


def parallelogram_points(b, w1, w2):
    D = det(w1, w2)
    corners = [b, add(b, w1), add(b, w2), add(add(b, w1), w2)]
    out = []
    for x in range(ceil(min(z[0] for z in corners)), floor(max(z[0] for z in corners))+1):
        for y in range(ceil(min(z[1] for z in corners)), floor(max(z[1] for z in corners))+1):
            z = (x, y)
            u = sub(z, b)
            t1, t2 = F(det(u, w2), D), F(det(w1, u), D)
            if 0 <= t1 < 1 and 0 <= t2 < 1:
                out.append(z)
    assert len(out) == abs(D)
    return out


def bezout(a, b):
    if b == 0:
        return (abs(a), 1 if a > 0 else -1, 0)
    g, x, y = bezout(b, a % b)
    return (g, y, x - (a // b)*y)


def cone_coefficients(vs, N, E, h):
    """A, B, C for k=E*t+h, including empty, point and segment fibers."""
    if not vs:
        return (F(0), F(0), F(0))
    if len(vs) == 1:
        b = at(vs[0], h)
        return (F(0), F(0), F(all(x.denominator == 1 for x in b)))
    if len(vs) == 2:
        v, u = vs
        w = primitive(sub(at(u, N), at(v, N)))
        nu = (-w[1], w[0])
        b0, b1 = at(v, h), at(u, h)
        level = dot(nu, b0)
        if level.denominator != 1:
            return (F(0), F(0), F(0))
        g, x, y = bezout(*nu)
        assert g == 1
        x0 = (x*int(level), y*int(level))
        j = 0 if w[0] else 1
        tau0, tau1 = F(b0[j]-x0[j], w[j]), F(b1[j]-x0[j], w[j])
        L = F(E*(u[0][j]-v[0][j]), w[j])
        assert L.denominator == 1 and L >= 0
        return (F(0), L, F(floor(tau1)-ceil(tau0)+1))
    directions = []
    for j, v in enumerate(vs):
        directions.extend(primitive(sub(at(vs[i], N), at(v, N)))
                          for i in ((j-1) % len(vs), (j+1) % len(vs)))
    T = 1
    while any(dot((1, T), w) == 0 for w in directions):
        T += 1
    c = (1, T)
    A = B = C = F(0)
    for j, v in enumerate(vs):
        w1, w2 = directions[2*j:2*j+2]
        d1, d2 = dot(c, w1), dot(c, w2)
        b = at(v, h)
        eta = dot(c, scale(v[0], E))
        assert eta.denominator == 1
        for z in parallelogram_points(b, w1, w2):
            x = dot(c, z)
            A += F(eta*eta, 2*d1*d2)
            B += F(eta*x, d1*d2) - F(eta, 2)*(F(1,d1)+F(1,d2))
            C += (F(x*x, 2*d1*d2) - F(x, 2)*(F(1,d1)+F(1,d2))
                  + F(1,4) + (F(d1,d2)+F(d2,d1))/12)
    return A, B, C


def count_slices(cs, k):
    ps = []
    for c1, c2 in combinations(cs, 2):
        v = intersection(c1, c2)
        if v is None:
            continue
        p = at(v, k)
        if all(a*p[0]+b*p[1] <= s*k+t for a,b,s,t in cs):
            ps.append(p)
    if not ps:
        return 0
    total = 0
    for x in range(ceil(min(p[0] for p in ps)), floor(max(p[0] for p in ps))+1):
        lower, upper = [], []
        good = True
        for a,b,s,t in cs:
            if b > 0:
                upper.append(F(s*k+t-a*x, b))
            elif b < 0:
                lower.append(F(s*k+t-a*x, b))
            elif a*x > s*k+t:
                good = False
        if good:
            assert lower and upper
            total += max(0, floor(min(upper))-ceil(max(lower))+1)
    return total


def divisors(n):
    return [d for d in range(1, n+1) if n % d == 0]


def trim(p):
    p = list(map(F, p))
    while p and p[-1] == 0:
        p.pop()
    return p


def divrem(p, q):
    p, q = trim(p), trim(q)
    assert q
    out = [F(0)] * max(0, len(p)-len(q)+1)
    while len(p) >= len(q):
        d = len(p)-len(q)
        v = p[-1]/q[-1]
        out[d] = v
        for j, x in enumerate(q):
            p[d+j] -= v*x
        p = trim(p)
    return trim(out), p


PHI = {}


def cyclotomic(d):
    if d not in PHI:
        p = [F(-1)] + [F(0)]*(d-1) + [F(1)]
        for e in divisors(d)[:-1]:
            p, r = divrem(p, cyclotomic(e))
            assert not r
        PHI[d] = p
    return PHI[d]


def period_fourier(word):
    present = [d for d in divisors(len(word)) if divrem(word, cyclotomic(d))[1]]
    return lcm(*present) if present else 1


def period_word(word):
    return next(d for d in divisors(len(word))
                if all(word[j] == word[(j+d) % len(word)] for j in range(len(word))))


def cell_audit():
    cases = [
        ('unit square', [(-1,0,0,0),(0,-1,0,0),(1,0,1,0),(0,1,1,0)]),
        ('full cross descriptor', [(-1,0,0,1),(0,-1,0,1),(2,0,1,0),(0,1,1,0)]),
        ('2x+3y triangle', [(-1,0,0,0),(0,-1,0,0),(2,3,1,0)]),
        ('nonprimitive rectangle', [(-1,0,0,0),(0,-1,0,0),(2,0,1,0),(0,3,1,0)]),
        ('rational translating triangle', [(-6,0,-3,-2),(0,-10,-2,-1),(5,5,7,3)]),
        ('thin slope one third', [(-1,0,0,0),(1,0,1,0),(1,-3,0,0),(-1,3,0,1)]),
        ('thin no lattice points', [(-1,0,0,0),(1,0,1,0),(0,-3,0,-1),(0,3,0,2)]),
        ('bounded translating box', [(-6,0,-3,-2),(6,0,3,5),(0,-3,0,-1),(0,3,0,4)]),
        ('integral oblique segment', [(1,1,1,0),(-1,-1,-1,0),(-3,0,0,-1),(2,0,2,1)]),
        ('nonintegral alternating segment', [(2,2,1,0),(-2,-2,-1,0),(-1,0,0,0),(1,0,1,0)]),
        ('vertical segment', [(2,0,1,1),(-2,0,-1,-1),(0,-1,0,0),(0,3,1,2)]),
        ('moving point', [(2,0,1,1),(-2,0,-1,-1),(0,3,1,0),(0,-3,-1,0)]),
        ('empty bounded cell', [(1,0,0,0),(-1,0,0,-1),(0,1,1,0),(0,-1,0,0)]),
        ('virtual negative side at parameter zero', [(-1,0,0,0),(0,-1,0,0),(1,0,1,-7),(0,1,1,-11)]),
    ]
    rng = random.Random(20261026)
    for j in range(24):
        cs = [(-1,0,0,rng.randrange(-3,4)),(0,-1,0,rng.randrange(-3,4)),
              (1,0,2,rng.randrange(-3,4)),(0,1,2,rng.randrange(-3,4))]
        for a,b,s in [(1,1,3),(2,1,4),(1,2,4)]:
            if rng.randrange(2):
                cs.append((a,b,s,rng.randrange(-3,4)))
        cases.append((f'random clipped box {j:02d}', cs))
    rows = []
    comparisons = 0
    for name, cs in cases:
        vs,N,E = stable_vertices(cs)
        tables = [[],[],[]]
        abc = []
        for h in range(E):
            A,B,C = cone_coefficients(vs,N,E,h)
            abc.append([str(A),str(B),str(C)])
            tables[0].append(A/E**2)
            tables[1].append(B/E - 2*A*h/E**2)
            tables[2].append(C - B*h/E + A*h*h/E**2)
            first = max(0, ceil(F(N-h,E)))
            for t in (first, first+1, first+3, first+7):
                k = E*t+h
                predicted = A*t*t+B*t+C
                actual = count_slices(cs,k)
                assert predicted == actual, (name,k,predicted,actual)
                comparisons += 1
        periods = [period_fourier(w) for w in tables]
        assert periods == [period_word(w) for w in tables]
        rows.append({'name':name,'dimension':min(2,len(vs)-1),
                     'certified_stable_from':N,'velocity_modulus':E,
                     'coefficient_periods':periods,
                     'A_B_C_in_t':abc})
        if name == '2x+3y triangle':
            assert tables[0] == [F(1,12)]*6
            assert tables[1] == [F(1,2)]*6
            assert tables[2] == [F(1),F(5,12),F(2,3),F(3,4),F(2,3),F(5,12)]
            assert periods == [1,1,6]
        if name == 'full cross descriptor':
            assert tables[0] == [F(1,2)]*2
            assert tables[1] == [F(3),F(5,2)]
            assert tables[2] == [F(4),F(3)]
            assert periods == [1,2,2]
    return {'families':len(cases),'exact_comparisons':comparisons,'cases':rows}


def complete_three_line_representations(n):
    sites = [(x,y) for x in range(n+1) for y in range(n+1-x)]
    sections = [dict() for _ in range(3)]
    sizes = [dict() for _ in range(3)]
    for j,(x,y) in enumerate(sites):
        for i in range(3):
            t = y-2*i*x
            sections[i][t] = sections[i].get(t,0) | (1 << j)
            sizes[i][t] = sizes[i].get(t,0)+1
    anchors = set()
    for i,j in combinations(range(3),2):
        for ti in sections[i]:
            for tj in sections[j]:
                delta = tj-ti-(j-i)
                denominator = 2*(j-i)
                if delta % denominator == 0:
                    X = delta//denominator
                    Y = i+2*i*X-ti
                    anchors.add((X,Y))
    multi = []
    for X,Y in anchors:
        masks = tuple(sections[i].get(i-Y+2*i*X,0) for i in range(3))
        assert sum(bool(b) for b in masks) >= 2
        assert not any(masks[i] & masks[j] for i,j in combinations(range(3),2))
        multi.append(masks)
    entries = multi + [(0,0,0)]
    for i in range(3):
        for b in sections[i].values():
            row = [0,0,0]
            row[i] = b
            entries.append(tuple(row))
    Fn = sum(1 for t in sections[1] if sizes[1][t] >= 3 and (t-n)%2
             and 2*t-n not in sections[0] and 2*t-n not in sections[2])
    return entries, Fn


def star_audit():
    rows = []
    totals = 0
    residuals = {c:[] for c in product((1,2),repeat=3)}
    classified_classes = 0
    for n in range(8,65):
        entries, Fn = complete_three_line_representations(n)
        counts = {}
        for colors in residuals:
            patterns = {}
            for masks in entries:
                one = two = 0
                for i in range(3):
                    if colors[i] == 1:
                        one |= masks[i]
                    else:
                        two |= masks[i]
                key = (one,two)
                if key in patterns:
                    count, first_masks = patterns[key]
                    long_ids = [i for i in range(3) if bin(first_masks[i]).count('1') >= 3]
                    new_long_ids = [i for i in range(3) if bin(masks[i]).count('1') >= 3]
                    assert long_ids == new_long_ids and len(long_ids) < 2
                    if long_ids:
                        i = long_ids[0]
                        assert first_masks[i] == masks[i]
                        off = (one | two) & ~masks[i]
                        assert 1 <= bin(off).count('1') <= 2
                        if bin(off).count('1') == 1:
                            # Site order starts (0,0),...,(0,n), hence V is bit n.
                            assert i == 1 and off == (1 << n)
                            assert colors[0] == colors[2] and count == 1
                    patterns[key] = (count+1, first_masks)
                else:
                    patterns[key] = (1, masks)
            linear_classes = 0
            for (one,two),(count,masks) in patterns.items():
                if count > 1:
                    classified_classes += 1
                    long_ids = [i for i in range(3) if bin(masks[i]).count('1') >= 3]
                    if len(long_ids) == 1:
                        off = (one | two) & ~masks[long_ids[0]]
                        if bin(off).count('1') == 1:
                            linear_classes += 1
            assert linear_classes == Fn * (colors[0] == colors[2])
            exact = len(patterns)
            beta = F(3)-F(colors[0] == colors[2],4)
            remainder = F(exact)-F(31,4)*n*n-beta*n
            residuals[colors].append(remainder)
            counts[''.join(map(str,colors))] = exact
            totals += 1
        # The infinite-scale coefficient is proved in the manuscript; this is
        # only an exact enumeration of its color-matching collision effect.
        gap = counts['121']-counts['122']
        assert abs(Fn-F(n,4)) <= 20
        assert abs(F(gap)+F(n,4)) <= 100
        if n in (8,9,16,32,48,64):
            rows.append({'n':n,'representations':len(entries),
                         'singleton_classes':Fn,'complexities':counts})
    return {'windows':57,'color_choices':8,'exact_pattern_counts':totals,
            'repeated_classes_checked_against_classification':classified_classes,
            'sample_rows':rows,
            'bounded_remainders_observed':{
                ''.join(map(str,c)):{'min':str(min(w)),'max':str(max(w))}
                for c,w in residuals.items()},
            'warning':'Finite checks do not establish boundedness, periods or the theorem.'}


def period_audit():
    rng = random.Random(2026102602)
    checks = 0
    for M in (1,2,3,4,6,8,12,18,24,30):
        for _ in range(12):
            d = rng.choice(divisors(M))
            base = [F(rng.randrange(-4,5),rng.randrange(1,5)) for _ in range(d)]
            word = (base*(M//d))
            assert period_fourier(word) == period_word(word)
            checks += 1
    # A supplied period can collapse completely; zero must have period one.
    assert period_fourier([F(0)]*18) == 1
    assert period_fourier([F(7)]*18) == 1
    # Two disjoint moving point cells have opposite parity and cancel periods.
    point_cases = [
        [(2,0,1,0),(-2,0,-1,0),(0,1,0,0),(0,-1,0,0)],
        [(2,0,1,1),(-2,0,-1,-1),(0,1,0,1),(0,-1,0,-1)],
    ]
    words = []
    for cs in point_cases:
        vs,N,E = stable_vertices(cs)
        assert E == 2
        w = [cone_coefficients(vs,N,E,h)[2] for h in range(E)]
        assert period_fourier(w) == 2
        words.append(w)
    total = [a+b for a,b in zip(*words)]
    assert total == [F(1),F(1)] and period_fourier(total) == 1
    return {'rational_word_checks':checks+2,'methods':'cyclotomic division and exact cyclic shifts',
            'disjoint_cell_cancellation':{'individual_periods':[2,2],'sum_period':1}}


def residue_reduction_audit():
    rng = random.Random(2026102603)
    for _ in range(1500):
        q = rng.randrange(1,8)
        r = rng.randrange(q)
        rho = (rng.randrange(q),rng.randrange(q))
        a = (rng.randrange(-5,6),rng.randrange(-5,6))
        b,c = rng.randrange(-6,7),rng.randrange(-7,8)
        k = rng.randrange(0,15)
        z = (rng.randrange(-10,11),rng.randrange(-10,11))
        u = add(scale(z,q),rho)
        original = dot(a,u) <= b*(q*k+r)+c
        reduced = dot(a,z) <= b*k+floor(F(b*r+c-dot(a,rho),q))
        assert original == reduced
    return {'exact_atom_checks':1500}


def main():
    result = {'arithmetic':'Python standard-library Fraction; deterministic',
              'cell_audit':cell_audit(),'residue_reduction_audit':residue_reduction_audit(),
              'period_audit':period_audit(),'star_audit':star_audit(),
              'scope':'No complete Presburger compiler is implemented; no coefficient is fitted.'}
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()

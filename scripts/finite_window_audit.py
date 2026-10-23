#!/usr/bin/env python3
"""Independent finite-window and certificate stress checks; stdlib only.
Uses complete quadrant threshold cells, including unbounded cells, and checks
all 511 nonempty subsets of a 3 by 3 window plus seeded additional windows.
Writes only finite_window_audit.json beside this script.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json
import random

HERE = Path(__file__).resolve().parent


def det(a, b):
    return a[0]*b[1]-a[1]*b[0]


def norm(z):
    return max(abs(z[0]), abs(z[1]))


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def ceiling(q):
    return -((-q.numerator)//q.denominator)


def hull(sites):
    points = sorted(set(sites))
    if len(points) <= 1:
        return points
    low, high = [], []
    for p in points:
        while len(low) >= 2 and det((low[-1][0]-low[-2][0], low[-1][1]-low[-2][1]),
                                    (p[0]-low[-1][0], p[1]-low[-1][1])) <= 0:
            low.pop()
        low.append(p)
    for p in reversed(points):
        while len(high) >= 2 and det((high[-1][0]-high[-2][0], high[-1][1]-high[-2][1]),
                                     (p[0]-high[-1][0], p[1]-high[-1][1])) <= 0:
            high.pop()
        high.append(p)
    return low[:-1]+high[:-1]


def avoids_square(poly, u, radius):
    # Exact separating-axis test, including point and segment hulls.
    axes = [(1, 0), (0, 1)]
    axes += [(b[1]-a[1], a[0]-b[0]) for a, b in zip(poly, poly[1:]+poly[:1])]
    for a, b in axes:
        if a == b == 0:
            continue
        extent = radius*(abs(a)+abs(b))
        values = [a*(x+u[0])+b*(y+u[1]) for x, y in poly]
        if min(values) > extent or max(values) < -extent:
            return True
    return False


def threshold_cells(cuts):
    cuts = sorted(set(cuts))
    result = [(cuts[0]-1, None, cuts[0]-1)]
    result += [(a, a, b-1) for a, b in zip(cuts, cuts[1:])]
    result += [(cuts[-1], cuts[-1], None)]
    return result


def quadrant(window):
    sites = sorted(window)
    occupied = set(sites)
    xs = threshold_cells([-x for x, y in sites])
    ys = threshold_cells([-y for x, y in sites])
    xmask = {a: sum(1 << bit for bit, (x, y) in enumerate(sites) if a+x >= 0)
             for a, lo, hi in xs}
    ymask = {b: sum(1 << bit for bit, (x, y) in enumerate(sites) if b+y >= 0)
             for b, lo, hi in ys}
    occurrence = defaultdict(list)
    for (a, lx, rx), (b, ly, ry) in product(xs, ys):
        occurrence[xmask[a] ^ ymask[b]].append((lx, rx, ly, ry))

    # Direct E_i, L_i and the resulting two-direction anchors.
    e1 = {s for s in occupied if add(s, (0, 1)) in occupied}
    e2 = {s for s in occupied if add(s, (1, 0)) in occupied}
    good = {(x, y) for x in {-1-s[0] for s in e2} for y in {-1-s[1] for s in e1}}
    f0 = {(0, 0), (-1, -1), (-1, 0), (0, -1)}
    central = set.intersection(*[{(f[0]-s[0], f[1]-s[1]) for s in occupied} for f in f0])
    poly = hull(sites)
    distant = {u for u in good if avoids_square(poly, u, 2)}
    counted = central | distant
    seen = set()
    for u in counted:
        mask = sum((int(u[0]+x >= 0) ^ int(u[1]+y >= 0)) << bit
                   for bit, (x, y) in enumerate(sites))
        cells = occurrence[mask]
        assert cells == [(u[0], u[0], u[1], u[1])], (window, u, cells)
        assert mask not in seen
        seen.add(mask)
    assert not central & distant
    assert len(counted) <= len(occurrence)
    return {'sites': len(sites), 'P': len(occurrence), 'A': len(good),
            'J': len(central), 'distant': len(distant), 'bound': len(counted)}


def basis(v):
    return next(e for e in product(range(-6, 7), repeat=2) if det(v, e) == 1)


def polynomial(shifts):
    terms = {(0, 0): 1}
    formal = {(0, 0)}
    for h in shifts:
        nxt = Counter()
        for z, c in terms.items():
            nxt[add(z, h)] += c
            nxt[z] -= c
        terms = {z: c for z, c in nxt.items() if c}
        formal |= {add(z, h) for z in list(formal)}
    return terms, formal


def certificate_trial(vectors, q, seed, empty=False):
    rng = random.Random(seed)
    p = 5
    m = len(vectors)
    bases = [basis(v) for v in vectors]
    kappas = [q*rng.choice([1, 2]) for _ in vectors]
    if m == 4:
        # Force H0+H1+H2=0: the fourth Q has no origin coefficient.
        kappas = [q]*m
    shifts = [(k*v[0], k*v[1]) for k, v in zip(kappas, vectors)]
    intervals = [(1, 0) if empty else (-1, 1) for v in vectors]
    tails = []
    transition = []
    for i in range(m):
        left = [[rng.randrange(p) for y in range(q)] for x in range(q)]
        right = [[(left[x][y]+1+rng.randrange(p-1)) % p for y in range(q)] for x in range(q)]
        tails.append((left, right))
        l, r = intervals[i]
        transition.append({(t, s): rng.randrange(p) for t in range(l, r+1) for s in range(q)})
    encoding = [5, -9, 0, 17, 2]

    def component(i, z):
        t = det(vectors[i], z)
        l, r = intervals[i]
        if t < l:
            return tails[i][0][z[0] % q][z[1] % q]
        if t > r:
            return tails[i][1][z[0] % q][z[1] % q]
        return transition[i][t, det(z, bases[i]) % q]

    @lru_cache(None)
    def eta(z):
        return encoding[sum(component(j, z) for j in range(m)) % p]

    polys = [polynomial([shifts[j] for j in range(m) if i != j]) for i in range(m)]

    def f(i, z):
        return sum(c*eta(add(z, a)) for a, c in polys[i][0].items())

    def g(i, sign, z):
        result = 0
        for a, c in polys[i][0].items():
            w = add(z, a)
            total = component(i, w)
            for j in range(m):
                if i == j:
                    continue
                side = int(sign*det(vectors[j], vectors[i]) > 0)
                total += tails[j][side][w[0] % q][w[1] % q]
            result += c*encoding[total % p]
        return result

    records = []
    samples = 0
    witnesses = []
    for i, (vi, ei) in enumerate(zip(vectors, bases)):
        formal = polys[i][1]
        betas = [det(vi, a) for a in formal]
        l, r = intervals[i]
        row_range = range(l-max(betas), r-min(betas)+1)
        end_rows = []
        for sign in [-1, 1]:
            nonzero = {t for t in row_range if any(g(i, sign, add((s*vi[0], s*vi[1]),
                                                                      (t*ei[0], t*ei[1])))
                                                       for s in range(kappas[i]))}
            assert nonzero
            end_rows.append((min(nonzero), max(nonzero)))
        assert end_rows[0] == end_rows[1]
        mu = end_rows[0]
        widths = []
        for k, vk in enumerate(vectors):
            margin = 1+max(abs(det(vk, a)) for a in formal)
            lk, rk = intervals[k]
            widths.append(max(abs(lk-margin), abs(rk+margin)))
        A = max(ceiling(F(norm(vectors[k])*widths[j]+norm(vectors[j])*widths[k],
                          abs(det(vectors[j], vectors[k])))) for j, k in combinations(range(m), 2))
        B = max(widths[i], (abs(vi[0])+abs(vi[1]))*A, abs(mu[0]), abs(mu[1]))
        T = 1+max(ceiling(F(B*abs(det(vectors[j], ei))+
                            max(abs(det(vectors[j], a)) for a in formal)+
                            max(abs(t) for t in intervals[j])+1,
                            abs(det(vectors[j], vi)))) for j in range(m) if i != j)
        R0i = max(norm(add((s*vi[0], s*vi[1]), (t*ei[0], t*ei[1])))
                  for s, t in product([-T, T], [-B, B]))
        # Rows near strip edges and profile extrema; far-end values test F9.
        ts = sorted({-B-1, -B, -B+1, B-1, B, B+1, mu[0], mu[1], 0})
        ss = sorted({-T-10, -T-1, -T, 0, T, T+1, T+10})
        for s, t in product(ss, ts):
            z = add((s*vi[0], s*vi[1]), (t*ei[0], t*ei[1]))
            value = f(i, z)
            samples += 1
            if abs(t) > B:
                assert value == 0
            elif abs(s) > T:
                assert value == g(i, 1 if s > 0 else -1, z)
        # Cartesian samples also attack bounded-core spill and non-unimodularity.
        for z in [(rng.randint(-3*B-3, 3*B+3), rng.randint(-3*B-3, 3*B+3))
                  for _ in range(300)]:
            value = f(i, z)
            samples += 1
            if abs(det(vi, z)) > B:
                assert value == 0
            elif abs(det(z, ei)) > T:
                assert value == g(i, 1 if det(z, ei) > 0 else -1, z)
        # Construct explicit nonzero f_i point, then finite strip-product marker.
        t = mu[0]
        s = next(s for s in range(kappas[i]) if g(i, 1, add((s*vi[0], s*vi[1]), (t*ei[0], t*ei[1]))))
        s += ((T-s)//kappas[i]+1)*kappas[i]
        witness = add((s*vi[0], s*vi[1]), (t*ei[0], t*ei[1]))
        assert f(i, witness) != 0
        witnesses.append(witness)
        records.append({'mu': mu, 'B': B, 'T': T, 'R0_i': R0i,
                        'origin_coefficient': polys[i][0].get((0, 0), 0),
                        'formal_size': len(formal), 'actual_support_size': len(polys[i][0])})
    cs = [r['B']+abs(r['mu'][0]) for r in records]
    K = max(ceiling(F(norm(vectors[j])*cs[i]+norm(vectors[i])*cs[j],
                      abs(det(vectors[i], vectors[j])))) for i, j in combinations(range(m), 2))
    d = witnesses[1][0]-witnesses[0][0], witnesses[1][1]-witnesses[0][1]
    marker_value = f(0, witnesses[0])*f(1, add(witnesses[0], d))
    assert marker_value != 0
    return {'directions': vectors, 'periods': kappas, 'tail_period': q,
            'empty_transition': empty, 'seed': seed, 'records': records,
            'R0': max(r['R0_i'] for r in records), 'K': K,
            'R': max(r['R0_i'] for r in records)+K+1,
            'pointwise_samples': samples, 'nonzero_marker_value': marker_value,
            'nonlinear_encoding': encoding}


def main():
    window_results = []
    base = [(x, y) for x in range(3) for y in range(3)]
    for bits in range(1, 1 << len(base)):
        window = [s for i, s in enumerate(base) if bits & (1 << i)]
        row = quadrant(window)
        window_results.append(row)
    rng = random.Random(20261023)
    random_results = []
    for trial in range(120):
        window = {(rng.randint(-12, 12), rng.randint(-12, 12)) for _ in range(rng.randint(1, 150))}
        random_results.append(quadrant(window))
    triangles = []
    for n in [4, 5, 10, 30, 31]:
        row = quadrant([(x, y) for x in range(n+1) for y in range(n+1-x)])
        row['n'] = n
        assert row['J'] == n*(n-1)//2
        assert row['distant'] == (n-4)*(n-3)//2
        assert row['bound'] == n*n-4*n+6
        triangles.append(row)
    vectors = [[(1, 0), (0, 1)], [(1, 0), (1, 2), (1, 4)],
               [(1, 0), (0, 1), (-1, -1), (1, -1)],
               [(2, 1), (-1, 2), (1, -1)]]
    certificates = [certificate_trial(v, q, 100+10*k+q, empty)
                    for k, v in enumerate(vectors) for q in [2, 3] for empty in [False, True]]
    cancelled, formal = polynomial([(1, 0), (0, 1), (-1, -1)])
    assert (0, 0) not in cancelled and (1, 0) in cancelled and (-1, 0) in cancelled
    assert (0, 0) in formal
    output = {'all_3x3_subsets': len(window_results),
              '3x3_nonzero_bounds': sum(r['bound'] > 0 for r in window_results),
              'random_windows': random_results, 'triangles': triangles,
              'certificate_trials': certificates,
              'certificate_samples': sum(r['pointwise_samples'] for r in certificates),
              'cancelled_origin_example': {'actual_support': sorted(cancelled),
                                           'origin_in_erosion_but_not_S': True,
                                           'origin_in_convex_hull': True},
              'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    dest = HERE/'finite_window_audit.json'
    dest.write_text(json.dumps(output, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: output[k] for k in ['all_3x3_subsets', '3x3_nonzero_bounds',
                                            'triangles', 'certificate_samples']}, ensure_ascii=False))


if __name__ == '__main__':
    main()

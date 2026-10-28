#!/usr/bin/env python3
"""Offline complete-plane checks for two-line stars; standard library only.

All writes, if requested, go beside this script.  This script imports no
workspace code, makes no network calls, and reads no environment files.
The manuscript proves the formulae independently of these finite checks.
"""

from collections import defaultdict
from fractions import Fraction
import json
from pathlib import Path


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def convex_sites(vertices, n):
    """Closed CCW convex polygon, using exact integer half-plane tests."""
    scaled = [(n * x, n * y) for x, y in vertices]
    return {
        (x, y)
        for x in range(min(a for a, b in scaled), max(a for a, b in scaled) + 1)
        for y in range(min(b for a, b in scaled), max(b for a, b in scaled) + 1)
        if all(
            det((b[0] - a[0], b[1] - a[1]), (x - a[0], y - a[1])) >= 0
            for a, b in zip(scaled, scaled[1:] + scaled[:1])
        )
    }


def rectangle(x0, x1, y0, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def union_sites(polys, n):
    return sorted(set().union(*(convex_sites(poly, n) for poly in polys)))


def pentagon(k, hole=False):
    """P_k, or its deeper-base version with one rectangular hole."""
    L = k + 1
    if hole:
        base = [rectangle(0, 1, -3, 0), rectangle(2, L, -3, 0),
                rectangle(1, 2, -3, -2), rectangle(1, 2, -1, 0)]
    else:
        base = [rectangle(0, L, -1, 0)]
    return base + [[(0, 0), (1, 0), (0, 1)],
                   [(1, 0), (L, 0), (L, 1)]]


def complete_count(sites, k, shear=0, literal=False):
    """theta=1[y=shear*x and x=0 mod k] XOR 1[x=0].

    Every anchor has independent integer parameters a=-X and
    h=shear*X-Y.  Active columns a are listed exactly.  An absent
    column has one representative for each a mod k.  Every actual
    section height h is listed; one absent height represents the rest.
    Thus these representatives exhaust the whole plane, not a box.
    """
    xs = sorted({x for x, y in sites})
    hs = sorted({y - shear * x for x, y in sites})
    inactive_h = hs[0] - 1
    inactive_as = [xs[0] - 1 - ((xs[0] - 1 - r) % k) for r in range(k)]
    aa = xs + inactive_as
    hh = hs + [inactive_h]
    patterns = set()
    if literal:
        for a in aa:
            for h in hh:
                X, Y = -a, -shear * a - h
                code = 0
                for j, (x, y) in enumerate(sites):
                    gx, gy = x + X, y + Y
                    c = int(gy == shear * gx and gx % k == 0) ^ int(gx == 0)
                    code |= c << j
                patterns.add(code)
        return len(patterns)
    columns = defaultdict(int)
    sections = defaultdict(int)
    for j, (x, y) in enumerate(sites):
        columns[x] |= 1 << j
        sections[(y - shear * x, x % k)] |= 1 << j
    for a in aa:
        col = columns[a]
        for h in hh:
            patterns.add(col ^ sections[(h, a % k)])
    return len(patterns)


def incidence_count(sites, k, shear=0):
    """Separate count using the proved long-column identity, no XOR sets."""
    columns = defaultdict(int)
    heights = [set() for _ in range(k)]
    for x, y in sites:
        columns[x] += 1
        heights[x % k].add(y - shear * x)
    assert min(columns.values()) >= 3
    phase_columns = [sum(x % k == r for x in columns) for r in range(k)]
    H = list(map(len, heights))
    return len(columns) + 1 + sum((c + 1) * h for c, h in zip(phase_columns, H))


def family_formula(k, n, hole=False):
    L = k + 1
    A = 1 + int(n % k != 0)
    return ((4 if hole else 2) * L * n * n
            + (5 if hole else 3) * L * n + 2
            + A * ((L * n) // k + 2))


def square_formula(n):
    return (4 * n * n + Fraction(25, 2) * n + 2 if n % 2 == 0
            else 4 * n * n + 13 * n + 2)


def rational(value):
    return str(Fraction(value))


def diagnostic_table(counts, modulus, leading):
    """Finite-difference diagnostic only, never proof of stabilization."""
    out = []
    for r in range(modulus):
        ns = [n for n in sorted(counts) if n % modulus == r]
        vals = [Fraction(counts[b] - counts[a], b - a) - leading * (a + b)
                for a, b in zip(ns, ns[1:])]
        out.append({"residue": r, "values": list(map(rational, vals)),
                    "stable_on_tested_range": len(set(vals)) <= 1})
    return out


def run():
    checks = {"description": "Finite checks only; see the independent mathematical proof.",
              "small_window_screen": [], "two_square_exact": [],
              "arbitrary_period_family": [], "hole_family": [],
              "literal_crosschecks": [], "assertion_count": 0}
    cases = {
        "separated_squares_axes": ([rectangle(0, 1, 0, 1), rectangle(0, 1, 2, 3)], 0, 2),
        "separated_squares_diagonal": ([rectangle(0, 1, 0, 1), rectangle(0, 1, 2, 3)], 1, 4),
        "square_frame": ([rectangle(0, 3, 0, 1), rectangle(0, 3, 2, 3),
                           rectangle(0, 1, 1, 2), rectangle(2, 3, 1, 2)], 0, 9),
        "L_region": ([rectangle(0, 3, 0, 1), rectangle(0, 1, 0, 3)], 0, 9),
        "concave_pentagon": (pentagon(2), 0, 6),
    }
    for name, (polys, shear, leading) in cases.items():
        counts = {}
        for n in range(2, 25):
            sites = union_sites(polys, n)
            c = complete_count(sites, 2, shear)
            assert c == incidence_count(sites, 2, shear)
            checks["assertion_count"] += 1
            counts[n] = c
            if name == "separated_squares_diagonal":
                assert c == square_formula(n)
                checks["two_square_exact"].append({"n": n, "count": c})
                checks["assertion_count"] += 1
            if n <= 4:
                lit = complete_count(sites, 2, shear, literal=True)
                assert c == lit
                checks["literal_crosschecks"].append({"case": name, "n": n, "count": c})
                checks["assertion_count"] += 1
        checks["small_window_screen"].append({"name": name, "word_period": 2,
            "leading": leading, "counts": counts,
            "beta_diagnostic": diagnostic_table(counts, 2, leading)})
    for k in range(2, 13):
        counts = {}
        for n in range(2, 3 * k + 5):
            sites = union_sites(pentagon(k), n)
            c = complete_count(sites, k)
            assert c == family_formula(k, n)
            assert c == incidence_count(sites, k)
            checks["assertion_count"] += 2
            counts[n] = c
            if k <= 5 and n <= 3:
                assert c == complete_count(sites, k, literal=True)
                checks["literal_crosschecks"].append({"case": "family", "k": k, "n": n, "count": c})
                checks["assertion_count"] += 1
        L = k + 1
        beta = [3 * L + Fraction(L, k) * (1 + int(s != 0)) for s in range(k)]
        gamma = [2 + (1 + int(s != 0)) * (2 - Fraction(s, k)) for s in range(k)]
        period = min(h for h in range(1, k + 1)
                     if all(beta[s] == beta[(s + h) % k] for s in range(k)))
        assert period == k
        checks["assertion_count"] += 1
        checks["arbitrary_period_family"].append({"k": k, "counts": counts,
            "beta_exact": list(map(rational, beta)), "gamma_exact": list(map(rational, gamma)),
            "minimal_beta_table_period": period})
    for k in range(2, 7):
        counts = {}
        for n in range(max(2, k), max(2, k) + 2 * k + 3):
            sites = union_sites(pentagon(k, hole=True), n)
            c = complete_count(sites, k)
            assert c == family_formula(k, n, hole=True)
            assert c == incidence_count(sites, k)
            counts[n] = c
            checks["assertion_count"] += 2
        checks["hole_family"].append({"k": k, "counts": counts})
    checks["status"] = "all checks passed"
    return checks


if __name__ == "__main__":
    report = run()
    payload = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    Path(__file__).with_suffix(".json").write_text(payload, encoding="utf-8")
    print(payload, end="")

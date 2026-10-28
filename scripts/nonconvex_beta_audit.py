#!/usr/bin/env python3
"""Independent, offline crosschecks of the nonconvex-beta counterexamples.

This program imports neither the previous search nor its geometry routines.
Literal pixel patterns use a proved complete list of whole-plane anchors.
No computation here replaces the manuscript's exact proof.
"""

from collections import defaultdict
from fractions import Fraction
import json
from pathlib import Path
import random


ROOT = Path(__file__).resolve().parent


def square_sites(n):
    return sorted((x, y) for x in range(n + 1)
                  for y in list(range(n + 1)) + list(range(2 * n, 3 * n + 1)))


def roof_sites(k, n, hole=False):
    """Integer inequalities for the outer polygon and optional open hole."""
    L = k + 1
    low = -3 * n if hole else -n
    out = []
    for x in range(L * n + 1):
        for y in range(low, n + 1):
            inside = y <= 0 or x + y <= n or x >= n + k * y
            in_hole = hole and n < x < 2 * n and -2 * n < y < -n
            if inside and not in_hole:
                out.append((x, y))
    return out


def statistics(S, k, shear):
    columns = defaultdict(int)
    sections = [set() for _ in range(k)]
    for x, y in S:
        columns[x] += 1
        sections[x % k].add(y - shear * x)
    assert min(columns.values()) >= 3
    c = [sum(x % k == r for x in columns) for r in range(k)]
    H = [len(h) for h in sections]
    count = len(columns) + 1 + sum((a + 1) * b for a, b in zip(c, H))
    return count, c, H


def literal_patterns(S, k, shear):
    """Evaluate theta on pixels using actual (X,Y), not the counting identity."""
    xs = sorted({x for x, y in S})
    heights = sorted({y - shear * x for x, y in S})
    # All observed columns, and one absent column in each residue class.
    aa = xs + [xs[0] - 1 - ((xs[0] - 1 - r) % k) for r in range(k)]
    hh = heights + [min(heights) - 1]
    patterns = set()
    for a in aa:
        for h in hh:
            X, Y = -a, -shear * a - h
            ones = []
            for j, (x, y) in enumerate(S):
                gx, gy = x + X, y + Y
                if (int(gx == 0) + int(gy == shear * gx and gx % k == 0)) % 2:
                    ones.append(j)
            patterns.add(tuple(ones))
    return len(patterns)


def roof_formula(k, n, hole=False):
    L, s = k + 1, n % k
    A = 1 if s == 0 else 2
    beta = (5 if hole else 3) * L + Fraction(L * A, k)
    gamma = 2 + A * (2 - Fraction(s, k))
    return (4 if hole else 2) * L * n * n + beta * n + gamma


def main():
    report = {"scope": "Independent finite crosschecks; see proof/integral-nonconvex-windows.md.",
              "square_checks": [], "roof_checks": [], "hole_checks": [],
              "random_long_column_checks": [], "period_table_checks": []}
    for n in range(2, 31):
        S = square_sites(n)
        count, c, H = statistics(S, 2, 1)
        expected = 4 * n * n + (Fraction(25, 2) if n % 2 == 0 else 13) * n + 2
        assert count == expected
        assert H == ([4 * n + 1, 4 * n - 2] if n % 2 == 0 else [4 * n, 4 * n])
        row = {"n": n, "count": count, "column_counts": c, "height_counts": H}
        if n <= 8:
            row["literal_count"] = literal_patterns(S, 2, 1)
            assert row["literal_count"] == count
        report["square_checks"].append(row)
    for hole in (False, True):
        for k in range(2, 17):
            for n in range(2, 2 * k + 4):
                S = roof_sites(k, n, hole)
                count, c, H = statistics(S, k, 0)
                s = n % k
                assert H == [(4 if hole else 2) * n + int(r in {0, s}) for r in range(k)]
                assert count == roof_formula(k, n, hole)
                row = {"k": k, "n": n, "count": count}
                if k <= 6 and n <= 4:
                    row["literal_count"] = literal_patterns(S, k, 0)
                    assert row["literal_count"] == count
                report["hole_checks" if hole else "roof_checks"].append(row)
    rng = random.Random(20261028)
    for j in range(180):
        k = rng.randint(2, 8)
        shear = rng.randint(-3, 3)
        xs = rng.sample(range(-5, 6), rng.randint(1, 5))
        S = sorted((x, y) for x in xs for y in rng.sample(range(-5, 7), rng.randint(3, 7)))
        count, _, _ = statistics(S, k, shear)
        direct = literal_patterns(S, k, shear)
        assert direct == count
        report["random_long_column_checks"].append({"case": j, "k": k, "shear": shear,
                                                    "sites": len(S), "count": count})
    for k in range(2, 129):
        L = k + 1
        table = [3 * L + Fraction(L * (1 if s == 0 else 2), k) for s in range(k)]
        periods = [h for h in range(1, k + 1)
                   if all(table[s] == table[(s + h) % k] for s in range(k))]
        assert periods == [k]
        report["period_table_checks"].append({"k": k, "minimal_period": periods[0]})
    report["status"] = "all independent checks passed"
    payload = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    Path(__file__).with_suffix(".json").write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()

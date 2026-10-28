#!/usr/bin/env python3
"""Offline exact diagnostics for the nonconvex-window manuscript.

Only the Python standard library is used. No network, Git, subprocesses or
environment files. The general theorems are proved in proof/nonconvex-windows.md;
these finite checks are not their proof.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
import json
from pathlib import Path


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def hull(points):
    pts = sorted(set(points))
    if len(pts) < 3:
        return pts
    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and cross(sub(lower[-1], lower[-2]), sub(p, lower[-1])) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(sub(upper[-1], upper[-2]), sub(p, upper[-1])) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def area(poly):
    if len(poly) < 3:
        return F(0)
    return abs(sum(cross(a, b) for a, b in zip(poly, poly[1:] + poly[:1]))) / F(2)


def clip(poly, normal, bound):
    """Intersect a convex polygon with normal dot x <= bound, exactly."""
    if not poly:
        return []
    out = []
    val = lambda z: normal[0] * z[0] + normal[1] * z[1] - bound
    for a, b in zip(poly, poly[1:] + poly[:1]):
        va, vb = val(a), val(b)
        if va <= 0:
            out.append(a)
        if (va < 0 < vb) or (vb < 0 < va):
            t = va / (va - vb)
            out.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
    return hull(out)


def intersect_polygons(a, b):
    out = a
    for x, y in zip(b, b[1:] + b[:1]):
        edge = sub(y, x)
        normal = (edge[1], -edge[0])
        out = clip(out, normal, normal[0] * x[0] + normal[1] * x[1])
    return out


def union_area(polys):
    ans = F(0)
    for count in range(1, len(polys) + 1):
        for ids in combinations(range(len(polys)), count):
            p = polys[ids[0]]
            for j in ids[1:]:
                p = intersect_polygons(p, polys[j])
            ans += (-1) ** (count + 1) * area(p)
    return ans


def rect(x0, x1, y0, y1):
    return [(F(x0), F(y0)), (F(x1), F(y0)), (F(x1), F(y1)), (F(x0), F(y1))]


def merge_intervals(intervals):
    out = []
    for a, b in sorted(intervals):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def projections(polys, v):
    return merge_intervals([(min(cross(v, p) for p in poly), max(cross(v, p) for p in poly)) for poly in polys])


def solve_forms(v, w, a, b):
    # det(v,z)=a and det(w,z)=b.
    det = cross(v, w)
    return ((F(a) * w[0] - F(b) * v[0]) / det,
            (F(a) * w[1] - F(b) * v[1]) / det)


def band_intersection_area(directions, interval_lists):
    ans = F(0)
    for intervals in product(*interval_lists):
        poly = hull([solve_forms(directions[0], directions[1], a, b)
                     for a in intervals[0] for b in intervals[1]])
        for v, (lo, hi) in zip(directions[2:], intervals[2:]):
            poly = clip(poly, (-v[1], v[0]), hi)
            poly = clip(poly, (v[1], -v[0]), -lo)
        ans += area(poly)
    return ans


def shadow_area(polys, directions):
    intervals = [projections(polys, v) for v in directions]
    ans = F(0)
    for k in range(2, len(directions) + 1):
        for ids in combinations(range(len(directions)), k):
            ar = band_intersection_area([directions[i] for i in ids], [intervals[i] for i in ids])
            ans += (-1) ** k * (k - 1) * ar
    return ans


def in_poly(p, poly):
    return all(cross(sub(b, a), sub(p, a)) >= 0 for a, b in zip(poly, poly[1:] + poly[:1]))


def lattice_window(polys, n):
    scaled = [[(n * x, n * y) for x, y in p] for p in polys]
    xs = [x for p in scaled for x, y in p]
    ys = [y for p in scaled for x, y in p]
    lo_x, hi_x = min(xs).__floor__(), max(xs).__ceil__()
    lo_y, hi_y = min(ys).__floor__(), max(ys).__ceil__()
    return sorted((x, y) for x in range(lo_x, hi_x + 1) for y in range(lo_y, hi_y + 1)
                  if any(in_poly((x, y), p) for p in scaled))


def graph_parameters(window):
    adj = defaultdict(set)
    for x, y in window:
        a, b = (0, y), (1, x)
        adj[a].add(b)
        adj[b].add(a)
    unseen, types = set(adj), Counter()
    while unseen:
        root = next(iter(unseen))
        todo, component = [root], set()
        while todo:
            v = todo.pop()
            if v in component:
                continue
            component.add(v)
            todo.extend(adj[v] - component)
        unseen -= component
        types[(sum(v[0] == 0 for v in component), sum(v[0] == 1 for v in component))] += 1
    rows, cols = len({y for x, y in window}), len({x for x, y in window})
    isolated = types[(1, 1)]
    small = sum(types[t] for t in [(1, 2), (2, 1), (2, 2)])
    formula = (rows + 1) * (cols + 1) - 2 * isolated - isolated * (isolated - 1) // 2 - 2 * small
    return rows, cols, isolated, small, formula


def cross_pattern_count(window):
    """Complete plane enumeration, not a sampled anchor box.

    Each x or y anchor either hits one listed coordinate or hits none;
    a single inactive representative exhausts the latter possibility.
    The color at every site is evaluated directly.
    """
    xs, ys = {x for x, y in window}, {y for x, y in window}
    anchors_x = [-x for x in xs] + [1 - min(xs)]
    anchors_y = [-y for y in ys] + [1 - min(ys)]
    patterns = set()
    for u, v in product(anchors_x, anchors_y):
        code = 0
        for bit, (x, y) in enumerate(window):
            if ((u + x == 0) + (v + y == 0)) % 2:
                code |= 1 << bit
        patterns.add(code)
    return len(patterns)


def axis_star_count(window, kind):
    """Exhaustive for three explicit eventually periodic axis stars.

    Transition coordinates are at 0 or 1. Inside the transition envelope
    enumerate every integer; outside it one representative per parity and
    side suffices. This classification covers all anchors in Z^2.
    """
    xs, ys = [x for x, y in window], [y for x, y in window]
    ax = range(-max(xs) - 4, -min(xs) + 6)
    ay = range(-max(ys) - 4, -min(ys) + 6)
    if kind == "unequal_tails":
        color = lambda x, y: (2 * (x % 2 == 0) + (x % 2 == 0 and y >= 0)
                              + (x % 2 == 1 and x >= 1 and y % 2 == 0)) % 5
    elif kind == "disjoint_periodic_lines":
        color = lambda x, y: int(y == 0 and x % 2 == 0) + int(x == 1 and y % 2 == 0)
    elif kind == "two_steps":
        color = lambda x, y: int(y >= 0) + 2 * int(x >= 0)
    else:
        raise ValueError(kind)
    return len({bytes(color(u + x, v + y) for x, y in window) for u, v in product(ax, ay)})


def full_line_count(window, lines, prime):
    """Exact enumeration by zero, one, and at least two visible lines."""
    projections_by_line = [set(cross(v, s) for s in window) for v, a, weight in lines]
    patterns = {bytes(len(window))}
    for (v, a, weight), heights in zip(lines, projections_by_line):
        for t in heights:
            patterns.add(bytes(weight if cross(v, s) == t else 0 for s in window))
    anchors = set()
    for i, j in combinations(range(len(lines)), 2):
        v, a, _ = lines[i]
        w, b, _ = lines[j]
        for s, t in product(projections_by_line[i], projections_by_line[j]):
            z = solve_forms(v, w, a - s, b - t)
            if all(coordinate.denominator == 1 for coordinate in z):
                anchors.add(tuple(map(int, z)))
    for u in anchors:
        patterns.add(bytes(sum(weight for v, a, weight in lines
                               if cross(v, (u[0] + x, u[1] + y)) == a) % prime
                           for x, y in window))
    return len(patterns)


def check_all_small_windows():
    grid = list(product(range(4), repeat=2))
    checked, configurations = 0, Counter()
    for mask in range(1, 1 << len(grid)):
        s = [p for j, p in enumerate(grid) if (mask >> j) & 1]
        r, c, k, t, prediction = graph_parameters(s)
        actual = cross_pattern_count(s)
        assert prediction == actual, (s, prediction, actual)
        assert actual >= len(s) + 1
        configurations[(k, t)] += 1
        checked += 1
    return {"windows": checked, "all_exact_formula_checks_pass": True,
            "distinct_small_component_parameter_pairs": len(configurations)}


def main():
    shapes = {
        "L_shape": [rect(0, 3, 0, 1), rect(0, 1, 0, 3)],
        "square_frame": [rect(0, 1, 0, 3), rect(2, 3, 0, 3), rect(1, 2, 0, 1), rect(1, 2, 2, 3)],
        "two_diagonal_squares": [rect(0, 1, 0, 1), rect(3, 4, 3, 4)],
        "four_product_squares": [rect(a, a + 1, b, b + 1) for a, b in product([0, 3], repeat=2)],
    }
    axes = [(1, 0), (0, 1)]
    three = [(1, 0), (0, 1), (1, 1)]
    expected = {"L_shape": (5, 9), "square_frame": (8, 9),
                "two_diagonal_squares": (2, 4), "four_product_squares": (4, 4)}
    geometry = {}
    for name, polys in shapes.items():
        ar, lam, lam3 = union_area(polys), shadow_area(polys, axes), shadow_area(polys, three)
        assert (ar, lam) == expected[name]
        assert lam3 > ar
        geometry[name] = {"area": str(ar), "axis_shadow_area": str(lam),
                          "three_direction_shadow_area": str(lam3)}

    exhaustive = check_all_small_windows()
    special = []
    for n in [2, 3, 5, 8, 12, 20]:
        for name, polys in shapes.items():
            s = lattice_window(polys, n)
            actual = cross_pattern_count(s)
            assert graph_parameters(s)[-1] == actual
            predicted = (3 * n + 2) ** 2 if name in ("L_shape", "square_frame") else (2 * n + 3) ** 2
            assert actual == predicted
            special.append({"shape": name, "n": n, "sites": len(s), "patterns": actual})

    thin = []
    for length in [2, 3, 7]:
        for n in [2, 3, 5, 8, 12]:
            s = sorted(set(product(range(n + 1), repeat=2)) | {(j, j) for j in range(n + 1, length * n + 1)})
            actual = cross_pattern_count(s)
            expected_count = F(length * length + 2 * length - 1, 2) * n * n + F(5 * length + 3, 2) * n + 4
            assert actual == expected_count == graph_parameters(s)[-1]
            thin.append({"length": length, "n": n, "sites": len(s), "patterns": actual,
                         "true_quadratic_coefficient": str(F(length * length + 2 * length - 1, 2)),
                         "incorrect_projection_coefficient": length * length})

    star_checks = []
    for name in ["L_shape", "two_diagonal_squares", "four_product_squares"]:
        for n in [2, 4, 7, 10]:
            s = lattice_window(shapes[name], n)
            lam = shadow_area(shapes[name], axes)
            for kind in ["unequal_tails", "disjoint_periodic_lines", "two_steps"]:
                p = axis_star_count(s, kind)
                star_checks.append({"shape": name, "star": kind, "n": n, "patterns": p,
                                    "coefficient": str(lam), "residual_divided_by_n": str((p - lam * n * n) / n)})

    line_checks = []
    lines = [((1, 0), 0, 1), ((0, 1), 0, 1), ((1, 1), 0, 1)]
    for name in ["L_shape", "two_diagonal_squares"]:
        lam = shadow_area(shapes[name], three)
        for n in [2, 4, 7]:
            s = lattice_window(shapes[name], n)
            p = full_line_count(s, lines, 5)
            line_checks.append({"shape": name, "n": n, "patterns": p,
                                "coefficient": str(lam), "residual_divided_by_n": str((p - lam * n * n) / n)})

    result = {"purpose": "Finite exact diagnostics, not a proof of the asymptotic theorem",
              "geometry": geometry, "exhaustive_small_windows": exhaustive,
              "exact_cross_examples": special, "thin_appendage_counterexamples": thin,
              "axis_star_checks": star_checks, "three_line_checks": line_checks}
    payload = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    Path(__file__).with_suffix(".json").write_text(payload, encoding="utf-8")
    print(payload, end="")



if __name__ == "__main__":
    main()

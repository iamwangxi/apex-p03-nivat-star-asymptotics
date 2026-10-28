#!/usr/bin/env python3
"""Exact full-plane checks for Theorem 4, including sparse words and det > 1."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json

import nonconvex_window_audit as base


def bezout(a, b):
    if b == 0:
        return (1 if a > 0 else -1), 0
    s, t = bezout(b, a % b)
    return t, s - (a // b) * t


def prepare(v, origin, word):
    a, b = bezout(*v)
    assert a * v[0] + b * v[1] == 1
    assert any(word)
    assert all(any(word[j] != word[(j + shift) % len(word)] for j in range(len(word)))
               for shift in range(1, len(word)))
    return {"v": v, "origin": origin, "word": word, "dual": (a, b),
            "height": base.cross(v, origin), "complement": (-b, a)}


def color(line, x, y):
    v = line["v"]
    if base.cross(v, (x, y)) != line["height"]:
        return 0
    a, b = line["dual"]
    ox, oy = line["origin"]
    return line["word"][(a * (x - ox) + b * (y - oy)) % len(line["word"])]


def patterns(window, lines, prime):
    heights = [set(base.cross(line["v"], s) for s in window) for line in lines]
    states = {bytes(len(window))}
    for i, line in enumerate(lines):
        ox, oy = line["origin"]
        wx, wy = line["complement"]
        vx, vy = line["v"]
        for t in heights[i]:
            for phase in range(len(line["word"])):
                u, v = ox - t * wx + phase * vx, oy - t * wy + phase * vy
                states.add(bytes(color(line, u + x, v + y) for x, y in window))
    anchors = set()
    for t, s in product(*heights):
        u = base.solve_forms(lines[0]["v"], lines[1]["v"], lines[0]["height"] - t, lines[1]["height"] - s)
        if all(x.denominator == 1 for x in u):
            anchors.add(tuple(map(int, u)))
    for u, v in anchors:
        states.add(bytes(sum(color(line, u + x, v + y) for line in lines) % prime for x, y in window))
    return len(states), len(anchors), list(map(len, heights))


def product_polygons(directions, intervals):
    return [base.hull([base.solve_forms(*directions, x, y) for x in a for y in b])
            for a, b in product(*intervals)]


def main():
    configurations = [
        ([(1, 0), (0, 1)], [(0, 0), (0, 0)], [(1, 0), (1, 0)], 2, [12, 20]),
        ([(1, 0), (1, 2)], [(0, 0), (0, 1)], [(1, 0), (0, 1, 2)], 5, [24, 36]),
        ([(1, 0), (2, 3)], [(1, 0), (0, 1)], [(1, 0), (1, 0, 0)], 3, [36, 48]),
        ([(1, 1), (1, -1)], [(0, 0), (1, 0)], [(1,), (1, 0)], 2, [16, 24]),
    ]
    intervals = [[(F(0), F(1)), (F(3), F(4))], [(F(0), F(1)), (F(2), F(3))]]
    rows = []
    for case, (directions, origins, words, prime, sizes) in enumerate(configurations):
        lines = [prepare(v, o, w) for v, o, w in zip(directions, origins, words)]
        polys = product_polygons(directions, intervals)
        for n in sizes:
            window = base.lattice_window(polys, n)
            actual, anchors, counts = patterns(window, lines, prime)
            assert counts == [2 * (n + 1)] * 2
            prediction = anchors + sum(len(line["word"]) * count for line, count in zip(lines, counts)) + 1
            assert actual == prediction, (case, n, actual, prediction)
            assert actual > len(window)
            rows.append({"case": case, "directions": directions, "origins": origins,
                         "words": words, "prime": prime, "n": n, "sites": len(window),
                         "two_line_anchors": anchors, "coset_count_difference": anchors - len(window),
                         "patterns": actual, "formula": prediction})
    out = {"purpose": "Exact finite diagnostics for Theorem 4; the proof is in the manuscript",
           "cases": rows, "all_checks_pass": True}
    payload = json.dumps(out, ensure_ascii=False, indent=2) + "\n"
    Path(__file__).with_suffix(".json").write_text(payload, encoding="utf-8")
    print(payload, end="")



if __name__ == "__main__":
    main()

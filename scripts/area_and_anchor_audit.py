#!/usr/bin/env python3
"""Exact rational band areas, finite anchor counts and complete pattern sets.
Independent enumeration scans a proved sufficient rectangle for multiple-hit
anchors, then adds every pure section and zero. No author code is imported.
Writes only area_and_anchor_audit.json beside this script; standard library only.
Run: python3 -B scripts/area_and_anchor_audit.py
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json


HERE = Path(__file__).resolve().parent
VECTORS = ((1, 0), (1, 2), (1, 4))


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def section_length(x, required):
    """Exact calculation with rational arithmetic and integer masks."""
    intervals = [(F(0), F(1)), (2*x-2, 2*x+1), (4*x-4, 4*x+1)]
    endpoints = sorted({a for pair in intervals for a in pair})
    ans = F(0)
    for a, b in zip(endpoints, endpoints[1:]):
        middle = (a+b)/2
        present = {i for i, (lo, hi) in enumerate(intervals) if lo < middle < hi}
        if required(present):
            ans += b-a
    return ans


def integrated_area(required):
    # Exact finite check; see the mathematical completeness argument in proof/.
    lines = [(0, 0), (0, 1), (2, -2), (2, 1), (4, -4), (4, 1)]
    knots = sorted({F(d-b, a-c) for (a, b), (c, d) in combinations(lines, 2) if a != c})
    # Exact finite check; see the mathematical completeness argument in proof/.
    assert section_length(knots[0]-1, required) == 0
    assert section_length(knots[-1]+1, required) == 0
    total = F(0)
    pieces = []
    for a, b in zip(knots, knots[1:]):
        la, lb = section_length(a, required), section_length(b, required)
        mid = section_length((a+b)/2, required)
        assert mid == (la+lb)/2, (a, b, la, lb, mid)
        piece = (b-a)*(la+lb)/2
        total += piece
        pieces.append([str(a), str(b), str(la), str(lb), str(piece)])
    return total, pieces


def triangle(n):
    return [(x, y) for x in range(n+1) for y in range(n+1-x)]


def configuration(z, concurrent=False):
    x, y = z
    return sum(y-2*i*x == (0 if concurrent else i) for i in range(3)) % 2


def four_shifts(i):
    a, b = [VECTORS[j] for j in range(3) if j != i]
    return [(0, 0), a, b, add(a, b)]


def diff_at(i, z):
    origin, a, b, ab = four_shifts(i)
    return (configuration(add(z, ab))-configuration(add(z, a))
            -configuration(add(z, b))+configuration(z))


def anchor_membership(n):
    sites = triangle(n)
    allowed = set(sites)
    profiles, erosions, layers, good_heights = [], [], [], []
    for i in range(3):
        shifts = four_shifts(i)
        # Exact finite check; see the mathematical completeness argument in proof/.
        candidates = {i-det(VECTORS[i], q) for q in shifts}
        profile = {t: diff_at(i, (0, t)) for t in sorted(candidates)}
        profile = {t: value for t, value in profile.items() if value}
        lo, hi = min(profile), max(profile)
        assert all(diff_at(i, (0, t)) == profile.get(t, 0) for t in range(-50, 51))
        eroded = {s for s in sites if all(add(s, q) in allowed for q in shifts)}
        heights = {det(VECTORS[i], s) for s in eroded}
        good = {lo-t for t in heights} & {hi-t for t in heights}
        profiles.append(profile)
        erosions.append(len(eroded))
        layers.append(sorted(heights))
        good_heights.append(good)
    # Exact finite check; see the mathematical completeness argument in proof/.
    anchors = []
    pair_counts = Counter()
    triple_by_y = Counter()
    for x in range(-3*n-2, 2*n+3):
        for y in range(-6*n-2, 5*n+3):
            flags = [det(v, (x, y)) in good_heights[i] for i, v in enumerate(VECTORS)]
            for i, j in combinations(range(3), 2):
                pair_counts[f'{i}{j}'] += flags[i] and flags[j]
            if all(flags):
                triple_by_y[y] += 1
            if sum(flags) >= 2:
                anchors.append((x, y))
    return {
        'profiles': profiles, 'erosion_sizes': erosions, 'layers': layers,
        'good_heights': [sorted(c) for c in good_heights],
        'pair_counts': dict(pair_counts), 'triple_count': sum(triple_by_y.values()),
        'triple_by_y': dict(sorted(triple_by_y.items())), 'good_count': len(anchors),
    }, anchors


def global_patterns(n, concurrent=False, good=()):
    """Exact calculation with rational arithmetic and integer masks."""
    sites = triangle(n)
    rows = [defaultdict(int) for _ in range(3)]
    for position, (x, y) in enumerate(sites):
        for i in range(3):
            rows[i][y-2*i*x] |= 1 << position

    def mask(u):
        x, y = u
        masks = [rows[i].get((0 if concurrent else i)-y+2*i*x, 0) for i in range(3)]
        return masks[0] ^ masks[1] ^ masks[2], sum(bool(a) for a in masks)

    # Exact finite check; see the mathematical completeness argument in proof/.
    # Exact finite check; see the mathematical completeness argument in proof/.
    # Exact finite check; see the mathematical completeness argument in proof/.
    occurrences = Counter()
    multiple_anchors = 0
    for x in range(-3*n-2, 2*n+3):
        for y in range(-6*n-2, 5*n+3):
            value, hits = mask((x, y))
            if hits >= 2:
                occurrences[value] += 1
                multiple_anchors += 1

    singles = set()
    single_witnesses = 0
    tangent = 20*(n+1)
    for i in range(3):
        for t, row in rows[i].items():
            witness = (tangent, 2*i*tangent+(0 if concurrent else i)-t)
            value, hits = mask(witness)
            assert hits == 1 and value == row, (n, i, t, witness)
            singles.add(row)
            single_witnesses += 1
    zero, hits = mask((tangent, -10*tangent))
    assert zero == 0 and hits == 0
    patterns = set(occurrences) | singles | {0}
    good_errors = []
    for u in good:
        value, hits = mask(u)
        if hits < 2 or occurrences[value] != 1 or value in singles or value == 0:
            good_errors.append([u, hits, occurrences[value]])

    # Exact finite check; see the mathematical completeness argument in proof/.
    width = (len(sites)+7)//8
    fingerprint = sha256(b''.join(p.to_bytes(width, 'little') for p in sorted(patterns))).hexdigest()
    return {
        'n': n, 'window_size': len(sites), 'pattern_count': len(patterns),
        'multiple_anchor_count': multiple_anchors,
        'distinct_multiple_patterns': len(occurrences),
        'distinct_single_patterns': len(singles), 'single_witnesses_checked': single_witnesses,
        'good_anchor_failures': good_errors, 'independent_fingerprint': fingerprint,
        'residual_from_31n2_over4': str(F(len(patterns))-F(31*n*n, 4)),
    }


def main():
    output = {}
    area, pieces = integrated_area(lambda present: len(present) >= 2)
    pair_areas = [integrated_area(lambda present, i=i, j=j: i in present and j in present)[0]
                  for i, j in combinations(range(3), 2)]
    triple_area = integrated_area(lambda present: len(present) == 3)[0]
    assert area == F(31, 4)
    output['area'] = {'union_at_least_two': str(area), 'pairs': list(map(str, pair_areas)),
                      'triple': str(triple_area), 'vertical_integration_pieces': pieces}
    finite, good = anchor_membership(30)
    assert finite['good_count'] == 4059
    assert list(finite['pair_counts'].values()) == [577, 499, 3975]
    assert finite['triple_count'] == 496
    output['finite_n30'] = finite
    output['global_disjoint_lines'] = []
    for n in [10, 20, 30, 40, 60, 100]:
        row = global_patterns(n, good=good if n == 30 else ())
        if n == 30:
            assert row['pattern_count'] == 7037
        assert row['good_anchor_failures'] == []
        output['global_disjoint_lines'].append(row)
    output['global_concurrent_case_A'] = [global_patterns(n, concurrent=True) for n in [10, 20, 30, 40]]
    output['script_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    dest = HERE/'area_and_anchor_audit.json'
    dest.write_text(json.dumps(output, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'area': str(area), 'good_n30': finite['good_count'],
                      'global': [(r['n'], r['pattern_count']) for r in output['global_disjoint_lines']],
                      'output': dest.name}, ensure_ascii=False))


if __name__ == '__main__':
    main()

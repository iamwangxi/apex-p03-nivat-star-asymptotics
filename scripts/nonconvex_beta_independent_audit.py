#!/usr/bin/env python3
"""Offline adversarial checks; standard library only; no author-code imports.

Geometry is rasterized from ordered polygon boundaries by integer winding.
The literal counter evaluates theta(x+U,y+V) at each pixel.  A second exact
counter enumerates support events and deduplicates whole binary bitmaps; it
does not use the cardinality identity (8) to obtain the answer.
"""

from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import random

HERE = Path(__file__).resolve().parent


def cross(a, b, z):
    return (b[0] - a[0]) * (z[1] - a[1]) - (b[1] - a[1]) * (z[0] - a[0])


def in_closed_polygon(z, vertices):
    """Integer winding number, with exact inclusion of all boundary edges."""
    winding = 0
    for a, b in zip(vertices, vertices[1:] + vertices[:1]):
        d = cross(a, b, z)
        if d == 0 and min(a[0], b[0]) <= z[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= z[1] <= max(a[1], b[1]):
            return True
        if a[1] <= z[1] < b[1] and d > 0:
            winding += 1
        elif b[1] <= z[1] < a[1] and d < 0:
            winding -= 1
    return winding != 0


def polygon_window(shape, k, n):
    L = k + 1
    if shape == 'squares':
        polygons = [((0, 0), (1, 0), (1, 1), (0, 1)),
                    ((0, 2), (1, 2), (1, 3), (0, 3))]
    else:
        bottom = -1 if shape == 'pentagon' else -3
        polygons = [((0, bottom), (L, bottom), (L, 1), (1, 0), (0, 1))]
    polygons = [tuple((n * x, n * y) for x, y in p) for p in polygons]
    xmin = min(x for p in polygons for x, y in p)
    xmax = max(x for p in polygons for x, y in p)
    ymin = min(y for p in polygons for x, y in p)
    ymax = max(y for p in polygons for x, y in p)
    S = []
    for x in range(xmin, xmax + 1):
        for y in range(ymin, ymax + 1):
            if not any(in_closed_polygon((x, y), p) for p in polygons):
                continue
            if shape == 'hole' and n < x < 2 * n and -2 * n < y < -n:
                continue
            S.append((x, y))
    return tuple(S)


def piece_window(shape, k, n):
    """Separate geometry check using the advertised convex-piece union."""
    L = k + 1
    if shape == 'squares':
        return tuple((x, y) for x in range(n + 1) for y in range(3 * n + 1)
                     if y <= n or y >= 2 * n)
    bottom = -1 if shape == 'pentagon' else -3
    S = []
    for x in range(L * n + 1):
        for y in range(bottom * n, n + 1):
            left_triangle = x >= 0 and y >= 0 and x + y <= n
            right_triangle = y >= 0 and x <= L * n and x - n >= k * y
            if shape == 'pentagon':
                base = -n <= y <= 0
            else:
                base = (-3 * n <= y <= 0 and
                        (x <= n or x >= 2 * n or y <= -2 * n or y >= -n))
            if base or left_triangle or right_triangle:
                S.append((x, y))
    return tuple(S)


def theta(x, y, k, s):
    return int(x == 0) ^ int(y == s * x and x % k == 0)


def literal_patterns(S, k, s):
    """A provably complete finite anchor strip; literal pixel evaluation.

    U=-x activates the vertical line.  U in [-maxX-k,-minX+k]
    contains every such U and absent U of every residue.  Every possible
    sparse-line intersection has V=s*U-(y-s*x).  At each U, enumerate the
    entire interval from s*U-maxH-1 through s*U-minH+1, including an empty
    sparse-line state.  Any outside anchor has an identical state inside.
    This coverage proof does not assume the three-point-column condition.
    """
    xs = [x for x, y in S]
    hs = [y - s * x for x, y in S]
    patterns = set()
    anchors = 0
    for U in range(-max(xs) - k, -min(xs) + k + 1):
        for V in range(s * U - max(hs) - 1, s * U - min(hs) + 2):
            bits = 0
            for j, (x, y) in enumerate(S):
                if theta(x + U, y + V, k, s):
                    bits |= 1 << j
            patterns.add(bits)
            anchors += 1
    return patterns, anchors


def event_patterns(S, k, s):
    """Direct whole-bitmap count at all distinct actual support events."""
    X = {x for x, y in S}
    U_values = {-x for x in X}
    # Absent vertical line, with every possible sparse tangential phase.
    for residue in range(k):
        U = -max(X) - 1
        while U % k != residue:
            U -= 1
        assert -U not in X
        U_values.add(U)
    patterns = set()
    anchors = 0
    for U in sorted(U_values):
        vertical = 0
        sparse = defaultdict(int)
        for j, (x, y) in enumerate(S):
            X_abs = x + U
            if X_abs == 0:
                vertical |= 1 << j
            if X_abs % k == 0:
                V = s * X_abs - y
                sparse[V] |= 1 << j
        # The empty event is realized at every sufficiently large |V|.
        patterns.add(vertical)
        anchors += 1
        for bitmap in sparse.values():
            patterns.add(vertical ^ bitmap)
            anchors += 1
    return patterns, anchors


def height_table(S, k, s):
    columns = {x for x, y in S}
    c = [sum(x % k == r for x in columns) for r in range(k)]
    H = [len({y - s * x for x, y in S if x % k == r}) for r in range(k)]
    return c, H


def claimed_lemma(S, k, s):
    c, H = height_table(S, k, s)
    return sum(c) + 1 + sum((a + 1) * b for a, b in zip(c, H))


def claimed_family(shape, k, n):
    if shape == 'squares':
        return 4 * n * n + (13 * n if n % 2 else 25 * n // 2) + 2
    L = k + 1
    A = 1 if n % k == 0 else 2
    if shape == 'pentagon':
        return 2 * L * n * n + 3 * L * n + 2 + A * (L * n // k + 2)
    return 4 * L * n * n + 5 * L * n + 2 + A * (L * n // k + 2)


def expected_heights(shape, k, n):
    if shape == 'squares':
        return [4 * n + 1, 4 * n - 2] if n % 2 == 0 else [4 * n, 4 * n]
    baseline = (2 if shape == 'pentagon' else 4) * n
    return [baseline + int(r in {0, n % k}) for r in range(k)]


def digest_patterns(patterns, sites):
    width = (sites + 7) // 8
    h = hashlib.sha256()
    for bitmap in sorted(patterns):
        h.update(bitmap.to_bytes(width, 'little'))
    return h.hexdigest()


def fit_quadratic(points):
    # Exact Newton interpolation, independent of the displayed coefficients.
    (x0, y0), (x1, y1), (x2, y2) = points
    d01 = Fraction(y1 - y0, x1 - x0)
    d12 = Fraction(y2 - y1, x2 - x1)
    a = (d12 - d01) / (x2 - x0)
    b = d01 - a * (x0 + x1)
    c = y0 - a * x0 * x0 - b * x0
    return a, b, c


def min_table_period(table):
    return next(d for d in range(1, len(table) + 1)
                if all(table[r] == table[(r + d) % len(table)] for r in range(len(table))))


def main():
    output = {'description': __doc__,
              'seed': 20261028, 'family': [], 'geometry': [], 'lemma_exhaustive': [],
              'lemma_random': [], 'periods_from_direct_counts': [], 'negative_controls': []}
    counts = Counter()

    # Winding-number raster vs independent convex-piece union, including n=1.
    for shape in ('squares', 'pentagon', 'hole'):
        for k in ((2,) if shape == 'squares' else range(2, 21)):
            for n in range(1, 11):
                S = polygon_window(shape, k, n)
                assert S == piece_window(shape, k, n), (shape, k, n, 'geometry')
                counts['geometry_comparisons'] += 1
        output['geometry'].append({'shape': shape, 'passed': True})

    # Direct full-plane counts.  Extra large k cases attack n << k.
    for shape in ('squares', 'pentagon', 'hole'):
        cases = ([(2, n) for n in range(2, 41)] if shape == 'squares' else
                 [(k, n) for k in range(2, 17) for n in range(2, 11)] +
                 [(k, n) for k in (17, 19, 31, 32, 63, 64, 127, 128) for n in (2, 3, 4)])
        for k, n in cases:
            s = 1 if shape == 'squares' else 0
            S = polygon_window(shape, k, n)
            patterns, events = event_patterns(S, k, s)
            result = len(patterns)
            assert result == claimed_family(shape, k, n), (shape, k, n, result)
            assert result == claimed_lemma(S, k, s), (shape, k, n, 'lemma')
            assert height_table(S, k, s)[1] == expected_heights(shape, k, n), (shape, k, n, 'heights')
            raw = None
            # Each family has a substantial actual-pixel crosscheck subset.
            if (shape == 'squares' and n <= 14) or (shape != 'squares' and k <= 10 and n <= 7):
                literal, raw = literal_patterns(S, k, s)
                assert literal == patterns, (shape, k, n, 'literal bitmap mismatch')
                counts['literal_family_cases'] += 1
                counts['literal_anchor_evaluations'] += raw
                counts['literal_pixel_evaluations'] += raw * len(S)
            output['family'].append({'shape': shape, 'k': k, 'n': n, 'sites': len(S),
                                     'complexity': result, 'event_anchors': events,
                                     'literal_anchors': raw,
                                     'patterns_sha256': digest_patterns(patterns, len(S))})
            counts['direct_family_cases'] += 1

    # Exhaustive lemma stress: every qualifying window in 3 specified columns
    # and 4 possible heights, with nonconsecutive x, including singleton phases.
    choices = [q for t in (3, 4) for q in itertools.combinations((-2, -1, 1, 3), t)]
    for width in (1, 2, 3):
        for ys in itertools.product(choices, repeat=width):
            S = tuple((x, y) for x, values in zip((-3, 0, 4), ys) for y in values)
            for k, s in itertools.product((2, 3, 4, 7), (-3, 0, 2)):
                patterns, events = event_patterns(S, k, s)
                assert len(patterns) == claimed_lemma(S, k, s), (S, k, s)
                # Literal check every qualifying shape at k=3, s=-3.
                if (k, s) == (3, -3):
                    raw_patterns, raw = literal_patterns(S, k, s)
                    assert raw_patterns == patterns
                    counts['literal_lemma_cases'] += 1
                    counts['literal_anchor_evaluations'] += raw
                    counts['literal_pixel_evaluations'] += raw * len(S)
                counts['exhaustive_lemma_cases'] += 1
    output['lemma_exhaustive'].append({'columns': [-3, 0, 4], 'heights': [-2, -1, 1, 3],
                                     'column_site_counts': [3, 4], 'k': [2, 3, 4, 7],
                                     's': [-3, 0, 2], 'passed': True})

    rng = random.Random(20261028)
    for trial in range(400):
        xs = rng.sample(range(-10, 11), rng.randint(1, 7))
        S = tuple(sorted((x, y) for x in xs for y in rng.sample(range(-12, 13), rng.randint(3, 9))))
        k, s = rng.randint(2, 17), rng.randint(-7, 7)
        patterns, events = event_patterns(S, k, s)
        assert len(patterns) == claimed_lemma(S, k, s), (trial, S, k, s)
        literal_count = None
        if trial < 60:
            raw_patterns, raw = literal_patterns(S, k, s)
            assert raw_patterns == patterns
            literal_count = len(raw_patterns)
            counts['literal_lemma_cases'] += 1
            counts['literal_anchor_evaluations'] += raw
            counts['literal_pixel_evaluations'] += raw * len(S)
        output['lemma_random'].append({'trial': trial, 'k': k, 's': s, 'S': S,
                                       'complexity': len(patterns), 'literal_count': literal_count})
        counts['random_lemma_cases'] += 1

    # Recover coefficients from directly deduplicated patterns, not Eq. (8).
    # This is only a finite crosscheck; the proof of all n is reviewed separately.
    for shape in ('squares', 'pentagon', 'hole'):
        for k in ((2,) if shape == 'squares' else range(2, 11)):
            coefficients = []
            for r in range(k):
                points = []
                ns = [n for n in range(2, 4 * k + 2) if n % k == r][:4]
                for n in ns:
                    S = polygon_window(shape, k, n)
                    patterns, _ = event_patterns(S, k, int(shape == 'squares'))
                    points.append((n, len(patterns)))
                a, b, c = fit_quadratic(points[:3])
                for n, value in points:
                    assert a * n * n + b * n + c == value
                L = k + 1
                if shape == 'squares':
                    expected = (Fraction(4), Fraction(13) if r else Fraction(25, 2), Fraction(2))
                else:
                    A = 1 if r == 0 else 2
                    expected = (Fraction((2 if shape == 'pentagon' else 4) * L),
                                Fraction((3 if shape == 'pentagon' else 5) * L) + Fraction(L * A, k),
                                2 + A * (2 - Fraction(r, k)))
                assert (a, b, c) == expected, (shape, k, r, (a, b, c), expected)
                coefficients.append((a, b, c))
                counts['coefficient_progressions'] += 1
            beta = [q[1] for q in coefficients]
            gamma = [q[2] for q in coefficients]
            pb, pg = min_table_period(beta), min_table_period(gamma)
            assert pb == k and pg == (1 if shape == 'squares' else k)
            output['periods_from_direct_counts'].append({'shape': shape, 'k': k,
                 'beta': list(map(str, beta)), 'gamma': list(map(str, gamma)),
                 'p_beta': pb, 'p_gamma': pg})

    # Detect a known failure if the essential column threshold is weakened.
    S = ((0, 0), (0, 1), (1, 0), (1, 1))
    literal, raw = literal_patterns(S, 2, 0)
    event, _ = event_patterns(S, 2, 0)
    assert literal == event
    assert len(literal) == 7 and claimed_lemma(S, 2, 0) == 11
    output['negative_controls'].append({'description': 'two sites per column is insufficient',
              'S': S, 'k': 2, 's': 0, 'direct_count': len(literal),
              'incorrectly_extended_lemma': claimed_lemma(S, 2, 0)})

    # Check exact rational coefficient tables at periods beyond enumeration.
    for k in range(2, 1025):
        L = k + 1
        b = [Fraction(3 * L) + Fraction(L * (1 if r == 0 else 2), k) for r in range(k)]
        g = [2 + (1 if r == 0 else 2) * (2 - Fraction(r, k)) for r in range(k)]
        assert b.count(min(b)) == 1 and g.count(min(g)) == 1
        # Unique minimum proves every cyclic period fixes zero; no floating point.
        assert min_table_period(b) == k and min_table_period(g) == k
        counts['rational_period_tables'] += 1

    output['totals'] = dict(counts)
    output['all_passed'] = True
    payload = json.dumps(output, ensure_ascii=False, indent=2) + "\n"
    Path(__file__).with_suffix(".json").write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == '__main__':
    main()

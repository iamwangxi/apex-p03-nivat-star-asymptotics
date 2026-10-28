#!/usr/bin/env python3
"""Offline, standard-library checks for finite star descriptions and examples.

Run from any directory; writes only the adjacent deterministic text output.
Finite checks support the manuscript; they do not prove eventual periodicity.
Line-pattern enumeration is exhaustive over Z^2, by the representation proof
in the manuscript, and is not a finite translation-box sample.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path

HERE = Path(__file__).resolve().parent


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def norm_key(u):
    return abs(u[0]) + abs(u[1]), u[0], u[1]


def mixed_description_check():
    """Finite-table formula versus direct fields, with all difficult features.

    Directions (2,1),(1,-1) have determinant -3. The first component has
    three transition rows and tangential period 3. The second has an empty
    transition strip and period 2. Both have unequal nonconstant tails.
    A nonzero doubly periodic background is added. All arithmetic is in F_5.
    """
    prime, modulus = 5, 6
    v = ((2, 1), (1, -1))
    e = ((1, 1), (0, 1))
    periods = (3, 2)
    bounds = ((-1, 1), (1, 0))
    rows = {-1: (0, 4, 1), 0: (3, 0, 2), 1: (1, 2, 4)}
    left = (lambda x, y: (x + 2 * y) % 3,
            lambda x, y: (x - y) % 2)
    right = (lambda x, y: (2 * x + y + 1) % 3,
             lambda x, y: (2 * (x % 2) + y % 2 + 1) % 5)
    background = lambda x, y: (x % 2 + 2 * (y % 3)) % 5

    def direct_component(i, z):
        # Closed forms in Cartesian coordinates, independent of the table
        # predicates below.
        x, y = z
        if i == 0:
            height, phase = 2 * y - x, (x - y) % 3
            if height < -1:
                return left[0](x, y)
            if height > 1:
                return right[0](x, y)
            return rows[height][phase]
        return left[1](x, y) if x + y <= 0 else right[1](x, y)

    def residue_table(fn):
        return {a: {(x, y) for x, y in product(range(modulus), repeat=2)
                    if fn(x, y) == a} for a in range(prime)}

    left_tables = tuple(residue_table(fn) for fn in left)
    right_tables = tuple(residue_table(fn) for fn in right)
    bg_table = residue_table(background)
    transition = ({(h, b): rows[h][b] for h in rows for b in range(3)}, {})

    def formula_component(i, z, a):
        height, phase = det(v[i], z), det(z, e[i])
        ell, rr = bounds[i]
        residue = (z[0] % modulus, z[1] % modulus)
        return ((height < ell and residue in left_tables[i][a]) or
                (height > rr and residue in right_tables[i][a]) or
                any(height == h and (phase - b) % periods[i] == 0
                    and value == a
                    for (h, b), value in transition[i].items()))

    def direct(z):
        return (sum(direct_component(i, z) for i in range(2)) +
                background(*z)) % prime

    assert det(*v) == -3
    assert all(det(vv, ee) == 1 for vv, ee in zip(v, e))
    count = 0
    for z in product(range(-24, 25), repeat=2):
        for i in range(2):
            matches = [a for a in range(prime) if formula_component(i, z, a)]
            assert matches == [direct_component(i, z)]
            step = tuple(periods[i] * x for x in v[i])
            assert direct_component(i, add(z, step)) == direct_component(i, z)
        residue = z[0] % modulus, z[1] % modulus
        bg_matches = [a for a in range(prime) if residue in bg_table[a]]
        assert bg_matches == [background(*z)]
        # The finite disjunction defining the sum has exactly this color.
        values = [a for a in range(prime)
                  if any((a0 + a1 + b) % prime == a
                         and formula_component(0, z, a0)
                         and formula_component(1, z, a1)
                         and residue in bg_table[b]
                         for a0, a1, b in product(range(prime), repeat=3))]
        assert values == [direct(z)]
        count += 1
    # Actual Case A witness using common tail/background periods H_i=6v_i.
    h1, h2 = (12, 6), (6, -6)
    witness = None
    for z in product(range(-15, 16), repeat=2):
        vals = [direct(z), direct(add(z, h1)), direct(add(z, h2)),
                direct(add(add(z, h1), h2))]
        for a in range(prime):
            difference = ((vals[3] == a) - (vals[1] == a) -
                          (vals[2] == a) + (vals[0] == a))
            if difference:
                witness = z, a, difference
                break
        if witness:
            break
    assert witness is not None
    print('Mixed all-star finite-table checks: PASS;', count, 'sites; Case A witness', witness)
    # A finite simulation checks only the ordering code, not existence of
    # global representatives. Global existence is proved in the manuscript.
    window = ((0, 0), (1, 0), (0, 1), (1, 1))
    classes = defaultdict(list)
    for u in product(range(-8, 9), repeat=2):
        classes[tuple(direct(add(u, s)) for s in window)].append(u)
    for cls in classes.values():
        representative = min(cls, key=norm_key)
        assert sum(not any(norm_key(vv) < norm_key(u) for vv in cls)
                   for u in cls) == 1
        assert not any(norm_key(vv) < norm_key(representative) for vv in cls)
    print('Finite ordering sanity check: PASS;', len(classes), 'classes (finite sample only)')


def line_patterns(window, prime, lines):
    """All patterns, not box sampling. lines=(v,e,level,periodic_word).

    Every anchor hitting >=2 supporting lines is listed by solving a pair of
    projection equations. Every remaining pattern is a pure phase or zero.
    Pure phases are realized far along the primitive component direction.
    """
    assert all(det(v, e) == 1 for v, e, _, _ in lines)
    sections = []
    for v, _, _, _ in lines:
        groups = defaultdict(list)
        for j, s in enumerate(window):
            groups[det(v, s)].append(j)
        sections.append(groups)
    anchors = set()
    for i, j in combinations(range(len(lines)), 2):
        vi, _, ai, _ = lines[i]
        vj, _, aj, _ = lines[j]
        aa, bb, cc, dd = -vi[1], vi[0], -vj[1], vj[0]
        denominator = aa * dd - bb * cc
        assert denominator
        for ti, tj in product(sections[i], sections[j]):
            h, k = ai - ti, aj - tj
            xx, yy = h * dd - bb * k, aa * k - h * cc
            if xx % denominator == 0 and yy % denominator == 0:
                anchors.add((xx // denominator, yy // denominator))
    representations = []
    for u in sorted(anchors):
        colors = bytearray(len(window))
        for i, (v, e, level, word) in enumerate(lines):
            height = level - det(v, u)
            for j in sections[i].get(height, ()):
                colors[j] = (colors[j] + word[det(add(u, window[j]), e) % len(word)]) % prime
        representations.append(bytes(colors))
    for i, (_, e, _, word) in enumerate(lines):
        for indices in sections[i].values():
            for phase in range(len(word)):
                colors = bytearray(len(window))
                for j in indices:
                    colors[j] = word[(phase + det(window[j], e)) % len(word)]
                representations.append(bytes(colors))
    representations.append(bytes(len(window)))
    return len(set(representations)), len(representations), sections, len(anchors)


def exact_examples():
    triangle_lines = [((1, k), (0, 1), k // 2, (1,)) for k in (0, 2, 4)]
    print('Three-line checks: n, count, count-(31n^2+11n)/4, other collision excess')
    for n in (8, 9, 12, 13, 20, 21, 32, 33):
        window = tuple((x, y) for x in range(n + 1) for y in range(n + 1 - x))
        count, raw, sections, anchors = line_patterns(window, 2, triangle_lines)
        assert tuple(map(len, sections)) == (n + 1, 3 * n, 5 * n - 5)
        assert raw == anchors + 9 * n - 3
        families = sum(len(sections[1][t]) >= 3 and (t - n) % 2 == 1
                       and 2 * t - n not in sections[0]
                       and 2 * t - n not in sections[2] for t in sections[1])
        assert abs(Q(families) - Q(n, 4)) <= 20
        other = raw - count - families
        assert other >= 0
        # This loose bound is the previously proved bound, not a fit.
        remainder = Q(count) - Q(31 * n * n + 11 * n, 4)
        assert abs(remainder) <= 20 * (2 * 41**2)**2 + 100
        print(n, count, str(remainder), other)

    for n in (8, 9, 10, 11, 16, 17):
        aa, bb, r = n // 2 + 1, n + 1, n % 2
        window = tuple(product(range(aa), range(bb)))
        sparse = [((1, 0), (0, 1), 0, (1, 0)),
                  ((0, 1), (-1, 0), -1, (1, 0))]
        count = line_patterns(window, 2, sparse)[0]
        assert count == aa * bb + 2 * aa + 2 * bb + 1
        assert count == Q(n * n, 2) + (Q(9, 2) - Q(r, 2)) * n + 6 - Q(3 * r, 2)
        cross = [((1, 0), (0, 1), 0, (1,)),
                 ((0, 1), (-1, 0), 0, (1,))]
        for prime in (2, 3):  # cancellation at the crossing for prime 2
            count = line_patterns(window, prime, cross)[0]
            assert count == (aa + 1) * (bb + 1)
            assert count == Q(n * n, 2) + (3 - Q(r, 2)) * n + 4 - r
    print('Half-width rectangle checks: PASS; sparse beta=(9/2,4), gamma=(6,9/2); full beta=(3,5/2), gamma=(4,3)')

    colored = [((1, 0), (0, 1), 0, (1, 2)),
               ((1, 2), (0, 1), 1, (1, 2, 3))]
    for n in (3, 4, 5, 8, 9, 12):
        # P=conv((0,0),(1,0),(2,2),(1,2)), with determinant 2.
        window = tuple((x, y) for x in range(2 * n + 1) for y in range(2 * n + 1)
                       if -2 * n <= y - 2 * x <= 0)
        count, raw, _, anchors = line_patterns(window, 5, colored)
        assert anchors == 2 * n * (n + 1)
        assert count == raw == 2 * n * n + 12 * n + 6
    print('Full-support determinant-2, periods (2,3) check: PASS; exact 2n^2+12n+6')


if __name__ == '__main__':
    from contextlib import redirect_stdout
    from io import StringIO
    log = StringIO()
    with redirect_stdout(log):
        mixed_description_check()
        exact_examples()
        print('ALL CHECKS PASSED. No modulus fitting; no proof of EQP is inferred from finite data.')
    content = log.getvalue()
    (HERE / 'presburger_examples_audit.txt').write_text(content, encoding='utf8')
    print(content, end='')

#!/usr/bin/env python3
"""Exact end-profile stress checks for four-direction stars.
Includes nonlinear integer encodings, unequal nonconstant periodic tails,
periods 1/2/3, empty transition strips and cancellation of an origin exponent.
Finite checks support, but do not replace, the general proof.
Writes only end_profile_audit.json beside this script; standard library only.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import json
import random

HERE = Path(__file__).resolve().parent
DIRECTIONS = ((1, 0), (0, 1), (-1, -1), (1, -1))
COMPLEMENTS = ((0, 1), (-1, 0), (1, 0), (0, 1))
ENCODING = (7, -3, 0, 19, 5)


def det(v, z):
    return v[0]*z[1]-v[1]*z[0]


def linear(s, v, t, e):
    return s*v[0]+t*e[0], s*v[1]+t*e[1]


def formal_terms(indices, period):
    terms = []
    for choices in product((0, 1), repeat=len(indices)):
        x = sum(period*DIRECTIONS[i][0] for i, yes in zip(indices, choices) if yes)
        y = sum(period*DIRECTIONS[i][1] for i, yes in zip(indices, choices) if yes)
        terms.append(((x, y), (-1)**(len(indices)-sum(choices))))
    return terms


def run_one(seed, period, empty_strip):
    rng = random.Random(seed)
    lo, hi = (0, -1) if empty_strip else (-1, 1)
    tables = []
    for i in range(4):
        left = [[rng.randrange(5) for _ in range(period)] for _ in range(period)]
        right = [[rng.randrange(5) for _ in range(period)] for _ in range(period)]
        # Exact finite check; see the mathematical completeness argument in proof/.
        right[0][0] = (left[0][0]+1) % 5
        middle = {t: [rng.randrange(5) for _ in range(period)] for t in range(lo, hi+1)}
        tables.append((left, right, middle))

    def tail(i, z, right):
        s = det(z, COMPLEMENTS[i])
        t = det(DIRECTIONS[i], z)
        table = tables[i][int(right)]
        return table[s % period][t % period]

    def component(i, z):
        s = det(z, COMPLEMENTS[i])
        t = det(DIRECTIONS[i], z)
        if t < lo:
            return tail(i, z, False)
        if t > hi:
            return tail(i, z, True)
        return tables[i][2][t][s % period]

    def eta(z):
        return ENCODING[sum(component(i, z) for i in range(4)) % 5]

    def isolated(i, sign, z):
        value = component(i, z)
        for j in range(4):
            if i != j:
                value += tail(j, z, sign*det(DIRECTIONS[j], DIRECTIONS[i]) > 0)
        return ENCODING[value % 5]

    def apply(terms, field, z):
        x, y = z
        return sum(c*field((x+dx, y+dy)) for (dx, dy), c in terms)

    results = []
    all_terms = formal_terms(list(range(4)), period)
    case_a = None
    for x, y in product(range(-8*period, 8*period+1), repeat=2):
        value = apply(all_terms, eta, (x, y))
        if value:
            case_a = {'z': [x, y], 'D_eta': value}
            break
    assert case_a is not None

    zero_constant_indices = []
    for i, (v, e) in enumerate(zip(DIRECTIONS, COMPLEMENTS)):
        terms = formal_terms([j for j in range(4) if j != i], period)
        actual = defaultdict(int)
        for point, coefficient in terms:
            actual[point] += coefficient
        actual = {point: coefficient for point, coefficient in actual.items() if coefficient}
        if (0, 0) not in actual:
            zero_constant_indices.append(i)

        # Exact finite check; see the mathematical completeness argument in proof/.
        a = next(t for t in range(lo, hi+period+2)
                 if any(component(i, linear(s, v, t, e)) != tail(i, linear(s, v, t, e), False)
                        for s in range(period)))
        b = next(t for t in range(hi, lo-period-2, -1)
                 if any(component(i, linear(s, v, t, e)) != tail(i, linear(s, v, t, e), True)
                        for s in range(period)))
        beta_min = sum(min(0, period*det(v, DIRECTIONS[j])) for j in range(4) if j != i)
        beta_max = sum(max(0, period*det(v, DIRECTIONS[j])) for j in range(4) if j != i)
        expected = [a-beta_max, b-beta_min]

        # Exact finite check; see the mathematical completeness argument in proof/.
        projected = [det(v, q) for q, _ in terms]
        scan_lo, scan_hi = lo-max(projected)-2, hi-min(projected)+2
        ends = {}
        end_values = {}
        for sign in (-1, 1):
            profile = {}
            for t in range(scan_lo, scan_hi+1):
                row = [apply(terms, lambda z: isolated(i, sign, z), linear(s, v, t, e))
                       for s in range(period)]
                if any(row):
                    profile[t] = row
            extrema = [min(profile), max(profile)]
            assert extrema == expected, (seed, period, i, sign, expected, extrema)
            ends[str(sign)] = extrema
            end_values[sign] = profile

        # Exact finite check; see the mathematical completeness argument in proof/.
        max_t = max(abs(scan_lo), abs(scan_hi))
        thresholds = []
        for j in range(4):
            if j == i:
                continue
            error = max_t*abs(det(DIRECTIONS[j], e))
            error += max(abs(det(DIRECTIONS[j], q)) for q, _ in terms)
            thresholds.append(F(max(abs(lo), abs(hi))+error+1,
                                abs(det(DIRECTIONS[j], v))))
        threshold = max(thresholds)
        cutoff = threshold.numerator//threshold.denominator+1
        comparisons = 0
        for sign in (-1, 1):
            for t in range(scan_lo, scan_hi+1):
                for s in range(cutoff+1, cutoff+period+1):
                    z = linear(sign*s, v, t, e)
                    assert apply(terms, eta, z) == apply(terms, lambda zz: isolated(i, sign, zz), z)
                    comparisons += 1
        results.append({'direction': v, 'a_b': [a, b], 'expected_extrema': expected,
                        'end_extrema': ends, 'ends_have_different_values': end_values[-1] != end_values[1],
                        'end_comparisons': comparisons,
                        'actual_constant_coefficient': actual.get((0, 0), 0)})
    assert 3 in zero_constant_indices
    return {'seed': seed, 'tangential_period': period, 'empty_transition_strip': empty_strip,
            'case_A_certificate': case_a, 'Q_with_zero_constant': zero_constant_indices,
            'directions': results}


def main():
    rows = [run_one(20261023+trial, period, empty)
            for trial, (period, empty) in enumerate(product((1, 2, 3), (False, True)))]
    path = HERE/'end_profile_audit.json'
    path.write_text(json.dumps(rows, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'cases': len(rows), 'direction_end_pairs': sum(len(r['directions']) for r in rows),
                      'comparisons': sum(d['end_comparisons'] for r in rows for d in r['directions']),
                      'all_passed': True, 'output': path.name}, ensure_ascii=False))


if __name__ == '__main__':
    main()

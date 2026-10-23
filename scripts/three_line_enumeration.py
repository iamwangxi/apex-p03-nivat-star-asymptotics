#!/usr/bin/env python3
"""Complete exact enumeration for the three-line triangle, n=3,...,164.
Every multiple-hit anchor lies in the scanned rectangle (see the proof);
every pure section is separately realized, as is zero. Integer masks decide
pattern equality exactly. Collision classes are compared with the formula.
Writes only three_line_enumeration.json beside this script; stdlib only.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import gc
import json

HERE = Path(__file__).resolve().parent


def popcount(mask):
    return bin(mask).count('1')


def audit(n):
    sites = [(x, y) for x in range(n + 1) for y in range(n + 1 - x)]
    rows = [defaultdict(int) for _ in range(3)]
    for bit, (x, y) in enumerate(sites):
        for i in range(3):
            rows[i][y - 2*i*x] |= 1 << bit
    rows = [dict(d) for d in rows]
    xmin, xmax = -3*n-2, 2*n+2
    ymin, ymax = -6*n-2, 5*n+2
    height = ymax-ymin+1
    classes = {}
    pairs = Counter()
    triples = 0
    multiple = 0

    def insert(mask, code):
        old = classes.get(mask)
        if old is None:
            classes[mask] = code
        elif isinstance(old, int):
            classes[mask] = [old, code]
        else:
            old.append(code)

    r0, r1, r2 = rows
    for x in range(xmin, xmax+1):
        a, b = 1+2*x, 2+4*x
        for y in range(ymin, ymax+1):
            m0, m1, m2 = r0.get(-y, 0), r1.get(a-y, 0), r2.get(b-y, 0)
            f0, f1, f2 = bool(m0), bool(m1), bool(m2)
            if f0+f1+f2 < 2:
                continue
            assert not (m0 & m1 or m0 & m2 or m1 & m2)
            multiple += 1
            pairs['01'] += f0 and f1
            pairs['02'] += f0 and f2
            pairs['12'] += f1 and f2
            triples += f0 and f1 and f2
            insert(m0 | m1 | m2, (x-xmin)*height+y-ymin)

    pure = []
    for i, d in enumerate(rows):
        for t, mask in sorted(d.items()):
            pure.append((i, t))
            # A concrete one-hit witness, independent of counting congruences.
            x = 20*(n+1)
            y = i+2*i*x-t
            parts = [rows[j].get(j-y+2*j*x, 0) for j in range(3)]
            assert sum(bool(v) for v in parts) == 1 and parts[i] == mask
            insert(mask, -len(pure))
    insert(0, -len(pure)-1)

    def decode(code):
        if code < 0:
            if code == -len(pure)-1:
                return {'type': 'zero', 'parameters': [None]*3}
            i, t = pure[-code-1]
            return {'type': 'pure', 'direction': i, 't': t,
                    'parameters': [t if j == i else None for j in range(3)]}
        x, y = code//height+xmin, code%height+ymin
        return {'type': 'anchor', 'anchor': [x, y],
                'parameters': [j-y+2*j*x for j in range(3)]}

    family_actual = set()
    family_formula = {t for t, mask in r1.items()
                      if popcount(mask) >= 3 and (t-n) % 2 == 1
                      and 2*t-n not in r0 and 2*t-n not in r2}
    excess = 0
    hist = Counter()
    categories = Counter()
    exceptions = []
    collision_details = []
    for mask, codes in classes.items():
        if isinstance(codes, int):
            continue
        reps = [decode(c) for c in codes]
        hist[len(codes)] += 1
        excess += len(codes)-1
        long_sets = []
        for rep in reps:
            long_sets.append(tuple(i for i, t in enumerate(rep['parameters'])
                                   if popcount(rows[i].get(t, 0)) >= 3))
        assert len(set(long_sets)) == 1
        long = long_sets[0]
        if len(long) == 0:
            category = 'no_long'
        elif len(long) == 1:
            i = long[0]
            t = reps[0]['parameters'][i]
            external = mask ^ rows[i][t]
            count = popcount(external)
            category = f'long_{i}_external_{count}'
            if count == 1:
                assert i == 1 and len(codes) == 2 and t in family_formula
                assert external == r0[n] == r2[n]
                family_actual.add(t)
            else:
                assert count == 2
        else:
            raise AssertionError(('collision with >=2 long directions', n, reps))
        categories[category] += len(codes)-1
        if category != 'long_1_external_1':
            exceptions.append({'category': category, 'representations': reps})
        if n in (30, 159, 160, 161, 162, 163, 164):
            collision_details.append({'category': category, 'representations': reps,
                                      'pattern_size': popcount(mask)})
    assert family_actual == family_formula
    assert len(pure) == 9*n-4
    assert multiple == pairs['01']+pairs['02']+pairs['12']-2*triples
    R = multiple+len(pure)+1
    assert R-len(classes) == excess
    row = {'n': n, 'sites': len(sites), 'P': len(classes),
           'multi_anchors': multiple, 'pair_counts': dict(pairs), 'triple': triples,
           'pure_representations': len(pure), 'R': R, 'excess': excess,
           'class_size_histogram': dict(hist), 'category_excess': dict(categories),
           'family_count': len(family_actual), 'family_parameters': sorted(family_actual),
           'other_excess': excess-len(family_actual),
           'residual': str(Fraction(len(classes)) - Fraction(31*n*n+11*n, 4)),
           'multi_residual': str(Fraction(multiple)-Fraction(31*n*n,4)+6*n),
           'family_residual': str(Fraction(len(family_actual))-Fraction(n,4))}
    if collision_details:
        row['collision_classes'] = collision_details
    del classes
    gc.collect()
    return row


def main():
    results = []
    for n in range(3, 165):
        result = audit(n)
        results.append(result)
        if n % 10 == 0 or n >= 159:
            print(json.dumps({k: result[k] for k in
                              ['n', 'P', 'multi_anchors', 'excess', 'family_count',
                               'other_excess', 'residual']}, ensure_ascii=False), flush=True)
    residuals = [Fraction(r['residual']) for r in results]
    output = {'method': 'complete integer-anchor rectangle + all pure sections + zero; exact bitmask classes',
              'range': [3, 164], 'rows': results,
              'residual_min': str(min(residuals)), 'residual_max': str(max(residuals)),
              'max_other_excess': max(r['other_excess'] for r in results),
              'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    dest = HERE/'three_line_enumeration.json'
    dest.write_text(json.dumps(output, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: output[k] for k in ['range', 'residual_min', 'residual_max',
                                            'max_other_excess']}, ensure_ascii=False))


if __name__ == '__main__':
    main()

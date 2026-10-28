#!/usr/bin/env python3
"""Audit the direct boundary formula against 157 frozen finite count tables.

The source's finite-difference comparison is a cross-check, never a proof of
stabilization. Formula values are also compared with the reviewed source values.
A small fixed coverage set is recounted from exhaustive all-anchor lists.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import lcm
import json
import integral_boundary_coefficients as g

HERE=Path(__file__).resolve().parent


def main():
    source=json.loads((HERE/'integral_boundary_cases.json').read_text(encoding='utf-8'))['cases']
    out=[]
    covered=set()
    for row in source:
        z=row['star'];star=g.s.s.Star(z['dirs'],z['q'],z['fields'],z['name'])
        assert g.s.s.is_genuine_star(star)
        poly=g.s.s.WINDOWS[row['window']]
        data=g.general_formula(star,poly)
        # JSON normalization also checks every auxiliary collar/catalog value.
        normalized=json.loads(json.dumps(data,default=str))
        assert normalized==row['expected_formula'],('source formula mismatch',row['case'])
        M=star.q*lcm(*(abs(g.g.det(a,b)) for a,b in combinations(star.dirs,2)))
        diagnostic=g.s.s.diagnostic(row['ns'],row['counts'],M)
        assert diagnostic is not None,('finite diagnostic absent',row['case'])
        assert F(diagnostic['A'])==data['Lambda']
        assert all(F(b)==data['beta'] for b in diagnostic['beta'])
        recounted=[]
        # Cover every polygon, component type, period, direction set, and the
        # nonzero-collision input without recounting the exploratory screens.
        features={('polygon',row['window']),('mode',star.name),('period',star.q),
                  ('directions',tuple(star.dirs))}
        if data['Gamma'] or not features<=covered:
            for n in (12,20):
                actual=g.s.s.exact_count(star,poly,n)
                expected=row['counts'][row['ns'].index(n)]
                assert actual==expected,('exact count mismatch',row['case'],n)
                recounted.append({'n':n,'patterns':actual})
            covered.update(features)
        out.append({'case':row['case'],'window':row['window'],'mode':star.name,'q':star.q,
                    'formula':data,'finite_diagnostic':diagnostic,
                    'fresh_exact_counts':recounted,'matches':True})
    result={'purpose':'Finite exact cross-checks; the infinite-range proof is in proof/integral-windows.md, Appendix B.',
            'cases_checked':len(out),'freshly_recounted_cases':sum(bool(x['fresh_exact_counts']) for x in out),
            'mismatches':[],'cases':out}
    payload=json.dumps(result,indent=2,default=str)+'\n'
    Path(__file__).with_suffix('.json').write_text(payload,encoding='utf-8')
    print(payload,end='')


if __name__=='__main__':
    main()

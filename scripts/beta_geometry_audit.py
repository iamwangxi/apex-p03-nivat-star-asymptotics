#!/usr/bin/env python3
"""Exact diagnostics for the separately reviewed representative-polygon theorem.
Finite checks do not prove the theorem. No general star/QE compiler is supplied.
Standard library only; writes only the adjacent deterministic JSON.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, lcm, ceil, floor
from pathlib import Path
import json
import random

def det(a, b):
    return a[0]*b[1]-a[1]*b[0]

def sub(a,b):
    return (a[0]-b[0],a[1]-b[1])

def solve(a,b,c,d):
    dd=det(a,b)
    return (F(c*b[1]-a[1]*d,dd),F(a[0]*d-c*b[0],dd))

def at(v,k):
    return tuple(v[0][j]*k+v[1][j] for j in range(2))

def hull(points):
    ps=sorted(set(points))
    if len(ps)<=1:return ps
    lo=[]
    for p in ps:
        while len(lo)>1 and det(sub(lo[-1],lo[-2]),sub(p,lo[-1]))<=0:lo.pop()
        lo.append(p)
    hi=[]
    for p in reversed(ps):
        while len(hi)>1 and det(sub(hi[-1],hi[-2]),sub(p,hi[-1]))<=0:hi.pop()
        hi.append(p)
    return lo[:-1]+hi[:-1]

def primitive(v):
    den=lcm(v[0].denominator,v[1].denominator)
    w=[int(x*den) for x in v]
    gg=gcd(abs(w[0]),abs(w[1]))
    assert gg
    return tuple(x//gg for x in w)

def geometry(rows):
    """Rows (a,b,s,t): a*x+b*y <= s*k+t; primitive normals."""
    assert all(gcd(abs(a),abs(b))==1 for a,b,s,t in rows)
    candidates=[]
    K=2
    for i,j in combinations(range(len(rows)),2):
        a,b=rows[i][:2],rows[j][:2]
        if not det(a,b):continue
        v=(solve(a,b,rows[i][2],rows[j][2]),
           solve(a,b,rows[i][3],rows[j][3]))
        candidates.append(v)
        for aa,bb,s,t in rows:
            aa1=F(s)-aa*v[0][0]-bb*v[0][1]
            bb1=F(t)-aa*v[1][0]-bb*v[1][1]
            if aa1:K=max(K,ceil(abs(bb1/aa1))+2)
    candidates=list(dict.fromkeys(candidates))
    for v,w in combinations(candidates,2):
        for j in range(2):
            a=v[0][j]-w[0][j];b=v[1][j]-w[1][j]
            if a:K=max(K,ceil(abs(b/a))+2)
    good=[]
    for v in candidates:
        ok=True
        for a,b,s,t in rows:
            sl=F(s)-a*v[0][0]-b*v[0][1]
            ct=F(t)-a*v[1][0]-b*v[1][1]
            if sl<0 or (not sl and ct<0):ok=False
        if ok:good.append(v)
    pos={at(v,K):v for v in good}
    vertices=[pos[p] for p in hull(pos)]
    if len(vertices)<=1:
        return vertices,F(0),F(0),K
    if len(vertices)==2:
        av=sub(vertices[1][0],vertices[0][0])
        dv=sub(vertices[1][1],vertices[0][1])
        direction=primitive(av if any(av) else dv)
        j=0 if direction[0] else 1
        length=av[j]/direction[j]
        return vertices,F(0),abs(length),K
    a2=a1=boundary=F(0)
    for v,w in zip(vertices,vertices[1:]+vertices[:1]):
        a2+=det(v[0],w[0])/2
        a1+=(det(v[0],w[1])+det(v[1],w[0]))/2
        av=sub(w[0],v[0]);dv=sub(w[1],v[1])
        direction=primitive(sub(at(w,K),at(v,K)))
        j=0 if direction[0] else 1
        lam=av[j]/direction[j]
        assert lam>=0
        boundary+=lam
    return vertices,a2,a1+boundary/2,K

def lattice_count(rows,k):
    candidates=[]
    for i,j in combinations(range(len(rows)),2):
        a,b=rows[i][:2],rows[j][:2]
        if not det(a,b):continue
        x,y=solve(a,b,rows[i][2]*k+rows[i][3],rows[j][2]*k+rows[j][3])
        if all(aa*x+bb*y<=s*k+t for aa,bb,s,t in rows):
            candidates.append((x,y))
    if not candidates:return 0
    total=0
    for x in range(ceil(min(p[0] for p in candidates)),floor(max(p[0] for p in candidates))+1):
        low=[];high=[];ok=True
        for a,b,s,t in rows:
            rhs=s*k+t-a*x
            if b>0:high.append(floor(F(rhs,b)))
            elif b<0:low.append(ceil(F(rhs,b)))
            elif rhs<0:ok=False
        if ok:
            assert low and high
            total+=max(0,min(high)-max(low)+1)
    return total

def periods(rows):
    normals=[r[:2] for r in rows]
    delta=lcm(1,*[abs(det(a,b)) for a,b in combinations(normals,2) if det(a,b)])
    j=lcm(1,*[abs(a[1]) for a in normals if a[1]])
    return delta,j

def check_cell(name,rows,expected=None):
    vs,a,b,K=geometry(rows)
    if expected is not None:assert (a,b)==tuple(map(F,expected)),(name,a,b,expected)
    delta,j=periods(rows);period=delta*j
    # Four samples on every candidate residue verify the predicted coefficient
    # and constant remainder without estimating either coefficient from data.
    remainders={}
    for h in range(period):
        k=K+(h-K)%period
        vals=[]
        for z in range(4):
            kk=k+z*period
            vals.append(F(lattice_count(rows,kk))-a*kk*kk-b*kk)
        assert len(set(vals))==1,(name,h,vals)
        remainders[str(h)]=str(vals[0])
    return dict(name=name,quadratic=str(a),linear=str(b),stable_geometry_from=K,
                delta=delta,J=j,sufficient_period=period,
                tested_scales=4*period,remainder_by_residue=remainders)

def normalize(rows,q,r,rho,M):
    out=[]
    for a,b,s,t in rows:
        g=gcd(abs(a),abs(b))
        assert g and M%(q*g)==0
        out.append((a//g,b//g,s*M//(q*g),(s*r+t-a*rho[0]-b*rho[1])//(q*g)))
    return out

def descriptor_checks():
    cases={}
    # These are proved bijective pattern descriptors, not arbitrary anchor boxes.
    for name,sparse in [('full_cross',False),('sparse_cross',True)]:
        if not sparse:
            cells=[[(-1,0,0,1),(2,0,1,0),(0,-1,0,1),(0,1,1,0)]]
        else:
            cells=[
              [(-1,0,0,0),(2,0,1,0),(0,-1,0,0),(0,1,1,0)],
              [(-1,0,0,2),(1,0,0,-1),(0,-1,0,0),(0,1,1,0)],
              [(-1,0,0,0),(2,0,1,0),(0,-1,0,2),(0,1,0,-1)],
              [(-1,0,0,1),(1,0,0,-1),(0,-1,0,1),(0,1,0,-1)]]
        M=lcm(1,*[gcd(abs(a),abs(b)) for cell in cells for a,b,s,t in cell])
        coeffs=[]
        for r in range(M):
            norm=[normalize(c,1,r,(0,0),M) for c in cells]
            gg=[geometry(c) for c in norm]
            A=sum(g[1] for g in gg);B=sum(g[2] for g in gg)
            leading=A/M**2;beta=B/M-2*r*leading
            assert leading==F(1,2)
            assert beta==([F(9,2),F(4)][r] if sparse else [F(3),F(5,2)][r])
            for k in range(max(8,*[g[3] for g in gg]),60):
                n=M*k+r
                actual=sum(lattice_count(c,k) for c in norm)
                aa=n//2+1;bb=n+1
                exact=aa*bb+2*aa+2*bb+1 if sparse else (aa+1)*(bb+1)
                assert actual==exact
            coeffs.append(dict(residue=r,A=str(A),B=str(B),leading=str(leading),beta=str(beta)))
        cases[name]=dict(M_beta=M,coefficients=coeffs)
    # A nontrivial coordinate coset: x+y == n mod 3 in an n-square.
    rows=[(-1,0,0,0),(1,0,1,0),(0,-1,0,0),(0,1,1,0)]
    for r in range(3):
        cells=[normalize(rows,3,r,rho,3) for rho in product(range(3),repeat=2) if (sum(rho)-r)%3==0]
        gg=[geometry(c) for c in cells]
        A=sum(g[1] for g in gg);B=sum(g[2] for g in gg)
        assert A/9==F(1,3) and B/3-2*r*A/9==F(2,3)
        for k in range(4,20):
            n=3*k+r
            actual=sum(lattice_count(c,k) for c in cells)
            direct=sum((x+y-n)%3==0 for x in range(n+1) for y in range(n+1))
            assert actual==direct
    cases['coordinate_cosets']=dict(modulus=3,leading='1/3',beta='2/3',scales=48)
    return cases

def eval_atom(atom,x,y):
    kind,a,b,c,d=atom
    v=a*x+b*y+c
    return v<=0 if kind=='le' else ((v%d==0) if kind=='div' else (v%d!=0))

def eliminate_at(clause,y):
    L=lcm(1,*[abs(a[1]) for a in clause if a[1]])
    low=[];high=[];mod=[('div',1,0,L)]
    for kind,a,b,c,d in clause:
        tt=b*y+c
        if not a:
            if not eval_atom((kind,0,0,tt,d),0,0):return False
            continue
        cc=L//abs(a);sign=1 if a>0 else -1
        if kind=='le':
            (high if a>0 else low).append(-cc*tt if a>0 else cc*tt)
        else:mod.append((kind,sign,cc*tt,d*cc))
    Q=lcm(*[m[3] for m in mod])
    if low:xs=range(max(low),max(low)+Q)
    elif high:xs=range(min(high)-Q+1,min(high)+1)
    else:xs=range(Q)
    for x in xs:
        if any(x<z for z in low) or any(x>z for z in high):continue
        if all(((a*x+t)%d==0)==(kind=='div') for kind,a,t,d in mod):return True
    return False

def elimination_checks():
    rng=random.Random(20261025);tests=0
    for trial in range(240):
        clause=[('le',1,0,-7,1),('le',-1,0,-6,1)]
        for _ in range(rng.randrange(1,6)):
            kind=rng.choice(['le','div','ndiv'])
            clause.append((kind,rng.randrange(-3,4),rng.randrange(-3,4),rng.randrange(-6,7),rng.randrange(2,5)))
        for y in range(-8,9):
            direct=any(all(eval_atom(a,x,y) for a in clause) for x in range(-6,8))
            assert direct==eliminate_at(clause,y),(clause,y)
            tests+=1
    one_sided=[
      [('le',1,-1,0,1),('div',2,0,1,4)],
      [('le',-1,1,0,1),('div',2,1,0,3)],
      [('div',2,1,0,4),('ndiv',3,0,0,3)],
      [('div',0,1,0,3)],
      []]
    for clause in one_sided:
        for y in range(-10,11):
            direct=any(all(eval_atom(a,x,y) for a in clause) for x in range(-100,101))
            assert direct==eliminate_at(clause,y)
            tests+=1
    return dict(evaluations=tests,random_seed=20261025,all_passed=True)

def main():
    cases=[
      ('nonunimodular_triangle',[(-1,0,0,0),(0,-1,0,0),(2,3,1,0)],('1/12','1/2')),
      ('thin_tilted_width_two',[(-2,-3,0,0),(2,3,0,2),(-1,0,0,0),(1,0,1,0)],('0','1')),
      ('thin_tilted_width_one',[(-2,-3,0,0),(2,3,0,1),(-1,0,0,0),(1,0,1,0)],('0','2/3')),
      ('moving_segment',[(2,3,1,0),(-2,-3,-1,0),(-1,0,0,0),(1,0,2,0)],('0','2/3')),
      ('moving_point',[(2,3,1,0),(-2,-3,-1,0),(-1,0,0,0),(1,0,0,0)],('0','0')),
      ('bounded_corner_cut',[(-1,0,0,0),(0,-1,0,0),(1,0,1,0),(0,1,1,0),(1,1,2,-2)],('1','2')),
      ('shifted_two_slopes',[(-1,0,0,0),(0,-1,0,0),(2,1,3,1),(-1,2,2,3)],None),
      ('bounded_moving_triangle',[(-1,0,-1,0),(0,-1,-2,0),(2,3,8,1)],('0','0')),
    ]
    result=dict(status='finite diagnostics passed; supplementary proof established in one fresh-context GPT-6.1 Sol AI adversarial review on 2026-10-25',
                cell_checks=[check_cell(*c) for c in cases],
                descriptors=descriptor_checks(),elimination=elimination_checks())
    # Nonzero rational recession need not imply a real cell is empty.
    # 2*x=1, y free has recession (0,1), but no integer point.
    assert all(2*x!=1 for x in range(-20,21))
    result['unbounded_empty_lattice_cell']='2*x=1 with y free; no integer point by parity'
    target=Path(__file__).with_suffix('.json')
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(json.dumps(dict(cell_cases=len(cases),cell_scales=sum(c['tested_scales'] for c in result['cell_checks']),
                         elimination_evaluations=result['elimination']['evaluations'],
                         descriptor_betas={k:v for k,v in result['descriptors'].items()},
                         status='PASS'),ensure_ascii=False))
if __name__=='__main__':main()

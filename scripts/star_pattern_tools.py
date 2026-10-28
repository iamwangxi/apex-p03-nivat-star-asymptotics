#!/usr/bin/env python3
"""Offline exact finite pattern counts for small finitely specified stars.

Coverage: all anchors with >=2 active strips are enumerated by height pairs;
exactly one active strip by height, tangent residue and inactive-tail choices;
zero active strips by all open direction sectors, after a certified common-
interior threshold. No finite anchor box is declared exhaustive. Finite differences, when used for comparison, are diagnostic only; they are
not certified eventual coefficient tables. Cone-based exact tables, when supplied, use the reviewed period code.
"""
from fractions import Fraction as F
from itertools import product, combinations
from math import gcd, lcm

import minimal_period_audit as period


def det(v, z):
    return v[0]*z[1]-v[1]*z[0]


def ceil(x):
    return -((-x.numerator)//x.denominator)


def basis(v):
    for x in range(-20,21):
        for y in range(-20,21):
            if det(v,(x,y)) == 1:
                return (x,y)
    raise ValueError(v)


def hull(vertices):
    return period.hull(vertices)


WINDOWS = {
    'triangle': [(0,0),(1,0),(0,1)],
    'square': [(0,0),(1,0),(1,1),(0,1)],
    'diamond': [(0,-1),(1,0),(0,1),(-1,0)],
    'hexagon': [(0,0),(1,0),(2,1),(2,2),(1,2),(0,1)],
    'skew_triangle': [(0,0),(2,0),(0,1)],
}


def sites_of(vertices,n):
    vs = hull(vertices)
    return [(x,y) for x in range(n*min(p[0] for p in vs),n*max(p[0] for p in vs)+1)
            for y in range(n*min(p[1] for p in vs),n*max(p[1] for p in vs)+1)
            if all(det((b[0]-a[0],b[1]-a[1]),(x-n*a[0],y-n*a[1])) >= 0
                   for a,b in zip(vs,vs[1:]+vs[:1]))]


class Star:
    def __init__(self,dirs,q,fields,name=''):
        self.dirs = list(map(tuple,dirs))
        self.q = q
        self.fields = fields
        self.name = name
        self.bases = [basis(v) for v in self.dirs]

    def color_component(self,i,z):
        f = self.fields[i]
        a = det(self.dirs[i],z)
        table = f['left'] if a<f['ell'] else f['right'] if a>f['r'] else f['rows'][a-f['ell']]
        return table[(z[0]%self.q)*self.q+z[1]%self.q]

    def color(self,z):
        out = 0
        for i in range(len(self.dirs)):
            out ^= self.color_component(i,z)
        return out

    def to_json(self):
        return {'name':self.name,'dirs':self.dirs,'q':self.q,'fields':self.fields}


def active_intervals(star,vertices,n):
    return [(f['ell']-n*max(det(v,p) for p in vertices),
             f['r']-n*min(det(v,p) for p in vertices))
            for v,f in zip(star.dirs,star.fields)]


def certified_threshold(star,vertices):
    c = (F(sum(p[0] for p in vertices),len(vertices)),
         F(sum(p[1] for p in vertices),len(vertices)))
    N=1
    for v,f in zip(star.dirs,star.fields):
        a,b = min(det(v,p) for p in vertices),max(det(v,p) for p in vertices)
        N=max(N,int(max(F(f['ell'],b-det(v,c)),F(-f['r'],det(v,c)-a)))+1)
    return N


def sector_rays(dirs):
    # A primitive fan, with explicit interior rays between consecutive lines.
    rays=set()
    for R in range(1,50):
        for x in range(-R,R+1):
            for y in (-R,R):
                if x or y:
                    z=(x,y)
                    s=tuple(1 if det(v,z)>0 else -1 for v in dirs)
                    if all(det(v,z) for v in dirs):
                        rays.add((s,z))
        for y in range(-R+1,R):
            for x in (-R,R):
                z=(x,y)
                if all(det(v,z) for v in dirs):
                    rays.add((tuple(1 if det(v,z)>0 else -1 for v in dirs),z))
        bysign={}
        for s,z in sorted(rays):
            bysign[s]=z
        if len(bysign)==2*len(dirs):
            return list(bysign.values())
    raise ValueError('ray search failed')


def representatives(star,vertices,n):
    assert n>=certified_threshold(star,vertices)
    dirs,q=star.dirs,star.q
    ints=active_intervals(star,vertices,n)
    anchors=set()
    # At least two active components.
    for i,j in combinations(range(len(dirs)),2):
        vi,vj=dirs[i],dirs[j]
        D=det(vi,vj)
        for a in range(ints[i][0],ints[i][1]+1):
            for b in range(ints[j][0],ints[j][1]+1):
                X=vj[0]*a-vi[0]*b
                Y=vj[1]*a-vi[1]*b
                if X%D==0 and Y%D==0:
                    anchors.add((X//D,Y//D))
    # Exactly one active component.
    for i,v in enumerate(dirs):
        e=star.bases[i]
        others=[j for j in range(len(dirs)) if j!=i]
        for a in range(ints[i][0],ints[i][1]+1):
            u0=(a*e[0],a*e[1])
            for signs in product((-1,1),repeat=len(others)):
                lo=hi=None
                for j,s in zip(others,signs):
                    B=det(dirs[j],v)
                    C=det(dirs[j],u0)
                    # lower tail: B*t+C <= low-1; upper: >= up+1.
                    if s==1:
                        B,C=-B,-C
                        target=-ints[j][1]-1
                    else:
                        target=ints[j][0]-1
                    if B>0:
                        H=(target-C)//B
                        hi=H if hi is None else min(hi,H)
                    else:
                        L=ceil(F(target-C,B))
                        lo=L if lo is None else max(lo,L)
                if lo is not None and hi is not None and lo>hi:
                    continue
                for r in range(q):
                    t=(lo+(r-lo)%q) if lo is not None else (hi-(hi-r)%q) if hi is not None else r
                    if (lo is None or t>=lo) and (hi is None or t<=hi):
                        anchors.add((u0[0]+t*v[0],u0[1]+t*v[1]))
    # All tails: every feasible choice has an open recession sector because
    # -n*c is strictly inside every active interval at the certified threshold.
    for d in sector_rays(dirs):
        for rho in product(range(q),repeat=2):
            T=1
            for v,(lo,hi) in zip(dirs,ints):
                B=q*det(v,d)
                C=det(v,rho)
                T=max(T,ceil(F(hi+1-C,B)) if B>0 else ceil(F(lo-1-C,B)))
            anchors.add((rho[0]+T*q*d[0],rho[1]+T*q*d[1]))
    return anchors,ints


def exact_count(star,vertices,n,direct=False):
    sites=sites_of(vertices,n)
    anchors,ints=representatives(star,vertices,n)
    if direct:
        return len({tuple(star.color((u[0]+s[0],u[1]+s[1])) for s in sites) for u in anchors})
    caches=[{} for _ in star.dirs]
    patterns=set()
    for u in anchors:
        bits=0
        rho=(u[0]%star.q,u[1]%star.q)
        for i,v in enumerate(star.dirs):
            a=det(v,u)
            a=max(ints[i][0]-1,min(a,ints[i][1]+1))
            key=(a,rho)
            if key not in caches[i]:
                # On a fixed congruence class, a clamped tail has its fixed
                # table; an active row needs the actual anchor height.
                f=star.fields[i]
                val=0
                for j,s in enumerate(sites):
                    h=a+det(v,s)
                    table=f['left'] if h<f['ell'] else f['right'] if h>f['r'] else f['rows'][h-f['ell']]
                    color=table[((rho[0]+s[0])%star.q)*star.q+(rho[1]+s[1])%star.q]
                    val|=color<<j
                caches[i][key]=val
            bits^=caches[i][key]
        patterns.add(bits)
    return len(patterns)


def line_field(q,word,v,level=0):
    e=basis(v)
    row=[word[det((x,y),e)%len(word)] for x,y in product(range(q),repeat=2)]
    z=[0]*(q*q)
    return {'ell':level,'r':level,'left':z,'right':z,'rows':[row]}


def random_star(rng,dirs,q,mode):
    fields=[]
    for v in dirs:
        z=[0]*(q*q)
        if mode=='sparse':
            word=[rng.randrange(2) for _ in range(q)]
            if not any(word): word[0]=1
            fields.append(line_field(q,word,v,rng.randrange(-1,2)))
        elif mode=='strip':
            fields.append({'ell':-1,'r':1,'left':z,'right':z,
                           'rows':[[rng.randrange(2) for _ in z] for _ in range(3)]})
        elif mode=='step':
            fields.append({'ell':0,'r':-1,
                           'left':[rng.randrange(2) for _ in z],
                           'right':[rng.randrange(2) for _ in z],'rows':[]})
        else:
            fields.append({'ell':0,'r':0,
                           'left':[rng.randrange(2) for _ in z],
                           'right':[rng.randrange(2) for _ in z],
                           'rows':[[rng.randrange(2) for _ in z]]})
    return Star(dirs,q,fields,mode)


def is_genuine_star(star):
    """Exact non-double-periodicity check on the finite input tables."""
    q=star.q
    for v,f in zip(star.dirs,star.fields):
        if f['left']!=f['right']:
            continue
        differs=False
        for a in range(f['ell'],f['r']+1):
            for x,y in product(range(q),repeat=2):
                if (det(v,(x,y))-a)%q==0:
                    k=x*q+y
                    if f['rows'][a-f['ell']][k]!=f['left'][k]:
                        differs=True
        if not differs:
            return False
    return True


def diagnostic(ns,counts,M):
    vals=dict(zip(ns,counts))
    out={}
    for r in range(M):
        xs=[x for x in ns if x%M==r]
        if len(xs)<4: return None
        x=xs[-4]
        A=F(vals[x+2*M]-2*vals[x+M]+vals[x],2*M*M)
        B=F(vals[x+M]-vals[x],M)-A*(2*x+M)
        C=F(vals[x])-A*x*x-B*x
        if vals[x+3*M]!=A*(x+3*M)**2+B*(x+3*M)+C:
            return None
        out[r]=(A,B,C)
    if len({v[0] for v in out.values()})!=1:
        return None
    bs=[out[r][1] for r in range(M)]
    cs=[out[r][2] for r in range(M)]
    return {'M':M,'A':str(out[0][0]),'beta':list(map(str,bs)),
            'gamma':list(map(str,cs)),'p_beta_diagnostic':period.period_word(bs),
            'p_gamma_diagnostic':period.period_word(cs)}


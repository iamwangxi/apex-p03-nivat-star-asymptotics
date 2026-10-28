#!/usr/bin/env python3
"""Exact general-tail boundary formula of proof/integral-windows.md, Appendix B.

Computes end-profile labels, pure-mode feasibility, invisible vertex clusters
and cluster collisions using exact rationals. It is not the Presburger compiler.
The proof is in the note; finite count comparisons are only cross-checks.
"""
from fractions import Fraction as F
from itertools import combinations,product
from math import lcm

import integral_sparse_audit as s
g=s.g


def linear_interval(inequalities,lo,hi):
    ok,a,b=g.scalar_interval(list(inequalities)+[(-F(1),-lo),(F(1),hi)])
    return (a,b) if ok and a<b else None


def general_formula(star,vertices,prime=2,collar_scale=1,gradient_shift=None):
    vs=star.dirs;q=star.q;m=len(vs);es=[g.basis(v) for v in vs]
    symbol_bits=(prime-1).bit_length()
    poly=g.hull(vertices);fs=g.facets(poly)
    assert all(a.denominator==1 for z in poly for a in z),'integral polygon required'
    lows=[min(g.det(v,x) for x in poly) for v in vs]
    highs=[max(g.det(v,x) for x in poly) for v in vs]
    def Fi(i,z):return star.color_component(i,z)
    def tail(i,side,z):return star.fields[i]['left' if side==0 else 'right'][(int(z[0])%q)*q+int(z[1])%q]
    amin=[];amax=[]
    for i,(v,e,f) in enumerate(zip(vs,es,star.fields)):
        low=[];high=[]
        for a in range(f['ell']-q,f['r']+q+1):
            for r in range(q):
                z=(a*e[0]+r*v[0],a*e[1]+r*v[1])
                if Fi(i,z)!=tail(i,0,z):low.append(a)
                if Fi(i,z)!=tail(i,1,z):high.append(a)
        assert low and high,'periodic input component'
        amin.append(min(low));amax.append(max(high))
    trial_L=[0]+[x for a in range(1,m+2) for x in (a,-a)]
    h=gradient_shift or next((q,q*L) for L in trial_L if all(g.det(v,(1,L)) for v in vs))
    assert all(z%q==0 for z in h) and all(g.det(v,h) for v in vs)
    gradrows=[]
    for i,(v,e) in enumerate(zip(vs,es)):
        d=g.det(v,h);ar=[]
        for a in range(amin[i]-max(0,d),amax[i]+max(0,-d)+1):
            if any(Fi(i,(a*e[0]+r*v[0]+h[0],a*e[1]+r*v[1]+h[1]))!=Fi(i,(a*e[0]+r*v[0],a*e[1]+r*v[1])) for r in range(q)):
                ar.append(a)
        assert ar
        gradrows.append(ar)
    N=sum(map(len,gradrows))
    center=(sum(x[0] for x in poly)/len(poly),sum(x[1] for x in poly)/len(poly))
    erosion=g.ceil(max(F(max(0,g.dot(a,h)),b-g.dot(a,center)) for a,b in fs))
    R=[g.interior_row_constants(poly,v,2*N*q)['C']+
       max(amax[i]-min(gradrows[i]),max(gradrows[i])-amin[i],0)+g.ceil(erosion*(highs[i]-lows[i]))+2
       for i,v in enumerate(vs)]
    R=[collar_scale*max(R[i],g.interior_row_constants(poly,v,q)['C']+
       max(amax[i]-star.fields[i]['ell']+q,star.fields[i]['r']+q-amin[i],0)+2)
       for i,v in enumerate(vs)]
    # One coordinate dictionary per vertex makes differently oriented clusters
    # directly comparable, while large colored corner regions stay as bitsets.
    geometries={};pointsets={x:set() for x in poly};sideinfo={}
    for j,v in enumerate(vs):
        for side,hgt in [(0,lows[j]),(1,highs[j])]:
            qs=[x for x in poly if g.det(v,x)==hgt]
            if len(qs)!=1:continue
            x=qs[0];sideinfo[j,x]=side
            for d in range(R[j]+1):
                geometry=s.cone_section(poly,v,x,d if side==0 else -d)
                geometries[j,x,d]=geometry;pointsets[x].update(geometry)
    coords={x:sorted(zs) for x,zs in pointsets.items()}
    indices={x:{z:k for k,z in enumerate(zs)} for x,zs in coords.items()}
    cache={};catalog={x:set() for x in poly};eta=[];collars={}
    for j,(v,e) in enumerate(zip(vs,es)):
        ee=[]
        for side,hgt in [(0,lows[j]),(1,highs[j])]:
            qs=[x for x in poly if g.det(v,x)==hgt]
            if len(qs)!=1:ee.append(F(0));continue
            x=qs[0];weight=F(0);collars[j,x]=[]
            # minimum face: bulk right tail; maximum face: bulk left tail.
            bulk=1-side
            for d in range(R[j]+1):
                c=amax[j]-d if side==0 else amin[j]+d
                collars[j,x].append(c)
                zs=[z for dd in range(d+1) for z in geometries[j,x,dd]]
                for phase in range(q):
                    w=(c*e[0]+phase*v[0],c*e[1]+phase*v[1]);bits=0
                    for z in zs:
                        p=(w[0]+z[0],w[1]+z[1])
                        value=(Fi(j,p)-tail(j,bulk,p))%prime
                        if value:bits|=value<<(symbol_bits*indices[x][z])
                    cache[j,x,c,phase]=bits
                    if bits:catalog[x].add(bits)
                    else:weight+=F(1,q)
            ee.append(weight)
        eta.append(ee)
    labels={};pure_intervals={};pure=F(0)
    for i,(v,e) in enumerate(zip(vs,es)):
        others=[k for k in range(m) if k!=i]
        for signs in product((0,1),repeat=len(others)):
            sigma=tuple(signs);smap=dict(zip(others,signs))
            for r in range(q):
                word=0;jj=0
                for a in range(star.fields[i]['ell']-q,star.fields[i]['r']+q+1):
                    for phase in range(q):
                        z=(a*e[0]+(phase+r)*v[0],a*e[1]+(phase+r)*v[1])
                        val=Fi(i,z)
                        for k,side in smap.items():val=(val+tail(k,side,z))%prime
                        word|=val<<(symbol_bits*jj);jj+=1
                labels[i,sigma,r]=word
            lowers=[];uppers=[]
            for k,side in smap.items():
                A=g.det(vs[k],e);B=g.det(vs[k],v)
                # B*tau <= D-A*xi.
                D=-highs[k]
                if side==1:A,B,D=-A,-B,lows[k]
                branch=(F(D,B),F(-A,B))
                (uppers if B>0 else lowers).append(branch)
            conditions=[(lo1-hi1,hi0-lo0) for lo0,lo1 in lowers for hi0,hi1 in uppers]
            pure_intervals[i,sigma]=linear_interval(conditions,-highs[i],-lows[i])
        cuts={-highs[i],-lows[i]}
        for (j,sg),interval in pure_intervals.items():
            if j==i and interval:cuts.update(interval)
        for a,b in zip(sorted(cuts),sorted(cuts)[1:]):
            xi=(a+b)/2;present=set()
            for signs in product((0,1),repeat=len(others)):
                interval=pure_intervals[i,tuple(signs)]
                if interval and interval[0]<xi<interval[1]:
                    present.update(labels[i,tuple(signs),r] for r in range(q))
            pure+=(b-a)*len(present)
    raw=pure;area=F(0)
    for size in range(2,m+1):
        coef=(-1)**size*(size-1)
        for I in combinations(range(m),size):
            hp=[]
            for j in I:
                n=(-vs[j][1],vs[j][0]);hp.extend([(n,-lows[j]),((-n[0],-n[1]),highs[j])])
            area+=coef*g.polygon_area(g.halfplane_polygon(hp))
            for j in I:
                for side,level,delta in [(0,-lows[j],amax[j]),(1,-highs[j],-amin[j])]:
                    ok,a,b=g.scalar_interval([(g.dot(n,vs[j]),hgt-level*g.dot(n,es[j])) for n,hgt in hp])
                    tau=max(F(0),b-a) if ok else F(0)
                    raw+=coef*tau*(F(delta)+F(1,2)-eta[j][side])
    cluster_rows={}
    for (j,vertex,c,phase),bits in cache.items():
        if bits:cluster_rows.setdefault((j,vertex,bits),set()).add(c)
    cluster_rows={key:sorted(values) for key,values in cluster_rows.items()}
    Gamma=F(0);nonzero=0
    for vertex,clusters in catalog.items():
        for C in sorted(clusters):
            delta=coords[vertex][((C&-C).bit_length()-1)//symbol_bits]
            for i,vi in enumerate(vs):
                T=q*lcm(*(abs(g.det(vi,vj)) for j,vj in enumerate(vs) if j!=i))
                ia,ib=g.det(vi,vertex)-highs[i],g.det(vi,vertex)-lows[i]
                density=F(0)
                for r in range(T):
                    branches=[]
                    for j,vj in enumerate(vs):
                        if j==i or (j,vertex) not in sideinfo:continue
                        bulk=1-sideinfo[j,vertex];D=g.det(vi,vj)
                        for c in cluster_rows.get((j,vertex,C),()):
                            w=(F(r*vj[0]-c*vi[0],D),F(r*vj[1]-c*vi[1],D))
                            if any(z.denominator!=1 for z in w):continue
                            phase=int(g.det(w,es[j]))%q
                            if cache[j,vertex,c,phase]!=C:continue
                            others=[k for k in range(m) if k not in (i,j)]
                            for sides in product((0,1),repeat=len(others)):
                                smap=dict(zip(others,sides));smap[j]=bulk;ineq=[]
                                for k,side in smap.items():
                                    if k==j:continue
                                    slope=F(g.det(vs[k],vj),D)
                                    if side==0:ineq.append((slope,g.det(vs[k],vertex)-highs[k]))
                                    else:ineq.append((-slope,-g.det(vs[k],vertex)+lows[k]))
                                interval=linear_interval(ineq,ia,ib)
                                if interval:
                                    sigma=tuple(smap[k] for k in range(m) if k!=i)
                                    label=labels[i,sigma,int(g.det(w,es[i]))%q]
                                    branches.append((interval,label))
                    cuts={ia,ib}
                    for interval,label in branches:cuts.update(interval)
                    for a,b in zip(sorted(cuts),sorted(cuts)[1:]):
                        x=(a+b)/2;ls=[lab for (lo,hi),lab in branches if lo<x<hi]
                        density+=(b-a)*F(len(ls)-len(set(ls)),T)
                Gamma+=density
                nonzero+=bool(density)
    return {'Lambda':area,'pure_beta':pure,'raw_beta':raw,'Gamma':Gamma,'beta':raw-Gamma,
            'inner_extremes_min':amin,'inner_extremes_max':amax,'eta':eta,'gradient_shift':h,
            'gradient_track_count':N,'collar_bounds':R,'catalog_sizes':[len(catalog[x]) for x in poly],
            'nonzero_collision_families':nonzero}


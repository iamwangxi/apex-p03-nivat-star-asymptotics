#!/usr/bin/env python3
"""Exact sparse-line coefficient construction and independent finite counts.

Reuses rational cone moments and exact minimal-period reduction from the
reviewed local code. No fit is used for the theorem's coefficient formula.
The general sparse-line formula is implemented for arbitrary many directions;
its finite tables are compared to exhaustive global pattern enumeration.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import lcm
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


import star_pattern_tools as s
import full_support_lines_audit as g
c=s.period


def cone_section(poly,v,q,d):
    """X={z in tangent cone at q: det(v,z)=d}; q is strict extreme."""
    e=g.basis(v)
    inc=[(a,b) for a,b in g.facets(poly) if g.dot(a,q)==b]
    ok,lo,hi=g.scalar_interval([(g.dot(a,v),-d*g.dot(a,e)) for a,b in inc])
    if not ok or g.ceil(lo)>g.floor(hi):
        return []
    assert lo is not None and hi is not None
    return [(j*v[0]+d*e[0],j*v[1]+d*e[1]) for j in range(g.ceil(lo),g.floor(hi)+1)]


def endpoint_profiles(poly,v,k):
    """Certified finite catalog containing every section of length <2*k."""
    low,high=min(g.det(v,q) for q in poly),max(g.det(v,q) for q in poly)
    C=g.interior_row_constants(poly,v,2*k)['C']
    out={'lower':[],'upper':[]}
    for side,h,sign in [('lower',low,1),('upper',high,-1)]:
        qs=[q for q in poly if g.det(v,q)==h]
        if len(qs)!=1:
            continue
        q=qs[0]
        for d in range(C+1):
            X=cone_section(poly,v,q,sign*d)
            if len(X)<2*k:
                out[side].append((q,sign*d,X))
    return out


def visible_loss(word,L):
    k=len(word)
    return F(sum(all(word[(r+j)%k]==0 for j in range(L)) for r in range(k)),k)


def sparse_data(vertices,vs,offsets,words,background=None):
    poly=g.hull(vertices)
    words=[list(w[:g.minimal_period(w)]) for w in words]
    ks=list(map(len,words))
    Q=background['q'] if background else 1
    def bg(z):
        return tuple(background['table'][((int(z[0])+x)%Q)*Q+(int(z[1])+y)%Q]
                     for x,y in product(range(Q),repeat=2)) if background else (0,)
    bg0=bg((0,0))
    tangent_bg_periods=[next(t for t in range(1,Q+1) if bg((t*v[0],t*v[1]))==bg0) for v in vs]
    pure_periods=[lcm(k,b) for k,b in zip(ks,tangent_bg_periods)]
    es=[g.basis(v) for v in vs]
    lows=[min(g.det(v,q) for q in poly) for v in vs]
    highs=[max(g.det(v,q) for q in poly) for v in vs]
    profiles=[endpoint_profiles(poly,v,k) for v,k in zip(vs,ks)]
    losses=[]
    catalog={q:set() for q in poly}
    for j,(word,pr) in enumerate(zip(words,profiles)):
        sides={}
        for side in ('lower','upper'):
            holes=sum(not X for q,d,X in pr[side])
            zero=F(0)
            for q,d,X in pr[side]:
                if X:
                    zero+=visible_loss(word,len(X))
                    catalog[q].update(X)
            sides[side]={'holes':holes,'zero_phase_weight':zero}
        losses.append(sides)
    raw=sum(F(k)*(hi-lo) for k,lo,hi in zip(pure_periods,lows,highs))
    area=F(0)
    for size in range(2,len(vs)+1):
        coef=(-1)**size*(size-1)
        for I in combinations(range(len(vs)),size):
            hp=[]
            for j in I:
                v=vs[j]; normal=(-v[1],v[0])
                hp.extend([(normal,-lows[j]),((-normal[0],-normal[1]),highs[j])])
            K=g.halfplane_polygon(hp)
            area+=coef*g.polygon_area(K)
            for j in I:
                v=vs[j];e=es[j]
                for side,level,delta in [('lower',-lows[j],F(offsets[j])),('upper',-highs[j],F(-offsets[j]))]:
                    ok,lo,hi=g.scalar_interval([(g.dot(a,v),b-level*g.dot(a,e)) for a,b in hp])
                    tau=max(F(0),hi-lo) if ok else F(0)
                    correction=F(1,2)-losses[j][side]['holes']-losses[j][side]['zero_phase_weight']
                    raw+=coef*tau*(delta+correction)
    Gamma=F(0)
    singleton_records=[]
    for q,points in catalog.items():
        for delta in sorted(points):
            for i,vi in enumerate(vs):
                T=lcm(*(abs(g.det(vi,vj))*lcm(ks[i],ks[j],Q) for j,vj in enumerate(vs) if i!=j))
                ia,ib=g.det(vi,q)-highs[i],g.det(vi,q)-lows[i]
                cuts={ia,ib}
                for j,vj in enumerate(vs):
                    if i==j:continue
                    Dij=g.det(vi,vj)
                    for k,vk in enumerate(vs):
                        if k in (i,j):continue
                        slope=F(g.det(vk,vj),Dij)
                        for h in (lows[k],highs[k]):
                            x=F(g.det(vk,q)-h,slope)
                            if ia<x<ib:cuts.add(x)
                intervals=list(zip(sorted(cuts),sorted(cuts)[1:]))
                density=F(0)
                for r in range(T):
                    candidates=[]
                    for j,vj in enumerate(vs):
                        if i==j:continue
                        # z=W_ij(r) lies on the physical j-line; u=z-(n*q+delta).
                        Dij=g.det(vi,vj)
                        z=((F(r*vj[0]-offsets[j]*vi[0],Dij)),
                           (F(r*vj[1]-offsets[j]*vi[1],Dij)))
                        if any(x.denominator!=1 for x in z):continue
                        d=g.det(vj,delta)
                        X=next((X for side in ('lower','upper') for qq,dd,X in profiles[j][side]
                                if qq==q and dd==d),None)
                        if X is None or delta not in X:continue
                        color=words[j][int(g.det(z,es[j]))%ks[j]]
                        if color==0:continue
                        if any(words[j][int(g.det((z[0]-delta[0]+x[0],z[1]-delta[1]+x[1]),es[j]))%ks[j]]!=0
                               for x in X if x!=delta):continue
                        label=(int(g.det(z,es[i]))%ks[i],color,bg(z))
                        candidates.append((j,label))
                    for a,b in intervals:
                        x=(a+b)/2
                        labels=[]
                        for j,label in candidates:
                            Dij=g.det(vi,vs[j])
                            if all(not(lows[k]<g.det(vs[k],q)-F(g.det(vs[k],vs[j]),Dij)*x<highs[k])
                                   for k in range(len(vs)) if k not in (i,j)):
                                labels.append(label)
                        density+=(b-a)*F(len(labels)-len(set(labels)),T)
                Gamma+=density
                if density:
                    singleton_records.append({'vertex':q,'offset':delta,'long_direction':i,'density':density})
    return {'Lambda':area,'raw_beta':raw,'Gamma':Gamma,'beta':raw-Gamma,
            'pure_periods':pure_periods,
            'losses':losses,'singleton_records':singleton_records}


def sparse_star(vs,offsets,words,background=None):
    q=lcm(*map(len,words),background['q'] if background else 1)
    fields=[]
    for v,a,w in zip(vs,offsets,words):
        e=g.basis(v)
        row=[w[g.det((x,y),e)%len(w)] for x,y in product(range(q),repeat=2)]
        z=[0]*(q*q)
        fields.append({'ell':a,'r':a,'left':z,'right':z,'rows':[row]})
    if background:
        Q=background['q']
        B=[background['table'][(x%Q)*Q+y%Q] for x,y in product(range(q),repeat=2)]
        fields[0]['left']=B
        fields[0]['right']=B
        fields[0]['rows']=[[a^b for a,b in zip(B,fields[0]['rows'][0])]]
    return s.Star(vs,q,fields,'sparse formula check')


def diagonal_exact_tables():
    tables=[[],[],[]]
    raw=[]
    # n=4*t+r. Each parity of the fine square produces a coarse rectangle.
    for r in range(4):
        A=B=C=F(0)
        for rx,ry in product((0,1),repeat=2):
            a0=(r-rx)//2;b0=(r-ry)//2
            cs=[(-1,1,2,a0),(1,-1,2,b0),(1,1,0,0),(-1,-1,4,a0+b0)]
            vs,N,E=c.stable_vertices(cs)
            assert E==1
            a,b,d=c.cone_coefficients(vs,N,E,0)
            # Count nonzero coarse patterns = both-active anchors+2*(a+b+1).
            A+=a; B+=b+8; C+=d+2*(a0+b0+1)
        C+=1
        raw.append([A,B,C])
        tables[0].append(A/16)
        tables[1].append(B/4-2*r*A/16)
        tables[2].append(C-r*B/4+r*r*A/16)
    assert tables==[[F(2)]*4,[F(8)]*4,[F(2),F(3),F(2),F(-1)]]
    periods=[c.period_word(x) for x in tables]
    assert periods==[c.period_fourier(x) for x in tables]==[1,1,4]
    return {'A_B_C_in_t':raw,'Lambda_beta_gamma':tables,'minimal_periods':periods}


def independent_sparse_count(poly,vs,offsets,words,prime,n,background=None):
    """Separate support-by-support enumeration; no transition representatives.

    Multi-line anchors are solved from pairs of nonempty geometric sections;
    pure entries are listed by their joint line/background tangent phases.
    Values are added pointwise modulo prime, so this also tests p>2.
    """
    sites=g.polygon_window(g.hull(poly),n)
    idx={z:j for j,z in enumerate(sites)}
    es=[g.basis(v) for v in vs]
    words=[tuple(w[:g.minimal_period(w)]) for w in words]
    sections=[]
    for v in vs:
        rows={}
        for z in sites:rows.setdefault(g.det(v,z),[]).append(z)
        sections.append(rows)
    Q=background['q'] if background else 1
    def B(z):
        return background['table'][(z[0]%Q)*Q+z[1]%Q] if background else 0
    periods=[]
    for v,w in zip(vs,words):
        b=next(t for t in range(1,Q+1)
               if all(B((x+t*v[0],y+t*v[1]))==B((x,y)) for x,y in product(range(Q),repeat=2)))
        periods.append(lcm(len(w),b))
    cache={}
    def component(i,h,r):
        key=(i,h,r)
        if key not in cache:
            w=words[i];e=es[i]
            cache[key]={idx[z]:w[(g.det(z,e)+r)%len(w)] for z in sections[i].get(h,[])
                        if w[(g.det(z,e)+r)%len(w)]}
        return cache[key]
    bgcache={}
    def make_key(edits,u):
        rho=(u[0]%Q,u[1]%Q)
        if rho not in bgcache:bgcache[rho]=[B((z[0]+u[0],z[1]+u[1])) for z in sites]
        vals=bgcache[rho].copy()
        for edit in edits:
            for j,val in edit.items():vals[j]=(vals[j]+val)%prime
        masks=[0]*(prime-1)
        for j,val in enumerate(vals):
            if val:masks[val-1]|=1<<j
        return tuple(masks)
    anchors=set()
    for i,j in combinations(range(len(vs)),2):
        D=g.det(vs[i],vs[j])
        for h1 in sections[i]:
            for h2 in sections[j]:
                a=offsets[i]-h1;b=offsets[j]-h2
                X=vs[j][0]*a-vs[i][0]*b;Y=vs[j][1]*a-vs[i][1]*b
                if X%D==0 and Y%D==0:anchors.add((X//D,Y//D))
    patterns=set();raw=0
    for u in anchors:
        edits=[component(i,offsets[i]-g.det(v,u),g.det(u,es[i])%len(words[i])) for i,v in enumerate(vs)]
        if sum(bool(e) for e in edits)>=2:
            patterns.add(make_key(edits,u));raw+=1
    for i,(v,e,K) in enumerate(zip(vs,es,periods)):
        for h in sections[i]:
            for t in range(K):
                u=((offsets[i]-h)*e[0]+t*v[0],(offsets[i]-h)*e[1]+t*v[1])
                patterns.add(make_key([component(i,h,g.det(u,e)%len(words[i]))],u));raw+=1
    backgrounds={make_key([],rho) for rho in product(range(Q),repeat=2)}
    patterns.update(backgrounds);raw+=len(backgrounds)
    return len(patterns),raw,raw-len(patterns)


def independent_rectangle_check():
    rows=[]
    for a,b in [(5,5),(5,6),(6,5),(6,6),(7,10),(10,7),(9,9)]:
        L=a+b+1
        N=F(L*L,2)+(F((-1)**a,2) if L%2 else 0)
        expected=int(N+2*L+1)
        actual,raw,excess=independent_sparse_count([(0,0),(a,0),(a,b),(0,b)],[(1,1),(1,-1)],[0,0],[(1,),(1,)],2,1)
        assert actual==raw==expected and excess==0,(a,b,actual,raw,expected)
        rows.append({'a':a,'b':b,'complexity':actual,'representation_excess':excess})
    return rows


def main():
    cases=[
        ('dotted perpendicular triangle',s.WINDOWS['triangle'],[(1,0),(0,1)],[0,0],[(1,0),(1,0)]),
        ('dotted diagonals square',s.WINDOWS['square'],[(1,1),(1,-1)],[0,0],[(1,0),(1,0)]),
        ('third-period oblique square',s.WINDOWS['square'],[(1,0),(1,1)],[0,0],[(1,0,0),(1,0,0)]),
        ('dotted skew triangle',s.WINDOWS['triangle'],[(1,0),(1,2)],[0,0],[(1,0),(1,0)]),
        ('three dotted lines triangle',s.WINDOWS['triangle'],[(1,0),(1,2),(1,4)],[0,1,2],[(1,0),(1,0),(1,0)]),
        ('mixed sparse hexagon',s.WINDOWS['hexagon'],[(1,0),(0,1),(1,1)],[0,1,-1],[(1,0),(1,),(1,0)]),
        ('sparse singleton collisions',s.WINDOWS['triangle'],[(1,0),(1,2),(1,4)],[0,1,2],[(1,0),(1,),(0,1)]),
        ('full cross checkerboard',s.WINDOWS['square'],[(1,0),(0,1)],[0,0],[(1,),(1,)],{'q':2,'table':[0,1,1,0]}),
        ('three lines checkerboard',s.WINDOWS['triangle'],[(1,0),(1,2),(1,4)],[0,1,2],[(1,),(1,),(1,)],{'q':2,'table':[0,1,1,0]}),
    ]
    out={'diagonal_exact':diagonal_exact_tables(),'independent_rectangles':independent_rectangle_check(),'cases':[]}
    for case in cases:
        name,poly,vs,offsets,words=case[:5]
        background=case[5] if len(case)>5 else None
        data=sparse_data(poly,vs,offsets,words,background)
        star=sparse_star(vs,offsets,words,background)
        counts=[]
        for n in [12,13,14,15,16,17,20,21,28,29,36,37,48,49]:
            exact=s.exact_count(star,poly,n)
            rem=F(exact)-data['Lambda']*n*n-data['beta']*n
            if name=='dotted diagonals square':
                assert exact==2*n*n+8*n+[2,3,2,-1][n%4]
            row={'n':n,'complexity':exact,'remainder':rem}
            if n in (12,20):
                independent,raw,excess=independent_sparse_count(poly,vs,offsets,words,2,n,background)
                assert exact==independent,(name,n,exact,independent)
                row.update({'independent_count':independent,'raw_representations':raw,'representation_excess':excess,
                            'excess_minus_Gamma_n':F(excess)-data['Gamma']*n})
            counts.append(row)
        out['cases'].append({'name':name,'vertices':poly,'directions':vs,'offsets':offsets,'words':words,
                             'background':background,
                             'computed_before_counting':data,'exact_counts':counts})
    payload=json.dumps(out,ensure_ascii=False,indent=2,default=str)+'\n'
    Path(__file__).with_suffix('.json').write_text(payload,encoding='utf-8')
    print(payload,end='')


if __name__=='__main__':
    main()

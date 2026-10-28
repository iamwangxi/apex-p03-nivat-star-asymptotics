#!/usr/bin/env python3
"""Full-support periodic-line exact checks; writes only adjacent deterministic JSON.

Geometry coefficients are computed from rational tangent-cone certificates,
not fitted from pattern counts. Pattern lists are exhaustive over Z^2.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import product, combinations
from math import gcd
from pathlib import Path
import json

OUT = Path(__file__).resolve().parent


def det(a, b):
    return a[0]*b[1]-a[1]*b[0]


def dot(a, b):
    return a[0]*b[0]+a[1]*b[1]


def ceil(q):
    q = F(q)
    return -((-q.numerator)//q.denominator)


def floor(q):
    q = F(q)
    return q.numerator//q.denominator


def lcm(a, b):
    return abs(a*b)//gcd(a, b)


def basis(v):
    a, b = abs(v[0]), abs(v[1])
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        q = old_r//r
        old_r, r = r, old_r-q*r
        old_s, s = s, old_s-q*s
        old_t, t = t, old_t-q*t
    assert old_r == 1
    x = old_s*(1 if v[0] >= 0 else -1)
    y = old_t*(1 if v[1] >= 0 else -1)
    e = (-y, x)
    assert det(v, e) == 1
    return e


def hull(points):
    points = sorted(set(tuple(map(F, p)) for p in points))
    def cross(o, a, b):
        return det((a[0]-o[0], a[1]-o[1]), (b[0]-o[0], b[1]-o[1]))
    low, high = [], []
    for p in points:
        while len(low) >= 2 and cross(low[-2], low[-1], p) <= 0:
            low.pop()
        low.append(p)
    for p in reversed(points):
        while len(high) >= 2 and cross(high[-2], high[-1], p) <= 0:
            high.pop()
        high.append(p)
    result = low[:-1]+high[:-1]
    assert len(result) >= 3
    return result


def facets(poly):
    result = []
    for p, q in zip(poly, poly[1:]+poly[:1]):
        edge = (q[0]-p[0], q[1]-p[1])
        a = (edge[1], -edge[0])
        result.append((a, dot(a, p)))
    return result


def scalar_interval(inequalities):
    """Return real s interval for A*s<=B, or infeasibility."""
    lo, hi = None, None
    for A, B in inequalities:
        if A == 0:
            if B < 0:
                return False, None, None
        elif A > 0:
            x = B/A
            hi = x if hi is None else min(hi, x)
        else:
            x = B/A
            lo = x if lo is None else max(lo, x)
    if lo is not None and hi is not None and lo > hi:
        return False, lo, hi
    return True, lo, hi


def integer_interval(inequalities):
    ok, lo, hi = scalar_interval(inequalities)
    if not ok:
        return None, lo, hi
    L = None if lo is None else ceil(lo)
    H = None if hi is None else floor(hi)
    if L is not None and H is not None and L > H:
        return None, lo, hi
    witness = L if L is not None else (H if H is not None else 0)
    return witness, lo, hi


def section_interval(fs, v, e, height, n):
    return scalar_interval([(dot(a, v), n*b-height*dot(a, e)) for a, b in fs])


def section_size(fs, v, e, height, n):
    ok, lo, hi = section_interval(fs, v, e, height, n)
    if not ok:
        return 0
    assert lo is not None and hi is not None
    return max(0, floor(hi)-ceil(lo)+1)


def interior_row_constants(poly, v, required=1):
    fs, e = facets(poly), basis(v)
    q = (sum(z[0] for z in poly)/len(poly), sum(z[1] for z in poly)/len(poly))
    low, high = min(det(v,z) for z in poly), max(det(v,z) for z in poly)
    a = det(v,q)
    ok, L, H = section_interval(fs, v, e, a, 1)
    assert ok and L is not None and H is not None and H > L
    h0 = H-L
    c = h0/max(a-low, high-a)
    C = ceil(F(required+1)/c)+2
    return {'C':C,'N':1,'c':c,'interior':q}


def cone_witness(poly, fs, v, e, vertex, residue, height, modulus):
    incident = [(a,b) for a,b in fs if dot(a,vertex)==b]
    inequalities = [(dot(a,v), residue*b-height*dot(a,e)) for a,b in incident]
    s, L, H = integer_interval(inequalities)
    if s is None:
        return {'feasible':False,'lower':None if L is None else str(L),
                'upper':None if H is None else str(H)}
    z = (s*v[0]+height*e[0], s*v[1]+height*e[1])
    k0 = 0
    for a,b in fs:
        if dot(a,vertex)==b:
            assert dot(a,z) <= residue*b
        else:
            gap = modulus*(b-dot(a,vertex))
            assert gap > 0
            k0 = max(k0, ceil((dot(a,z)-residue*b)/gap))
    k0 = max(0,k0)
    scaled_vertex = (int(modulus*vertex[0]),int(modulus*vertex[1]))
    for k in [k0,k0+1]:
        n = modulus*k+residue
        point = (k*scaled_vertex[0]+z[0],k*scaled_vertex[1]+z[1])
        assert all(dot(a,point) <= n*b for a,b in fs)
    return {'feasible':True,'witness':list(z),'stable_k':k0}


def projection_certificate(poly, v):
    poly, e = hull(poly), basis(v)
    fs = facets(poly)
    low, high = min(det(v,z) for z in poly), max(det(v,z) for z in poly)
    qlo = next(z for z in poly if det(v,z)==low)
    qhi = next(z for z in poly if det(v,z)==high)
    modulus = 1
    for z in poly:
        for a in z:
            modulus = lcm(modulus,a.denominator)
    row_constants = interior_row_constants(poly,v)
    C = row_constants['C']
    records, max_k = [], 0
    for r in range(modulus):
        holes_low, holes_high, witnesses = [], [], []
        for side, vertex in [('lower',qlo),('upper',qhi)]:
            for h in range(C+1):
                t = ceil(r*low)+h if side=='lower' else floor(r*high)-h
                record = cone_witness(poly,fs,v,e,vertex,r,t,modulus)
                record.update({'side':side,'offset':h,'height':t})
                witnesses.append(record)
                if record['feasible']:
                    max_k = max(max_k,record['stable_k'])
                elif side=='lower':holes_low.append(h)
                else:holes_high.append(h)
        b = floor(r*high)-ceil(r*low)+1-r*(high-low)-len(holes_low)-len(holes_high)
        records.append({'residue':r,'holes_lower':holes_low,'holes_upper':holes_high,
                        'b':str(b),'cone_certificates':witnesses})
    N = max(2,modulus*max_k,ceil(F(2*C+4)/(high-low))+2)
    return {'direction':list(v),'basis':list(e),'low':str(low),'high':str(high),
            'width':str(high-low),'modulus':modulus,'C':C,'stable_from_n':N,
            'row_constants':{'c':str(row_constants['c'])},
            'residues':records}


def predicted_projection(poly, cert, n):
    low, high = F(cert['low']), F(cert['high'])
    a,b = ceil(n*low),floor(n*high)
    r = cert['residues'][n%cert['modulus']]
    return set(range(a,b+1))-{a+h for h in r['holes_lower']}-{b-h for h in r['holes_upper']}


def actual_projection(poly, v, n):
    fs, e = facets(poly), basis(v)
    lo, hi = min(det(v,z) for z in poly),max(det(v,z) for z in poly)
    return {t for t in range(ceil(n*lo),floor(n*hi)+1) if section_size(fs,v,e,t,n)}


def polygon_window(poly,n):
    fs=facets(poly)
    x0,x1=ceil(n*min(z[0] for z in poly)),floor(n*max(z[0] for z in poly))
    y0,y1=ceil(n*min(z[1] for z in poly)),floor(n*max(z[1] for z in poly))
    return [(x,y) for x in range(x0,x1+1) for y in range(y0,y1+1)
            if all(dot(a,(x,y)) <= n*b for a,b in fs)]


def minimal_period(word):
    for k in range(1,len(word)+1):
        if len(word)%k==0 and all(word[i]==word[(i+k)%len(word)] for i in range(len(word))):
            return k
    raise AssertionError('no period')


def global_pair_count(poly,vs,offsets,words,prime,n):
    words=[tuple(w[:minimal_period(w)]) for w in words]
    assert all(all(0<a<prime for a in w) for w in words)
    kappas=list(map(len,words)); es=[basis(v) for v in vs]
    window=polygon_window(poly,n); index={z:i for i,z in enumerate(window)}
    symbol_bits=(prime-1).bit_length(); field_mask=(1<<symbol_bits)-1
    patterns, lengths=[],[]
    for v,e,w in zip(vs,es,words):
        table=defaultdict(lambda:[0]*len(w)); sizes=Counter()
        for bit,z in enumerate(window):
            t,s=det(v,z),det(z,e);sizes[t]+=1
            for phase in range(len(w)):
                table[t][phase] |= w[(phase+s)%len(w)] << (symbol_bits*bit)
        patterns.append(dict(table));lengths.append(sizes)
    d=det(vs[0],vs[1]);D=abs(d)
    znum=(offsets[0]*vs[1][0]-offsets[1]*vs[0][0],
          offsets[0]*vs[1][1]-offsets[1]*vs[0][1])
    center=None
    if znum[0]%d==znum[1]%d==0:
        center=(znum[0]//d,znum[1]//d)
    center_value=None if center is None else sum(w[det(center,e)%len(w)] for w,e in zip(words,es))%prime
    def pattern(u):
        answer=0
        for i in range(2):
            t=offsets[i]-det(vs[i],u);phase=det(u,es[i])%kappas[i]
            answer |= patterns[i].get(t,[0]*kappas[i])[phase]
        if center is not None:
            rel=(center[0]-u[0],center[1]-u[1])
            if rel in index:
                shift=symbol_bits*index[rel]
                answer=(answer & ~(field_mask<<shift)) | (center_value<<shift)
        return answer
    multi={}
    for t1 in patterns[0]:
        A=offsets[0]-t1
        for t2 in patterns[1]:
            B=offsets[1]-t2
            nums=(A*vs[1][0]-B*vs[0][0],A*vs[1][1]-B*vs[0][1])
            if nums[0]%d or nums[1]%d:continue
            u=(nums[0]//d,nums[1]//d)
            multi[u]=(t1,t2)
    counts=Counter(); exception=Counter(); good=[]
    for u,ts in multi.items():
        p=pattern(u);counts[p]+=1
        if any(lengths[i][ts[i]] >= 3*kappas[i] for i in range(2)):
            good.append((u,p))
        else:exception[p]+=1
    short_limit=3*sum(kappas)
    for i in range(2):
        for t,row_patterns in patterns[i].items():
            for p in row_patterns:
                counts[p]+=1
                if lengths[i][t] <= short_limit:exception[p]+=1
    counts[0]+=1;exception[0]+=1
    assert all(counts[p]==1 for u,p in good)
    excess=sum(counts.values())-len(counts)
    exception_excess=sum(exception.values())-len(exception)
    assert excess==exception_excess
    assert set(counts)-set(exception) == {p for p,c in counts.items() if p not in exception and c==1}
    # Every pure phase in the representation list is globally realizable.
    for i in range(2):
        j=1-i
        for t,phase_patterns in patterns[i].items():
            for phase,p in enumerate(phase_patterns):
                u0=((offsets[i]-t)*es[i][0]+phase*vs[i][0],
                    (offsets[i]-t)*es[i][1]+phase*vs[i][1])
                s=0
                while offsets[j]-det(vs[j],u0) in patterns[j]:
                    s+=1
                    u0=(u0[0]+kappas[i]*vs[i][0],u0[1]+kappas[i]*vs[i][1])
                    assert s <= len(patterns[j])+1
                assert pattern(u0)==p
    return {'n':n,'window_size':len(window),'patterns':len(counts),
            'both_line_anchors':len(multi),'pure_representations':sum(k*len(t) for k,t in zip(kappas,patterns)),
            'representation_excess':excess,'exception_representations':sum(exception.values()),
            'exception_patterns':len(exception),'long_multi_globally_unique':len(good),
            'projection_sizes':[len(t) for t in patterns],'kappas':kappas,
            'intersection_lattice_point':None if center is None else list(center),
            'intersection_color':center_value}


def family_case(name,vertices,vs,offsets,words,prime,ns):
    poly=hull(vertices);certs=[projection_certificate(poly,v) for v in vs]
    modulus=certs[0]['modulus'];assert certs[1]['modulus']==modulus
    widths=[F(c['width']) for c in certs];D=abs(det(*vs))
    kappas=[minimal_period(w) for w in words]
    leading=widths[0]*widths[1]/D
    betas=[]
    for r in range(modulus):
        b1,b2=[F(c['residues'][r]['b']) for c in certs]
        betas.append((widths[0]*b2+widths[1]*b1)/D+kappas[0]*widths[0]+kappas[1]*widths[1])
    geometry_checks=[]
    N=max(c['stable_from_n'] for c in certs)
    for n in range(N,N+2*modulus+2):
        for v,cert in zip(vs,certs):
            actual=actual_projection(poly,v,n)
            assert actual==predicted_projection(poly,cert,n)
            b=F(cert['residues'][n%modulus]['b'])
            assert len(actual)==n*F(cert['width'])+b
        geometry_checks.append(n)
    records=[]
    for n in ns:
        rec=global_pair_count(poly,vs,offsets,words,prime,n)
        rec['residue']=n%modulus
        rec['predicted_beta']=str(betas[n%modulus])
        rec['bounded_remainder']=str(F(rec['patterns'])-leading*n*n-betas[n%modulus]*n)
        records.append(rec)
    return {'name':name,'vertices':[[str(a) for a in z] for z in poly],
            'directions':vs,'offsets':offsets,'words':words,'prime':prime,
            'leading':str(leading),'beta_by_residue':list(map(str,betas)),
            'modulus':modulus,'projection_certificates':certs,
            'projection_checks_from_proved_stability':geometry_checks,'global_checks':records}


def halfplane_polygon(inequalities):
    points=[]
    for (a,b),(c,d) in combinations(inequalities,2):
        den=det(a,c)
        if den==0:continue
        z=((b*c[1]-a[1]*d)/den,(a[0]*d-b*c[0])/den)
        if all(dot(e,z)<=f for e,f in inequalities):points.append(z)
    return hull(points)


def polygon_area(poly):
    return abs(sum(det(a,b) for a,b in zip(poly,poly[1:]+poly[:1])))/2


def multi_geometry(poly,vs,offsets):
    """Geometric coefficients by finite inclusion-exclusion, no fitted counts."""
    modulus=1
    for z in poly:
        for a in z:modulus=lcm(modulus,a.denominator)
    lows=[min(det(v,z) for z in poly) for v in vs]
    highs=[max(det(v,z) for z in poly) for v in vs]
    subsets=[];leading=F(0);betas=[F(0)]*modulus
    for k in range(2,len(vs)+1):
        for inds in combinations(range(len(vs)),k):
            coefficient=(-1)**k*(k-1)
            planes=[];labels=[]
            for i in inds:
                eta=(-vs[i][1],vs[i][0])
                planes.extend([(eta,-lows[i]),(tuple(-a for a in eta),highs[i])])
                labels.extend([(i,'upper'),(i,'lower')])
            region=halfplane_polygon(planes)
            area=polygon_area(region);leading+=coefficient*area
            edges=[];subset_beta=[F(0)]*modulus
            for x,y in zip(region,region[1:]+region[:1]):
                midpoint=((x[0]+y[0])/2,(x[1]+y[1])/2)
                matches=[j for j,(a,b) in enumerate(planes) if dot(a,midpoint)==b]
                assert len(matches)==1
                i,side=labels[matches[0]];v=vs[i]
                lattice_length=abs((y[0]-x[0])/v[0] if v[0] else (y[1]-x[1])/v[1])
                deltas=[]
                for r in range(modulus):
                    delta=offsets[i]+r*lows[i]-ceil(r*lows[i]) if side=='upper' else -offsets[i]+floor(r*highs[i])-r*highs[i]
                    subset_beta[r]+=(delta+F(1,2))*lattice_length
                    deltas.append(str(delta))
                edges.append({'direction_index':i,'side':side,'lattice_length':str(lattice_length),'delta_by_residue':deltas})
            for r in range(modulus):betas[r]+=coefficient*subset_beta[r]
            subsets.append({'indices':inds,'coefficient':coefficient,'area':str(area),'edges':edges,'anchor_beta_by_residue':list(map(str,subset_beta))})
    return {'modulus':modulus,'leading':leading,'anchor_betas':betas,'subsets':subsets,
            'lows':lows,'highs':highs}


def global_multi_count(poly,vs,offsets,words,prime,n):
    """Exhaustive global patterns in the class with long extreme sections."""
    m=len(vs);words=[tuple(w[:minimal_period(w)]) for w in words]
    kappas=list(map(len,words));es=[basis(v) for v in vs]
    window=polygon_window(poly,n);index={z:i for i,z in enumerate(window)}
    bits=(prime-1).bit_length();mask=(1<<bits)-1
    tables=[];sizes=[]
    for v,e,w in zip(vs,es,words):
        table=defaultdict(lambda:[0]*len(w));size=Counter()
        for bit,z in enumerate(window):
            t,s=det(v,z),det(z,e);size[t]+=1
            for phase in range(len(w)):
                table[t][phase] |= w[(s+phase)%len(w)] << (bits*bit)
        assert min(size.values()) >= max(2*m,len(w))
        tables.append(dict(table));sizes.append(size)
    anchors=set();centers=set()
    for i,j in combinations(range(m),2):
        d=det(vs[i],vs[j])
        def solve(A,B):
            nums=(A*vs[j][0]-B*vs[i][0],A*vs[j][1]-B*vs[i][1])
            return None if nums[0]%d or nums[1]%d else (nums[0]//d,nums[1]//d)
        center=solve(offsets[i],offsets[j])
        if center is not None:centers.add(center)
        for t,s in product(tables[i],tables[j]):
            u=solve(offsets[i]-t,offsets[j]-s)
            if u is not None:anchors.add(u)
    center_colors={z:sum(words[i][det(z,es[i])%kappas[i]] for i in range(m) if det(vs[i],z)==offsets[i])%prime for z in centers}
    def pattern(u):
        p=0
        for i in range(m):
            t=offsets[i]-det(vs[i],u);phase=det(u,es[i])%kappas[i]
            p |= tables[i].get(t,[0]*kappas[i])[phase]
        for z,color in center_colors.items():
            rel=(z[0]-u[0],z[1]-u[1])
            if rel in index:
                shift=bits*index[rel];p=(p & ~(mask<<shift)) | (color<<shift)
        return p
    patterns={pattern(u) for u in anchors};assert len(patterns)==len(anchors)
    pure_count=sum(k*len(t) for k,t in zip(kappas,tables))
    for table in tables:
        for phases in table.values():
            for p in phases:
                assert p not in patterns;patterns.add(p)
    assert 0 not in patterns;patterns.add(0)
    assert len(patterns)==len(anchors)+pure_count+1
    # Finite inclusion-exclusion of all geometric active-line intersections.
    ie=0;subset_counts=[]
    for k in range(2,m+1):
        for inds in combinations(range(m),k):
            count=sum(all(offsets[i]-det(vs[i],u) in tables[i] for i in inds) for u in anchors)
            ie+=(-1)**k*(k-1)*count
            subset_counts.append({'indices':inds,'anchors':count})
    assert ie==len(anchors)
    return {'n':n,'window_size':len(window),'patterns':len(patterns),'multi_anchors':len(anchors),
            'pure_representations':pure_count,'projection_sizes':[len(t) for t in tables],
            'minimum_section_sizes':[min(s.values()) for s in sizes],
            'zero_pattern':1,'all_representations_distinct':True,'subset_counts':subset_counts}


def multi_case(name,vs,offsets,words,prime,scale,ns):
    vertices=[(F(0),F(0))]
    for v in vs:
        vertices=[(z[0]+s*scale*v[0],z[1]+s*scale*v[1]) for z in vertices for s in [-1,1]]
    poly=hull(vertices);geometry=multi_geometry(poly,vs,offsets)
    kappas=[minimal_period(w) for w in words]
    pure_beta=sum(k*(b-a) for k,a,b in zip(kappas,geometry['lows'],geometry['highs']))
    betas=[b+pure_beta for b in geometry['anchor_betas']]
    records=[]
    for n in ns:
        rec=global_multi_count(poly,vs,offsets,words,prime,n)
        rec['predicted_beta']=str(betas[n%geometry['modulus']])
        rec['bounded_remainder']=str(F(rec['patterns'])-geometry['leading']*n*n-betas[n%geometry['modulus']]*n)
        records.append(rec)
    return {'name':name,'directions':vs,'offsets':offsets,'words':words,'prime':prime,
            'vertices':[[str(a) for a in z] for z in poly],'modulus':geometry['modulus'],
            'leading':str(geometry['leading']),'beta_by_residue':list(map(str,betas)),
            'geometric_subsets':geometry['subsets'],'global_checks':records}


def singleton_stability_certificate(poly,vs,certs):
    """One/two cone witnesses certify every eventual singleton section."""
    fs=facets(poly);M=certs[0]['modulus'];records=[];max_k=0
    for j,(v,cert) in enumerate(zip(vs,certs)):
        e=basis(v)
        for r in range(M):
            for side in ['lower','upper']:
                extreme=F(cert['low'] if side=='lower' else cert['high'])
                q=next(z for z in poly if det(v,z)==extreme)
                incident=[(a,b) for a,b in fs if dot(a,q)==b]
                for h in range(cert['C']+1):
                    t=ceil(r*extreme)+h if side=='lower' else floor(r*extreme)-h
                    witness,L,H=integer_interval([(dot(a,v),r*b-t*dot(a,e)) for a,b in incident])
                    ss=[]
                    if witness is not None:
                        ss=[witness]
                        if H is None or witness+1<=floor(H):ss.append(witness+1)
                        elif L is None or witness-1>=ceil(L):ss.append(witness-1)
                    zs=[];k0=0
                    for s in ss:
                        z=(s*v[0]+t*e[0],s*v[1]+t*e[1]);zs.append(z)
                        for a,b in fs:
                            if dot(a,q)==b:assert dot(a,z)<=r*b
                            else:k0=max(k0,ceil((dot(a,z)-r*b)/(M*(b-dot(a,q)))))
                    max_k=max(max_k,k0)
                    records.append({'direction':j,'residue':r,'side':side,'offset':h,'height':t,
                                    'vertex':[str(x) for x in q],'witnesses':zs,'stable_k':max(0,k0),
                                    'singleton':len(ss)==1})
    N=max(max(c['stable_from_n'] for c in certs),M*max_k,1)
    return {'stable_from_n':N,'witness_records':records}


def singleton_vertex_types(poly,vs,certs):
    """Exact eventual singleton positions, in each integral vertex translate."""
    fs=facets(poly);M=certs[0]['modulus'];records=[]
    for q in poly:
        for r in range(M):
            types=defaultdict(set)
            for i,(v,cert) in enumerate(zip(vs,certs)):
                e=basis(v);height=det(v,q)
                lo,hi=F(cert['low']),F(cert['high'])
                if height==lo and sum(det(v,z)==lo for z in poly)==1:side='lower'
                elif height==hi and sum(det(v,z)==hi for z in poly)==1:side='upper'
                else:continue
                incident=[(a,b) for a,b in fs if dot(a,q)==b]
                for h in range(cert['C']+1):
                    t=ceil(r*lo)+h if side=='lower' else floor(r*hi)-h
                    witness,L,H=integer_interval([(dot(a,v),r*b-t*dot(a,e)) for a,b in incident])
                    if witness is None:continue
                    assert L is not None and H is not None
                    if floor(H)-ceil(L)+1!=1:continue
                    s=ceil(L);z=(s*v[0]+t*e[0],s*v[1]+t*e[1])
                    types[z].add(i)
            for z,J in sorted(types.items()):
                records.append({'vertex':q,'residue':r,'offset':z,'directions':sorted(J)})
    return records


def singleton_collision_density(q,J,i,vs,offsets,words,lows,highs):
    """Finite interval/phase integral for one singleton-position type."""
    m=len(vs);kappas=list(map(minimal_period,words));es=[basis(v) for v in vs]
    T=1
    for j in range(m):
        if j!=i:T=lcm(T,abs(det(vs[i],vs[j]))*lcm(kappas[i],kappas[j]))
    left=det(vs[i],q)-highs[i];right=det(vs[i],q)-lows[i]
    endpoints={left,right};J=[j for j in J if j!=i]
    if len(J)<2:return {'density':'0','phase_modulus':T,'interval_cells':[]}
    for j in J:
        dij=det(vs[i],vs[j])
        for k in range(m):
            if k in [i,j]:continue
            slope=F(det(vs[k],vs[j]),dij);assert slope
            for h in [lows[k],highs[k]]:
                x=(det(vs[k],q)-h)/slope
                if left<x<right:endpoints.add(x)
    endpoints=sorted(endpoints);total=F(0);cells=[]
    for a,b in zip(endpoints,endpoints[1:]):
        mid=(a+b)/2
        geometrically_valid=[]
        for j in J:
            dij=det(vs[i],vs[j])
            if all(not lows[k] < det(vs[k],q)-F(det(vs[k],vs[j]),dij)*mid < highs[k]
                   for k in range(m) if k not in [i,j]):geometrically_valid.append(j)
        sum_excess=0;phase_records=[]
        for t in range(T):
            labels=[];candidates=[]
            for j in geometrically_valid:
                d=det(vs[i],vs[j])
                w=(F(t*vs[j][0]-offsets[j]*vs[i][0],d),F(t*vs[j][1]-offsets[j]*vs[i][1],d))
                if any(x.denominator!=1 for x in w):continue
                w=tuple(map(int,w))
                label=(det(w,es[i])%kappas[i],words[j][det(w,es[j])%kappas[j]])
                labels.append(label);candidates.append(j)
            excess=len(labels)-len(set(labels));sum_excess+=excess
            if excess:phase_records.append({'t_residue':t,'directions':candidates,'labels':labels,'excess':excess})
        total+=(b-a)*F(sum_excess,T)
        cells.append({'left':str(a),'right':str(b),'geometrically_valid':geometrically_valid,
                      'sum_excess_over_phases':sum_excess,'nonzero_phases':phase_records})
    return {'density':str(total),'phase_modulus':T,'interval_cells':cells}


def global_case_B_count(poly,vs,offsets,words,prime,n,allow_intersections=False):
    """All global representations, including all finite-pattern multiplicities."""
    m=len(vs);kappas=list(map(minimal_period,words));es=[basis(v) for v in vs]
    window=polygon_window(poly,n);index={z:i for i,z in enumerate(window)}
    bits=(prime-1).bit_length();mask=(1<<bits)-1
    tables=[];sizes=[];positions=[]
    for v,e,w in zip(vs,es,words):
        table=defaultdict(lambda:[0]*minimal_period(w));size=Counter();points=defaultdict(list)
        for bit,z in enumerate(window):
            t,s=det(v,z),det(z,e);size[t]+=1;points[t].append(z)
            for phase in range(minimal_period(w)):
                table[t][phase] |= w[(s+phase)%minimal_period(w)] << (bits*bit)
        tables.append(dict(table));sizes.append(size);positions.append(points)
    anchors=set();centers=set()
    for i,j in combinations(range(m),2):
        d=det(vs[i],vs[j])
        nums=(offsets[i]*vs[j][0]-offsets[j]*vs[i][0],offsets[i]*vs[j][1]-offsets[j]*vs[i][1])
        if nums[0]%d==nums[1]%d==0:
            assert allow_intersections,'not Case B full support'
            centers.add((nums[0]//d,nums[1]//d))
        for t,s in product(tables[i],tables[j]):
            A,B=offsets[i]-t,offsets[j]-s
            nums=(A*vs[j][0]-B*vs[i][0],A*vs[j][1]-B*vs[i][1])
            if nums[0]%d or nums[1]%d:continue
            anchors.add((nums[0]//d,nums[1]//d))
    counts=Counter();template_entries=defaultdict(list);other_entries=Counter();bad_anchors=0;bad_patterns=set()
    center_colors={z:sum(words[i][det(z,es[i])%kappas[i]] for i in range(m) if det(vs[i],z)==offsets[i])%prime for z in centers}
    threshold=[max(2*m if allow_intersections else m,k) for k in kappas]
    for u in anchors:
        active=[];p=0
        for i in range(m):
            t=offsets[i]-det(vs[i],u)
            if t in tables[i]:
                active.append((i,t,sizes[i][t]));p |= tables[i][t][det(u,es[i])%kappas[i]]
        short_ids={i for i,t,L in active if L<threshold[i]}
        overlap_with_short=False
        for z,color in center_colors.items():
            rel=(z[0]-u[0],z[1]-u[1])
            if rel in index:
                shift=bits*index[rel];p=(p & ~(mask<<shift)) | (color<<shift)
                if any(det(vs[i],z)==offsets[i] for i in short_ids):overlap_with_short=True
        counts[p]+=1
        long=[x for x in active if x[2]>=threshold[x[0]]]
        if len(active)-len(long)>=2 or overlap_with_short:
            bad_anchors+=1;bad_patterns.add(p)
        if len(active)==2 and len(long)==1:
            i=long[0][0];j,t,L=next(x for x in active if x[0]!=i)
            if L==1:
                q=positions[j][t][0];w=(u[0]+q[0],u[1]+q[1])
                if det(vs[i],w)!=offsets[i]:
                    key=(i,q,det(vs[i],w),det(u,es[i])%kappas[i],words[j][det(w,es[j])%kappas[j]])
                    template_entries[p].append(key)
                else:other_entries[p]+=1
            else:other_entries[p]+=1
        else:other_entries[p]+=1
    H=sum(threshold)+m
    for i,table in enumerate(tables):
        for t,phases in table.items():
            for p in phases:
                counts[p]+=1;other_entries[p]+=1
                if sizes[i][t]<=H:bad_patterns.add(p)
    counts[0]+=1;other_entries[0]+=1;bad_patterns.add(0)
    classified_excess=0
    for p,keys in template_entries.items():
        assert len(set(keys))==1,'singleton template mismatch'
        classified_excess+=len(keys)-1
    excess=sum(counts.values())-len(counts)
    assert classified_excess<=excess
    closure_size=sum(counts[p] for p in bad_patterns)
    closure_excess=closure_size-len(bad_patterns)
    family_excess_outside=sum(len(keys)-1 for p,keys in template_entries.items() if p not in bad_patterns)
    assert all(counts[p]==len(template_entries[p]) for p,c in counts.items() if c>1 and p not in bad_patterns)
    assert sum(counts.values())-closure_excess-family_excess_outside==len(counts),'exact formula (55) failed'
    # Independent enumeration of the finite one-parameter families, using
    # singleton sites and physical line intersections rather than pair anchors.
    singleton_sites=defaultdict(set)
    for j in range(m):
        for t,ps in positions[j].items():
            if len(ps)==1:singleton_sites[ps[0]].add(j)
    direct_family_excess=0
    for q,J in singleton_sites.items():
        if len(J)<2:continue
        for i in range(m):
            for height,L in sizes[i].items():
                if L<threshold[i]:continue
                t=offsets[i]+det(vs[i],q)-height;labels=[]
                if t==offsets[i]:continue
                for j in J-{i}:
                    d=det(vs[i],vs[j])
                    w=(F(t*vs[j][0]-offsets[j]*vs[i][0],d),F(t*vs[j][1]-offsets[j]*vs[i][1],d))
                    if any(x.denominator!=1 for x in w):continue
                    w=tuple(map(int,w));u=(w[0]-q[0],w[1]-q[1])
                    if any(offsets[k]-det(vs[k],u) in tables[k] for k in range(m) if k not in [i,j]):continue
                    labels.append((det(w,es[i])%kappas[i],words[j][det(w,es[j])%kappas[j]]))
                direct_family_excess+=len(labels)-len(set(labels))
    assert direct_family_excess==classified_excess
    exceptional_reps=bad_anchors+sum(kappas[i] for i in range(m) for t,L in sizes[i].items() if L<=H)+1
    assert excess-classified_excess <= (2*m+sum(kappas)+1)*exceptional_reps
    return {'n':n,'window_size':len(window),'patterns':len(counts),'multi_anchors':len(anchors),
            'projection_sizes':[len(t) for t in tables],'raw_representations':sum(counts.values()),
            'total_excess':excess,'singleton_family_excess':classified_excess,
            'other_excess':excess-classified_excess,
            'exceptional_representations':exceptional_reps,
            'bounded_collision_closure_size':closure_size,'closure_patterns':len(bad_patterns),
            'singleton_excess_outside_closure':family_excess_outside,'exact_formula_55_check':'passed',
            'independent_singleton_family_count':direct_family_excess,
            'singleton_families':sum(len(v)>1 for v in template_entries.values()),
            'singleton_key_checks':'passed'}


def large_singleton_family_check(poly,vs,offsets,words,n,certs,allow_intersections=False,singleton_certificate=None,singleton_types=None):
    """O(n) exact family count without storing a window or pattern bitsets."""
    m=len(vs);es=[basis(v) for v in vs];kappas=list(map(minimal_period,words));fs=facets(poly)
    sizes=[];singletons=defaultdict(set)
    for j,(v,e,cert) in enumerate(zip(vs,es,certs)):
        table={}
        for t in range(ceil(n*F(cert['low'])),floor(n*F(cert['high']))+1):
            ok,L,H=section_interval(fs,v,e,t,n)
            if not ok:continue
            length=max(0,floor(H)-ceil(L)+1)
            if not length:continue
            table[t]=length
            if length==1:
                s=ceil(L);q=(s*v[0]+t*e[0],s*v[1]+t*e[1]);singletons[q].add(j)
        if n>=cert['stable_from_n']:assert set(table)==predicted_projection(poly,cert,n)
        sizes.append(table)
    if singleton_certificate is not None and n>=singleton_certificate['stable_from_n']:
        expected=defaultdict(set);M=certs[0]['modulus'];k,r=divmod(n,M)
        for typ in singleton_types:
            if typ['residue']!=r:continue
            q,z=typ['vertex'],typ['offset']
            site=(int(k*M*q[0])+z[0],int(k*M*q[1])+z[1])
            expected[site].update(typ['directions'])
        assert dict(expected)==dict(singletons),'singleton cone classification incomplete'
    total=0;nonzero_classes=0
    for q,J in singletons.items():
        if len(J)<2:continue
        for i in range(m):
            for height,length in sizes[i].items():
                if length<max(2*m if allow_intersections else m,kappas[i]):continue
                t=offsets[i]+det(vs[i],q)-height;labels=[]
                if t==offsets[i]:continue
                for j in J-{i}:
                    d=det(vs[i],vs[j])
                    w=(F(t*vs[j][0]-offsets[j]*vs[i][0],d),F(t*vs[j][1]-offsets[j]*vs[i][1],d))
                    if any(x.denominator!=1 for x in w):continue
                    w=tuple(map(int,w));u=(w[0]-q[0],w[1]-q[1])
                    if any(offsets[k]-det(vs[k],u) in sizes[k] for k in range(m) if k not in [i,j]):continue
                    labels.append((det(w,es[i])%kappas[i],words[j][det(w,es[j])%kappas[j]]))
                e=len(labels)-len(set(labels));total+=e;nonzero_classes+=bool(e)
    return {'n':n,'exact_singleton_excess':total,'nonzero_classes':nonzero_classes,
            'common_singleton_sites':sum(len(J)>1 for J in singletons.values())}


def case_B_case(name,vertices,vs,offsets,words,prime,ns,allow_intersections=False):
    poly=hull(vertices);certs=[projection_certificate(poly,v) for v in vs]
    geometry=multi_geometry(poly,vs,offsets);M=geometry['modulus']
    base=geometry['anchor_betas'][:];kappas=list(map(minimal_period,words))
    for r in range(M):
        for subset in geometry['subsets']:
            for edge in subset['edges']:
                i=edge['direction_index'];record=certs[i]['residues'][r]
                holes=len(record['holes_lower' if edge['side']=='upper' else 'holes_upper'])
                base[r]-=subset['coefficient']*F(edge['lattice_length'])*holes
        base[r]+=sum(k*(b-a) for k,a,b in zip(kappas,geometry['lows'],geometry['highs']))
    types=singleton_vertex_types(poly,vs,certs)
    singleton_cert=singleton_stability_certificate(poly,vs,certs)
    gammas=[F(0)]*M;density_records=[]
    for typ in types:
        if len(typ['directions'])<2:continue
        for i in range(len(vs)):
            density=singleton_collision_density(typ['vertex'],typ['directions'],i,vs,offsets,words,geometry['lows'],geometry['highs'])
            gammas[typ['residue']]+=F(density['density'])
            density_records.append({'vertex':[str(x) for x in typ['vertex']],
                                    'residue':typ['residue'],'offset':typ['offset'],
                                    'singleton_directions':typ['directions'],'long_direction':i,**density})
    betas=[a-b for a,b in zip(base,gammas)];records=[]
    for n in ns:
        rec=global_case_B_count(poly,vs,offsets,words,prime,n,allow_intersections)
        rec['predicted_beta']=str(betas[n%M]);rec['predicted_collision_density']=str(gammas[n%M])
        rec['bounded_remainder']=str(F(rec['patterns'])-geometry['leading']*n*n-betas[n%M]*n)
        rec['collision_remainder']=str(F(rec['singleton_family_excess'])-gammas[n%M]*n)
        records.append(rec)
    large_checks=[]
    for n in sorted({max(100,singleton_cert['stable_from_n'])+r for r in range(M)}|{200+r for r in range(M)}|{500+r for r in range(M)}):
        rec=large_singleton_family_check(poly,vs,offsets,words,n,certs,allow_intersections,singleton_cert,types)
        rec['predicted_density']=str(gammas[n%M])
        rec['bounded_collision_remainder']=str(F(rec['exact_singleton_excess'])-n*gammas[n%M])
        large_checks.append(rec)
    return {'name':name,'vertices':[[str(a) for a in z] for z in poly],'directions':vs,'offsets':offsets,
            'words':words,'prime':prime,'allow_intersections':allow_intersections,'modulus':M,'leading':str(geometry['leading']),
            'base_beta_by_residue':list(map(str,base)),'collision_density_by_residue':list(map(str,gammas)),
            'beta_by_residue':list(map(str,betas)),'projection_certificates':certs,
            'singleton_stability_certificate':singleton_cert,
            'geometric_subsets':geometry['subsets'],'singleton_type_integrals':density_records,'global_checks':records,
            'large_exact_singleton_checks':large_checks}


def counterexample_check(n):
    # Sparse orthogonal Case B star, old configuration but new rational window.
    M=n//2+1;N=n+1
    window=[(x,y) for x in range(M) for y in range(N)]
    hm,vm=defaultdict(int),defaultdict(int)
    for bit,(x,y) in enumerate(window):
        hm[(y,x%2)] |= 1<<bit
        vm[(x,y%2)] |= 1<<bit
    xs=set(range(2-M,2));ys=set(range(1-N,1))
    def p(x,y):return hm.get((-y,x%2),0)|vm.get((1-x,y%2),0)
    ps={p(x,y) for x in xs for y in ys}
    for y in ys:
        for phase in range(2):ps.add(p(4*M+phase,y))
    for x in xs:
        for phase in range(2):ps.add(p(x,4*N+phase))
    ps.add(0)
    exact=M*N+2*M+2*N+1;r=n%2
    claimed=F(n*n,2)+(F(9,2)-F(r,2))*n+6-F(3*r,2)
    assert M>=5 and N>=5 and len(ps)==exact==claimed
    return {'n':n,'patterns':len(ps),'leading':'1/2',
            'beta':str(F(9,2)-F(r,2)),'constant':str(6-F(3*r,2))}


def main():
    cases=[
        family_case('integral triangle, determinant 2, disjoint lines',
                    [(0,0),(1,0),(0,1)],[(1,0),(1,2)],[0,1],[[1],[1,2]],3,[20,30,40,60]),
        family_case('integral triangle, six missing extreme layers',
                    [(0,0),(1,0),(0,1)],[(1,2),(1,4)],[1,2],[[1],[1]],2,[20,30,40,60]),
        family_case('rational rectangle, full intersecting lines',
                    [(0,0),(F(1,2),0),(F(1,2),1),(0,1)],[(1,0),(0,1)],[0,0],[[1],[1]],3,[10,11,20,21,40,41]),
        family_case('integral triangle, periods two/three and cancellation',
                    [(0,0),(1,0),(0,1)],[(1,0),(0,1)],[0,0],[[1,4],[4,2,3]],5,[20,30,40,60]),
        family_case('rational trapezoid, determinant three',
                    [(0,0),(F(5,3),0),(1,F(4,3)),(0,1)],[(1,2),(2,1)],[0,0],[[1,2],[3]],5,[24,25,26,36,37,38]),
        family_case('rational rectangle, determinant two Case B',
                    [(0,0),(F(1,2),0),(F(1,2),1),(0,1)],[(1,0),(1,2)],[0,1],[[1],[1]],2,[20,21,30,31,40,41]),
    ]
    for rec in cases[2]['global_checks']:
        n=rec['n'];assert rec['patterns']==(n//2+2)*(n+2)
    assert cases[0]['beta_by_residue']==['17/2']
    assert cases[1]['beta_by_residue']==['1/2']
    assert cases[2]['beta_by_residue']==['3','5/2']
    assert cases[3]['beta_by_residue']==['7']
    assert cases[5]['beta_by_residue']==['9/2','4']
    multi_cases=[
        multi_case('three directions, integral hexagon, periodic colors and triple point',
                   [(1,0),(0,1),(1,1)],[0,0,0],[[1],[1,4],[2,3,4]],5,F(1),[5,6,9,10]),
        multi_case('three directions, rational third-hexagon',
                   [(1,0),(0,1),(1,1)],[0,1,2],[[1],[1,2],[3]],5,F(1,3),[9,10,11,12]),
        multi_case('four directions, nonunimodular octagon',
                   [(1,0),(1,2),(0,1),(2,1)],[0,1,-1,2],[[1],[1,2],[2],[3,4]],5,F(1),[4,5,6,7]),
    ]
    case_B_cases=[
        case_B_case('old three-line triangle recovered by general singleton integral',
                    [(0,0),(1,0),(0,1)],[(1,0),(1,2),(1,4)],[0,1,2],[[1],[1],[1]],2,[20,21,30,31,40,41,60]),
        case_B_case('three-line triangle with periodic middle line',
                    [(0,0),(1,0),(0,1)],[(1,0),(1,2),(1,4)],[0,1,2],[[1],[1,2],[1]],3,[20,21,30,31,40,41,60]),
        case_B_case('three-line rational half-triangle',
                    [(0,0),(F(1,2),0),(0,F(1,2))],[(1,0),(1,2),(1,4)],[0,1,2],[[1],[1],[1]],2,[30,31,40,41,60,61,80,81]),
        case_B_case('four nonconcurrent disjoint lines, periodic colors, triangle',
                    [(0,0),(1,0),(0,1)],[(1,0),(1,2),(1,4),(1,6)],[0,1,2,7],[[1],[1,2],[2],[2,1,1]],3,[20,21,30,31,40,41,60]),
        case_B_case('four disjoint lines, equal colors and multiple reassignments',
                    [(0,0),(1,0),(0,1)],[(1,0),(1,2),(1,4),(1,6)],[0,1,2,7],[[1],[1],[1],[1]],2,[30,31,40,41,60]),
        case_B_case('integer-translated triangle and reversed primitive orientations',
                    [(-2,3),(-1,3),(-2,4)],[(-1,0),(-1,-2),(1,4)],[0,-1,2],[[1],[1],[1]],2,[30,31,40,41]),
    ]
    assert case_B_cases[0]['leading']=='31/4'
    assert case_B_cases[0]['beta_by_residue']==['11/4']
    assert case_B_cases[0]['collision_density_by_residue']==['1/4']
    assert case_B_cases[-1]['leading']=='31/4'
    assert case_B_cases[-1]['beta_by_residue']==['11/4']
    assert [r['patterns'] for r in case_B_cases[-1]['global_checks']]==[7037,7512,12487,13120]
    intersecting_cases=[
        case_B_case('three intersecting lines, arbitrary triangle, binary cancellation',
                    [(0,0),(1,0),(0,1)],[(1,0),(1,2),(1,4)],[0,0,2],[[1],[1],[1]],2,[30,31,40,41,60],True),
        case_B_case('four concurrent lines, triangle, periods and total intersection cancellation',
                    [(0,0),(1,0),(0,1)],[(1,0),(0,1),(1,1),(1,2)],[0,0,0,0],[[1],[1,2],[1],[2,3,4]],5,[20,21,30,31,40,41,60],True),
        case_B_case('three intersecting lines, rational trapezoid',
                    [(0,0),(F(5,3),0),(1,F(4,3)),(0,1)],[(1,0),(1,2),(2,1)],[0,0,1],[[1,2],[3],[1]],5,[24,25,26,36,37,38],True),
    ]
    result={'version':'v1.2','source_manuscript_sha256':'6a0974713cb34194d5892b98284c20c5c7eae56326c69466b48565dcd8ba9a7f',
            'scope':'arbitrary full-support nonzero periodically colored affine lines and arbitrary rational polygons; not all CN stars',
            'rational_case_B_counterexample':[counterexample_check(n) for n in [8,9,20,21,40,41]],
            'cases':cases,'long_extreme_section_cases':multi_cases,'case_B_arbitrary_polygon_cases':case_B_cases,
            'intersecting_arbitrary_polygon_cases':intersecting_cases,'all_assertions':'passed'}
    target=OUT/'full_support_lines_audit.json'
    if target.exists():
        previous=json.loads(target.read_text());assert previous.get('version')=='v1.2'
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'all_assertions':'passed','counterexample_betas':['9/2','4'],
                      'families':[{'name':x['name'],'leading':x['leading'],'betas':x['beta_by_residue'],
                                   'N_projection':max((c['stable_from_n'] for c in x.get('projection_certificates',[])),default=None),
                                   'checks':[(r['n'],r['patterns'],r['bounded_remainder']) for r in x['global_checks']]}
                                  for x in cases+multi_cases+case_B_cases+intersecting_cases]},ensure_ascii=False,indent=2))


if __name__=='__main__':main()

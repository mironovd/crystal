from sage_import import sage_import
sage_import('crystal_op')

import numpy as np

def mutate_node(I,js,C,w,node,badroot,t,accumulator,sinks,stats):
    
    if not ('type_two_mutation' in stats.keys()):
        stats['type_two_mutation']=[]

    m = [node.degree(tt) for tt in t]
        
    result = []

    b=None
    if node in accumulator:
        b=accumulator[node]
    else:
        b=crystal_op.bn(C,w,js,t,node)
        accumulator[node]=b

#    print(b)
    if len([kk for kk in b[1:] if kk < 0])>0:
        print("[WARN] Found negative b_i", b,node)
    for l in I:
        for k in crystal_op.r(w,l):
            if crystal_op.plus(w,l,k)<oo:
                if m[k]>0 and b[crystal_op.plus(w,l,k)]>0:
                    if m[k]>m[crystal_op.plus(w,l,k)]:
                        nn = node*crystal_op.Ainv(C,w,l,k,t)
                        result.append([l, nn, k,k])
                        bn=b.copy()
                        bn[k]+=1
                        bn[crystal_op.plus(w,l,k)]-=1
                        accumulator[nn]=bn
                    if m[k]==m[crystal_op.plus(w,l,k)]:
                        mm=2
                        while crystal_op.pplus(w,l,k,mm)<oo and m[crystal_op.pplus(w,l,k,mm)]==0 and b[crystal_op.pplus(w,l,k,mm)]==0:
                            mm+=1
                        if mm!=None:
                            if crystal_op.pplus(w,l,k,mm)<oo and m[crystal_op.pplus(w,l,k,mm)]==-1 and b[crystal_op.pplus(w,l,k,mm)]==1:
                                nn = node*crystal_op.Ainv(C,w,l,k,t)
                                result.append([l, nn, k,k])
                                bn=b.copy()
                                bn[k]+=1
                                bn[crystal_op.plus(w,l,k)]-=1
                                accumulator[nn]=bn
                                stats['type_two_mutation'].append([node,nn,mm])

    if len(result)==0:
        sinks.append(node)
    return result


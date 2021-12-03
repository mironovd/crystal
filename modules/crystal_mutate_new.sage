from sage_import import sage_import
sage_import('crystal_op')

import numpy as np

def mutate_node(I,js,C,w,node,badroot,t,accumulator):
    
    m = [node.degree(tt) for tt in t]
        
    result = []
    for l in I:
        b=None
        if t in accumulator:
            b=accumulator[t]
        else:
            b=crystal_op.bn(C,w,js,t)
            accumulator[t]=b

        for k in crystal_op.r(w,l):
            if crystal_op.plus(w,l,k)<oo:
                if m[k]>0 and b[crystal_op.plus(w,l,k)]>0:
                    if m[k]>m[crystal_op.plus(w,l,k)]:
                        tn = node*crystal_op.Ainv(C,w,l,k,t)
                        result.append([l, tn, k,k])
                        bn=b.copy()
                        bn[k]-=1
                        bn[crystal_op.plus(w,l,k)]+=1
                        accumulator[tn]=bn
                    if m[k]==m[crystal_op.plus(w,l,k)]:
                        mm=2
                        while crystal_op.pplus(w,l,k,mm)<oo and m[crystal_op.pplus(w,l,k,mm)]==0 and b[crystal_op.pplus(w,l,k,mm)]==0:
                            mm+=1
                        if mm!=None:
                            if crystal_op.pplus(w,l,k,mm)<oo and m[crystal_op.pplus(w,l,k,mm)]==-1 and b[crystal_op.pplus(w,l,k,mm)]==1:
                                tn = node*crystal_op.Ainv(C,w,l,k,t)
                                result.append([l, tn, k,k])
                                bn=b.copy()
                                bn[k]-=1
                                bn[crystal_op.plus(w,l,k)]+=1
                                accumulator[tn]=bn

    return result


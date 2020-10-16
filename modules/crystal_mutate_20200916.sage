from sage_import import sage_import
sage_import('crystal_op')

def mutate_node(I,C,w,node,badroot,t):
    m = [node.degree(tt) for tt in t]
    N=len(w)
    result = []
    for l in I:
        for k in crystal_op.r(w,l): #range(1,N+1):
            if m[k]==1 and crystal_op.plus(w,l,k)<oo and m[crystal_op.plus(w,l,k)]==0:
                s=[]
                kmax=k
                tt=1
                while crystal_op.pplus(w,l,k,tt)<oo:
                    kmax=crystal_op.pplus(w,l,k,tt)
                    tt+=1
                for ss in range(k,kmax+1):
                    if C[w[ss-1]-1][w[k-1]-1] < 0  and m[ss]==-1 :
                        s.append(ss)
                if len(s)>0:
                    for ss in s:
                        le=0
                        for lle in range(k+1,ss):
                            if m[lle]!=0 and C[w[lle-1]-1][w[k-1]-1]!=0:
                                le=lle
                        if le!=0:
                            result.append([l,node*crystal_op.Ainv(C,w,l,k,t), k, k])
                        if le==0:
                            if badroot or (crystal_op.minus(w,w[ss-1],ss)>0 and crystal_op.minus(w,w[ss-1],ss)<k) :
                                tt = min([ttt for ttt in range(ss+1,N+1) if w[ttt-1]==w[k-1]])
#                                print("t ",tt, " t- ", minus(w,l,tt))
                                result.append([l,node*crystal_op.Ainv(C,w,l,crystal_op.minus(w,l,tt),t), k, crystal_op.minus(w,l,tt)])
                else:
                    result.append([l,node*crystal_op.Ainv(C,w,l,k,t),k,k])
            elif m[k] == 1 and crystal_op.plus(w,l,k)<oo and m[crystal_op.plus(w,l,k)] == 1 and crystal_op.pplus(w,l,k,2)<oo and m[crystal_op.pplus(w,l,k,2)]== -1:
                result.append([l, node*crystal_op.Ainv(C,w,l,k,t),k,k])
                result.append([l, node*crystal_op.Ainv(C,w,l,crystal_op.plus(w,l,k),t), k, crystal_op.plus(w,l,k)])
            elif m[k] == 2 and crystal_op.plus(w,l,k)<oo and m[crystal_op.plus(w,l,k)]== 1:
                tt=2
                while crystal_op.pplus(w,l,k,tt)<oo:
                    if m[crystal_op.pplus(w,l,k,tt)]==-1:
                        result.append([l, node*crystal_op.Ainv(C,w,l,crystal_op.pplus(w,l,k,tt-1),t),k,crystal_op.pplus(w,l,k,tt-1)])
                    tt+=1
            elif m[k] == 2 and crystal_op.plus(w,l,k)<oo and m[crystal_op.plus(w,l,k)]==0:
                result.append([l, node*crystal_op.Ainv(C,w,l,k,t),k,k ])
    return result


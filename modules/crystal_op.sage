
def pad(some_list, target_len):
     return some_list[:target_len] + [0]*(target_len - len(some_list))

def roots_from_reduced_word(w, ct):
     al = RootSystem(ct).root_lattice().simple_roots()
     ret = []
     for i in range(len(w)):
         ret.append(al[w[i]].weyl_action(w[:i]))
     return ret


def r(w,jj):
    return [k+1 for k,i in enumerate(w) if i==jj]


def plus(w,j,k):
    return next((x for x in r(w,j) if x>k), oo)

def minus(w,j,k):
#    print(w,j,k)
#    print(r(w,j))
    m=[x for x in r(w,j) if plus(w,j,x)==k]
#    print(m)
    if len(m)>0:
        return m[0]
    else:
        return -1
#    return [x for x in r(w,j) if plus(w,j,x)==k][0]

def pplus(w,j,k,power):
    if power>0:
        ret=plus(w,j,k)
        for i in range(power-1):
            ret=plus(w,j,ret)
    else:
        ret=k
    return ret


def RR(w,m,j,t):
    explist = [m.degree(tt) for tt in t]
    return [k for k in r(w,j) if explist[k]>0 ]


def Ainv(C,w,j,k,t):
    end=plus(w,j,k)
    tp=1
    if end == oo:
        end = len(w)+1
        tp=1
    else:
        tp=t[end]
    return (prod([ t[m]^(-C[w[m-1]-1][j-1])   for m in range(k+1,end) ]))/(t[k]*tp)



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


def bn(C,ww,jj,t,node):


    R=C.root_system()

    j=jj
    w=ww

    if R.cartan_type()[0]=='E' and R.cartan_type()[1]==6:
        if j == 1 :
            j = 6
        elif j == 6 :
            j = 1
        elif j == 3 :
            j = 5
        elif j == 5 :
            j = 3
    

    if  R.cartan_type()[0]=='A':
        j = R.cartan_type()[1]-j+1
    if  R.cartan_type()[0]=='D' and is_odd(R.cartan_type()[1]):
        if j == R.cartan_type()[1] :
            j = R.cartan_type()[1]-1
        elif j == R.cartan_type()[1] - 1 :
            j = R.cartan_type()[1]
#    if R.cartan_type()[0]=='C':
#        w=list(reversed(w))
#        if j == R.cartan_type()[1]:
#            j = R.cartan_type()[1]-1
#        elif j == R.cartan_type()[1]-1:
#            j = R.cartan_type()[1]
#        j=R.cartan_type()[1]-j+1



    P=R.weight_lattice()
    h = P.simple_coroots()
    La = P.fundamental_weights()
    m = [node.degree(tt) for tt in t]
    N=len(w)

    b=[None for iiii in range(N+1)]
    b[N]=m[N]+La[j].weyl_action([j]).scalar(h[w[N-1]])
    for tt in range(N-1,0,-1):
        b[tt]=m[tt]
        b[tt]+=(La[j].weyl_action([j]).scalar(h[ w[tt-1] ]))
        for ll in range(tt,N):
            b[tt]-=b[ll+1]*C[w[(ll+1)-1]-1][w[tt-1]-1]

    return b

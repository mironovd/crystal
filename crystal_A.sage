from numpy import prod

import argparse

parser=argparse.ArgumentParser(description='Clrystal calculations for A_n', prog='sage cl.sage')
parser.add_argument('degree', metavar='deg', type=int,
        help='rank of A_n algebra')
parser.add_argument('-i', '--inversions', metavar='inversions', dest='inversions',type=int, default=-1,
        help='process only decompositions with set number of inversions')
parser.add_argument('-n', '--num', metavar='randdec', dest='num',action='store', type=int, default=0,
        help='process only N random decompositions')

args = parser.parse_args()
#print(args.num)

def mergeSortInversions(arr):
    if len(arr) == 1:
        return arr, 0
    else:
        a = arr[:len(arr)//2]
        b = arr[len(arr)//2:]
        a, ai = mergeSortInversions(a)
        b, bi = mergeSortInversions(b)
        c = []
        i = 0
        j = 0
        inversions = 0 + ai + bi
    while i < len(a) and j < len(b):
        if a[i] <= b[j]:
            c.append(a[i])
            i += 1
        else:
            c.append(b[j])
            j += 1
            inversions += (len(a)-i)
    c += a[i:]
    c += b[j:]
    return c, inversions

def Inversions(arr):
    b,c=mergeSortInversions(arr)
    return c

type=["A",args.degree]
#js=3

def pad(some_list, target_len):
     return some_list[:target_len] + [0]*(target_len - len(some_list))

def roots_from_reduced_word(w, ct):
     al = RootSystem(ct).root_lattice().simple_roots()
     ret = []
     for i in range(len(w)):
         ret.append(al[w[i]].weyl_action(w[:i]))
     return ret

W=WeylGroup(type)
al=RootSystem(type).root_lattice().simple_roots()
I=[i.leading_support() for i in al]
w=W.long_element(as_word=True)

#print(w)
C=CartanMatrix(type)
N=len(w)

ww=W.long_element().reduced_words()
ww.sort()

if args.inversions>=0:
    ww=[ws[0] for ws in [[ws,Inversions(ws)] for ws in ww] if ws[1]==args.inversions]
    
if args.num>0:
    shuffle(ww)
    ww=ww[0:args.num-1]



R=LaurentPolynomialRing(QQ,'t',N+1)

t=list(var(R.variable_names()))
#t.insert(0,1)
#print(t)

#start=t[[p for p in [[i+1,x] for i,x in enumerate(roots_from_reduced_word(w,type)) if len(x.monomials())==1] if p[1]==al[js]][0][0]]

def r(w,jj):
    return [k+1 for k,i in enumerate(w) if i==jj]


def plus(w,j,k):
    return next((x for x in r(w,j) if x>k), oo)

def minus(w,j,k):
    return [x for x in r(w,j) if plus(w,j,x)==k][0]

def pplus(w,j,k,power):
    if power>0:
        ret=plus(w,j,k)
        for i in range(power-1):
            ret=plus(w,j,ret)
    else:
        ret=k
    return ret


def RR(w,m,j):
    explist = [m.degree(tt) for tt in t]
    return [k for k in r(w,j) if explist[k]>0 ]


def Ainv(w,j,k):
    end=plus(w,j,k)
    tp=1
    if end == oo:
        end = N+1
        tp=1
    else:
        tp=t[end]
    return (prod([ t[m]^(-C[w[m-1]-1][j-1])   for m in range(k+1,end) ]))/(t[k]*tp)


def mutate_node(w,node):
    poss= [[x[0],x[1][0]] for x in [ [l,[k for k in RR(w,node,l) if not (plus(w,l,k) in RR(w,node,l))]] for l in I] if len(x[1])>0]
    return [[x[0], node*Ainv(w,x[0],x[1]) ]  for x in poss]


def Step(n,wc):
    St.append([])
    nodes=St[n]
    res=[]
    new_nodes = [[node,mutate_node(wc,node)] for node in nodes]
#    print(new_nodes)
    for r in new_nodes:
        node=r[0]
        for new_node in r[1]:
            if not (new_node[1] in res):
                res.append(new_node[1])
                print(new_node)
#            print(node)
#            print(new_node)
            Gs[(str(node),str(new_node[0]))]=r[1][0][1]
    St[n+1]=res

#xx=true
#st=0

#print(St)

qt=1

for wx in ww:
    print("\n==== Case ",qt," ====\n")
    print("Long element decomposition: ",wx,"\n")
    qt+=1
    for js in I:
        print("==== Begin ====")
        print("Simple root: ",js)
        start=t[[p for p in [[i+1,x] for i,x in enumerate(roots_from_reduced_word(wx,type)) if len(x.monomials())==1] if p[1]==al[js]][0][0]]
        xx= true
        st = 0
        St=[[start]]
        Gs={}
        print(St)
        while xx:
            print("===\n Step: ",str(st))
            Step(st,wx)
#            print(St[st+1])
            if len(St[st+1])==0:
               xx=false
            st+=1

        print(str(St))
        print(str(Gs))
        print("==== End ====")

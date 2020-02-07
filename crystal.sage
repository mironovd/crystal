from numpy import prod

import argparse

parser=argparse.ArgumentParser(description='Crystal calculations for A_n', prog='sage crystal_A.sage')

required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('-t', '--type', metavar='type', dest='type', default="",
        help='type of algebra', required=True)
required.add_argument('-r', '--rank', metavar='rank', dest='rank', type=int, default="",
        help='rank of algebra',  required=True)
optional.add_argument('-i', '--inversions', metavar='inversions', dest='inversions',type=int, default=-1,
        help='process only decompositions with set number of inversions')
optional.add_argument('-n', '--num', metavar='N', dest='num',action='store', type=int, default=0,
        help='process only N random decompositions')

parser.add_argument('word',type=int,metavar='word',default=[],nargs='*',
        help='start calculations from word (needs to be reduced word corresponding to longest element of Weyl group)')

args = parser.parse_args()
#print(args.num)

print(args.word)

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




def BraidOrbit(word,rels):
 def pattern_match (L, i, X, l):
  for ind in range(l):
   if L[i+ind] != X[ind]:
    return False
  return True
  
 l=len(word)
 words = set(tuple(word))
 test_words = [ tuple(word) ]

 rels = rels + [ [b,a] for a,b in rels ]
 rels = [ [tuple(a), tuple(b), len(a) ]  for a,b in rels ]

 loop_ind = 0
 list_len = 1
 yield word
 while loop_ind < list_len:
  test_word = test_words[loop_ind]
  loop_ind += 1
  for rel in rels:
   left = rel[0]
   right = rel[1]
   rel_l = rel[2]
   for i in range(l-rel_l+1):
    if pattern_match(test_word, i, left, rel_l):
     new_word=test_word[:i]+right+test_word[i+rel_l:]
     if new_word not in words:
      words.add(new_word)
      test_words.append(new_word)
      list_len+=1
      yield new_word
 return

def coxeter_braid_orbit(coxeter_group, word):
 word=list(word)
 braid_rels=coxeter_group.braid_relations()
 return BraidOrbit(word, braid_rels)


type=[args.type,args.rank]
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
w=[]
if len(args.word)>0:
 w=args.word
else:
 w=W.long_element(as_word=True)


#print(w)
C=CartanMatrix(type)
N=len(w)

#ww=W.long_element().reduced_words()
#ww.sort()

#if args.inversions>=0:
#    ww=[ws[0] for ws in [[ws,Inversions(ws)] for ws in ww] if ws[1]==args.inversions]
#    
#if args.num>0:
#    shuffle(ww)
#    ww=ww[0:args.num-1]

def generator(wg,word):
 gen=coxeter_braid_orbit(wg,word)
 n=0
 for el in gen:
  y=True
  if args.inversions>=0 and Inversions(el)!=args.inversions:
   y=False
  if args.num>0 and n>=args.num:
   break
  if y:
   n+=1
   yield el
 return


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

for wx in generator(W,w):
    print("\n==== Case ",qt," ====\n")
    print("Long element decomposition: ",wx,"\n")
    qt+=1
    for js in I:
        print("==== Begin ====")
        print("Simple root: ",js)
#        if len([q for q in C[js-1] if q!=0])>3:
#             print("not implemented yet")
#             continue
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
        print("Number of monomials: ",sum([len(z) for z in St]))
        print(str(Gs))
        print("==== End ====")

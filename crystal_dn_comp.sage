from numpy import prod
import os
import argparse
import glob
import re

from sage.symbolic.expression_conversions import laurent_polynomial

from sage_import import sage_import
sys.modules['sage_import'].__dict__['extendpath_rel']('modules')

sage_import('braid_orbit')
sage_import('crystal_op')
sage_import('crystal_step')
sage_import('crystal_generate')

parser=argparse.ArgumentParser(description='Crystal calculations for D_n - comparison of combinatorial method and determinant method', prog='sage crystal_dn_comp.sage')

required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('-m', '--method', metavar='mutation_method', dest='mutation_method', default="",
        help='mutation method, one of: '+(", ".join(sorted([re.sub(r'^crystal_mutate_(.*)\.sage',r'\1',os.path.basename(i)) for i  in glob.glob( os.path.dirname(os.path.realpath(__file__))+'/modules/crystal_mutate_*.sage' ) ] ))), 
        required=True)


#required.add_argument('-t', '--type', metavar='type', dest='type', default="",
#        help='type of algebra', required=True)
required.add_argument('-r', '--rank', metavar='rank', dest='rank', type=int, default="",
        help='rank of algebra',  required=True)
required.add_argument('-k', '--root', metavar='root', dest='root', type=int, default=0,
        help='root of algebra')
optional.add_argument('-i', '--inversions', metavar='inversions', dest='inversions',type=int, default=-1,
        help='process only decompositions with set number of inversions')
optional.add_argument('-n', '--num', metavar='N', dest='num',action='store', type=int, default=0,
        help='process only N random decompositions')



parser.add_argument('word',type=int,metavar='word',default=[],nargs='*',
        help='start calculations from word (needs to be reduced word corresponding to longest element of Weyl group)')

args = parser.parse_args()


type=['D',args.rank]

W=WeylGroup(type)
al=RootSystem(type).root_lattice().simple_roots()
I=[i.leading_support() for i in al]
IN=I 
if args.root>0 :
     IN=[args.root]
# can't work for last 2 roots
IN=[kk for kk in IN if kk<(args.rank-1)]

w=[]
if len(args.word)>0:
 w=args.word
else:
 w=W.long_element(as_word=True)

C=CartanMatrix(type)
N=len(w)

R=LaurentPolynomialRing(QQ,'t',N+1)

t=list(var(R.variable_names()))

qt=1


def mutation_import(mod):
    sage_import('crystal_mutate_'+mod,None,None,'mutate_crystal')
    return mutate_crystal.mutate_node

mutate_method = mutation_import(args.mutation_method)

def is_badroot(C,j):
    return False
    


def f(w,r,k,x,t):
     i=w[k-1]
     res=[]
     if i<r:
         res=tuple([0]+[x[j] for j in range(1,(i-1)+1)]+
                   [(1/t[k])*x[i],t[k]*x[i+1]+x[i]] + 
                   [x[j] for j in range(i+2,(2*r-i-1)+1)] + 
                   [(1/t[k])*x[2*r-i],t[k]*x[2*r-i+1]+x[2*r-i] ]+
                   [x[j] for j in range(2*r-i+2,2*r+1)])
     else:
         res=tuple([0]+[x[j] for j in range(1,(r-2)+1)]+
                   [(1/t[k])*x[r-1],(1/t[k])*x[r],x[r-1]+t[k]*x[r+1],x[r]+t[k]*x[r+2]] +
                   [x[j] for j in range(r+3,2*r+1)])
     return res

def v(w,r,j,t):
    N=len(w)
    z=tuple([1 if jj==j else 0 for jj in range(0,2*r+1)])
    for kk in reversed(range(1,N+1)):
        z=f(w,r,kk,z,t)
    return z
    

def vmat(w,r,t):
    return Matrix([v(w,r,j,t)[1:(2*r+1)] for j in range(1,r+1)])
    
    
def explist(m,t):
    def tindex(tt,t):
        return [str(tk) for tk in t].index(str(tt))
    ex=dict([(tindex(tt,t),m.degree(tt)) for tt in m.variables()])
#    print(ex)
    return [ ex[tt] if tt in ex.keys() else 0 for tt in range(len(t)) ]

         
for wx in braid_orbit.braid_orbit_generator(W,w,args.num,args.inversions):
     print("\n==== Case ",qt," ====\n")
     print("Long element decomposition: ",wx,"\n")
     qt+=1
     vm=vmat(wx,args.rank,t)
     for js in IN:
         print("==== Begin ====")
         print("Simple root: ",js)

         start=t[[p for p in [[i+1,x] for i,x in enumerate(crystal_op.roots_from_reduced_word(wx,type)) if len(x.monomials())==1] if p[1]==al[js]][0][0]]
         St,Gs = crystal_generate.generate_crystal(I,js,C,wx,start,mutate_method,is_badroot,t,False)
         nm = sum([len(z) for z in St])
         print("Number of monomials by method: ",nm)

         subvm=vm[[j for j in range(0,js-1)]+[js],[j for j in range((2*args.rank-js+1)-1,2*args.rank)]]
         ans=laurent_polynomial(subvm.det(),R)
#         print("Answer:\n",ans,"\n",St)
         nd = ans.number_of_terms()
#         nd = len(ans.monomials())
         print("Number of monomials by det: ", nd)
         print("Number of monomials equals: ",nm==(nd+1))
         dexp=[ explist(m,t) for m in ans.monomials()]
         mexp=[[m.degree(tt) for tt in t] for m in [z for zz in St for z in zz] if not m==1]
         diffdm = [prod([t[i]^p for i,p in enumerate(x)]) for x in dexp if not x in mexp]
         diffmd = [prod([t[i]^p for i,p in enumerate(x)]) for x in mexp if not x in dexp]
         print("Monomials in det, not in method: ", diffdm )
         print("Monomials in method, not in det: ", diffmd )
         print("Sets of monomials equal: ",diffdm==[] and diffmd==[])
         print("==== End ====")
#         import code; code.interact(local=dict(globals(), **locals()))

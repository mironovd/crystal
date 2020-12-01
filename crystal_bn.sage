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

parser=argparse.ArgumentParser(description='Crystal calculations for D_n', prog='sage crystal.sage')

required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

#required.add_argument('-m', '--method', metavar='mutation_method', dest='mutation_method', default="",
#        help='mutation method, one of: '+(", ".join(sorted([re.sub(r'^crystal_mutate_(.*)\.sage',r'\1',os.path.basename(i)) for i  in glob.glob( os.path.dirname(os.path.realpath(__file__))+'/modules/crystal_mutate_*.sage' ) ] ))), 
#        required=True)


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


type=['B',args.rank]

W=WeylGroup(type)
al=RootSystem(type).root_lattice().simple_roots()
I=[i.leading_support() for i in al]
IN=I 
if args.root>0 :
     IN=[args.root]
# can't work for last 2 roots
IN=[kk for kk in IN if kk<(args.rank)]

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



def f(w,r,k,x,t):
     i=w[k-1]
     res=[]
     if i<r:
         res=tuple([0]+[x[j] for j in range(1,(i-1)+1)]+
                   [(1/t[k])*x[i],t[k]*x[i+1]+x[i]] + 
                   [x[j] for j in range(i+2,(2*r-i)+1)] + 
                   [(1/t[k])*x[2*r-i+1],t[k]*x[2*r-i+2]+x[2*r-i+1] ]+
                   [x[j] for j in range(2*r-i+3,(2*r+1)+1)])
     else:
         res=tuple([0]+[x[j] for j in range(1,(r-1)+1)]+
                   [(1/(t[k]^2))*x[r], x[r+1]+(1/t[k])*x[r], x[r]+2*t[k]*x[r+1]+(t[k]^2)*x[r+2] ] +
                   [x[j] for j in range(r+3,(2*r+1)+1)])
     return res

def v(w,r,j,t):
    N=len(w)
    z=tuple([1 if jj==j else 0 for jj in range(0,2*r+2)])
    for kk in reversed(range(1,N+1)):
        z=f(w,r,kk,z,t)
    return z
    

def vmat(w,r,t):
    return Matrix([v(w,r,j,t)[1:(2*r+2)] for j in range(1,2*r+1)])
    
         
for wx in braid_orbit.braid_orbit_generator(W,w,args.num,args.inversions):
     print("\n==== Case ",qt," ====\n")
     print("Long element decomposition: ",wx,"\n")
     qt+=1
     vm=vmat(wx,args.rank,t)
     for js in IN:
         print("==== Begin ====")
         print("Simple root: ",js)
         subvm=vm[[j for j in range(0,js-1)]+[js],[j for j in range((2*args.rank-js+2)-1,2*args.rank+1)]]
         ans=laurent_polynomial(subvm.det(),R)
         print("Answer:\n",ans)
         print("Number of monomials: ", len(ans.monomials()))
         print("==== End ====")


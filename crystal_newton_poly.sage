from numpy import prod
import os
import argparse
import glob
import re

from itertools import chain

from sage_import import sage_import
sys.modules['sage_import'].__dict__['extendpath_rel']('modules')

sage_import('braid_orbit')
sage_import('crystal_op')
sage_import('crystal_step')
sage_import('crystal_generate')

parser=argparse.ArgumentParser(description='Crystal calculations for A_n', prog='sage crystal.sage')

required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('-m', '--method', metavar='mutation_method', dest='mutation_method', default="",
        help='mutation method, one of: '+(", ".join(sorted([re.sub(r'^crystal_mutate_(.*)\.sage',r'\1',os.path.basename(i)) for i  in glob.glob( os.path.dirname(os.path.realpath(__file__))+'/modules/crystal_mutate_*.sage' ) ] ))), 
        required=True)


required.add_argument('-t', '--type', metavar='type', dest='type', default="",
        help='type of algebra', required=True)
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


type=[args.type,args.rank]

W=WeylGroup(type)
al=RootSystem(type).root_lattice().simple_roots()
I=[i.leading_support() for i in al]
IN=I 
if args.root>0 :
     IN=[args.root]
w=[]
if len(args.word)>0:
 w=args.word
else:
 w=W.long_element(as_word=True)

C=CartanMatrix(type)
N=len(w)

R=LaurentPolynomialRing(QQ,'t',N+1)

t=list(var(R.variable_names()))

exec((",".join(["t"+str(tt) for tt in t[1:]]))+' = polygen(QQ,'+"'"+(",".join(["t"+str(tt) for tt in t[1:]]))+"'"+')')

qt=1


def mutation_import(mod):
    sage_import('crystal_mutate_'+mod,None,None,'mutate_crystal')
    return mutate_crystal.mutate_node

mutate_method = mutation_import(args.mutation_method)

def is_badroot(C,j):
#    len([q for q in C[j-1] if q!=0])>3
    sum(C[j-1])<0

def is_badroot_a(C,j):
    if C.cartan_type()[0]=='A':
       return False
    if C.cartan_type()[0]=='G' or C.cartan_type()[0]=='F':
       return False
    if C.cartan_type()[0]=='B':
       return not(j==C.cartan_type[1])
    if C.cartan_type()[0]=='C':
       return not(j==1)
    if C.cartan_type()[0]=='D':
       return not(j>=(C.cartan_type()[1]-1) or j==1)
    if C.cartan_type()[0]=='E':
       if C.cartan_type()[1]>=8:
          return False
       if C.cartan_type()[1] == 6:
          return not(j==1 or j==6)
       if C.cartan_type()[1] == 7:
          return not(j==7)


         
for wx in braid_orbit.braid_orbit_generator(W,w,args.num,args.inversions):
     print("\n==== Case ",qt," ====\n")
     print("Long element decomposition: ",wx,"\n")
     qt+=1
     pt=0
     for js in IN:
         print("==== Begin ====")
         print("Simple root: ",js)
         #jsbad = is_badroot(C,js)
         start=t[[p for p in [[i+1,x] for i,x in enumerate(crystal_op.roots_from_reduced_word(wx,type)) if len(x.monomials())==1] if p[1]==al[js]][0][0]]
         print("Start: ",start)
         St,Gs,sinks = crystal_generate.generate_crystal(I,js,C,wx,start,mutate_method,is_badroot,t,True)
#         print(str(St))
         print("mutation done")
         print(*St, sep='\n')
         print("Number of monomials: ",sum([len(z) for z in St]))
         print(str(Gs))
         print("Sinks (",len(sinks),"):")
         print(str(sinks))
         
         ans=sum(list(chain.from_iterable(St)))
         pt+=ans
         print("Half potential: ",ans)
         exec(preparse('ans_pp='+re.sub('t','tt',str(ans))))
         ans_num=(ans_pp*(ans_pp.denominator())).numerator()
         nw_pol=ans_num.newton_polytope()
         print("Vertices: ",len(nw_pol.vertices()))
         mi=[(tuple(ii.exponents()[0]), ans_num.coefficient(ii)) for ii in ans_num.monomials()]
         nw_pol_int_points = nw_pol.integral_points()

         if not any([nw_pol.interior_contains(ii) for ii in nw_pol_int_points]) :
             print("NW OK: newton polytope is empty")
         else:
             print("NW NOT OK: newton polytope contains integral interior points")

         if sorted([v[0] for v in mi]) == sorted([tuple(xx) for xx in nw_pol_int_points]):
             print("NW OK: newton polytope has only points from half-potential (saturated)")
         else:
             print("NW NOT OK: newton polytope has more points outside of half-potential")

         print("==== End ====")

     print("Full potential: ",pt)
     exec(preparse('ans_pp='+re.sub('t','tt',str(pt))))
     ans_num=(ans_pp*(ans_pp.denominator())).numerator()
     nw_pol=ans_num.newton_polytope()
     print("Full potentical vertices: ",len(nw_pol.vertices()))
     mi=[(tuple(ii.exponents()[0]), ans_num.coefficient(ii)) for ii in ans_num.monomials()]
     nw_pol_int_points = nw_pol.integral_points()

     if not any([nw_pol.interior_contains(ii) for ii in nw_pol_int_points]) :
         print("NW OK: Full potentical newton polytope is empty")
     else:
         print("NW NOT OK: Full potentical newton polytope contains integral interior points")

     if sorted([v[0] for v in mi]) == sorted([tuple(xx) for xx in nw_pol_int_points]):
         print("NW OK: Full potentical newton polytope has only points from half-potential (saturated)")
     else:
         print("NW NOT OK: Full potentical newton polytope has more points outside of half-potential")


from numpy import prod
import os
import argparse
import glob
import re
import inspect

from sage_import import sage_import
sys.modules['sage_import'].__dict__['extendpath_rel']('modules')

sage_import('braid_orbit')
sage_import('crystal_op')
sage_import('crystal_step')
sage_import('crystal_generate')

parser=argparse.ArgumentParser(description='Crystal calculations for A_n', prog='sage crystal.sage')

required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('-m', '--method', metavar='mutation_methods', dest='mutation_method', default="",nargs='*',
        help='two or more mutation methods: '+(", ".join(sorted([re.sub(r'^crystal_mutate_(.*)\.sage',r'\1',os.path.basename(i)) for i  in glob.glob( os.path.dirname(os.path.realpath(__file__))+'/modules/crystal_mutate_*.sage' ) ] ))), 
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

print("Comparing mutation methods: "+", ".join(args.mutation_method))

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

qt=1

def mutation_import(mod):
    sage_import('crystal_mutate_'+mod,None,None,'mutate_crystal')
    return mutate_crystal.mutate_node

mutation_methods = [ ( mutation_import(m), m)  for m in args.mutation_method ]
#print([inspect.getfile(m) for m in mutation_methods])

flatten = lambda l: [item for sublist in l for item in sublist]

uniq2d = lambda l : list(map(list, set(map(tuple, map(set, l)))))

def is_badroot(C,j):
#    len([q for q in C[j-1] if q!=0])>3
    sum(C[j-1])<0


for wx in braid_orbit.braid_orbit_generator(W,w,args.num,args.inversions):
     print("\n==== Case ",qt," ====\n")
     print("Long element decomposition: ",wx,"\n")
     qt+=1
     for js in IN:
#         print("==== Begin ====")
         print("Simple root: ",js)
         jsbad = len([q for q in C[js-1] if q!=0])>3
         start=t[[p for p in [[i+1,x] for i,x in enumerate(crystal_op.roots_from_reduced_word(wx,type)) if len(x.monomials())==1] if p[1]==al[js]][0][0]]
         mutation_results = [ (m[1],  crystal_generate.generate_crystal(I,js,C,wx,start,m[0],is_badroot,t,False) )for m in mutation_methods ]
#         print(mutation_results)
         
         if len(set([len(flatten(s[1][0])) for s in mutation_results]))>1:
             print("Number of monomials differ")
             print(mutation_results)

         if len(uniq2d([flatten(s[1][0]) for s in mutation_results] ))>1:
             print("Set of monomials differ")
             print(mutation_results)

#         if len(uniq2d([s[1][0] for s in mutation_results] ))>1:
#             print("Step structure of monomials differ")
#             print(mutation_results)

         if len(uniq2d([s[1][1] for s in mutation_results] ))>1:
             print("Graph structure of monomials differ")
             print(mutation_results)
             

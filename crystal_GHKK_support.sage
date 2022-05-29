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
sage_import('braid_orbit_classes_walk')
sage_import('polymake')

parser=argparse.ArgumentParser(description='Crystal calculations', prog='sage crystal_GHKK_support.sage')

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
optional.add_argument('-p', '--polymake',  dest='polymake',action='store_true', default=true,
        help='use Polymake to check Newton polytopes')
optional.add_argument('-P', '--no-polymake',  dest='polymake',action='store_false',  default=true,
        help='do not use Polymake to check Newton polytopes')
optional.add_argument('-v', '--verbose',  dest='verbose',action='store_true', default=false,
        help='use Polymake to check Newton polytopes')



parser.add_argument('word',type=int,metavar='word',default=[],nargs='*',
        help='start calculations from word (needs to be reduced word corresponding to longest element of Weyl group)')

args = parser.parse_args()


type=[args.type,args.rank]

pm=None

if args.polymake:
    pm=polymake.polymake_start()

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
YR=LaurentPolynomialRing(QQ,'Y',N+1)

t=list(var(R.variable_names()))
Y=list(var(YR.variable_names()))

exec((",".join(["t"+str(tt) for tt in t[1:]]))+' = polygen(QQ,'+"'"+(",".join(["t"+str(tt) for tt in t[1:]]))+"'"+')')

qt=1

GV=VectorSpace(QQ,N)

def contains_vector(VS,vec):
    return VS.ambient_vector_space().subspace(VS.basis()+[vec]).rank()==VS.rank()

def mutation_import(mod):
    sage_import('crystal_mutate_'+mod,None,None,'mutate_crystal')
    return mutate_crystal.mutate_node

mutate_method = mutation_import(args.mutation_method)


braid_orbit_method = braid_orbit.braid_orbit_generator
if args.type in ['A','D','E']:
    braid_orbit_method = braid_orbit_classes_walk.braid_orbit_classes_walk_generator

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


def check_poly(mons_list,tv):
    mexp=[]
    if 1 in mons_list:
         mexp += [[0 for tt in tv]]
    mexp += [[mo.degree(tt) for tt in tv] for mo in mons_list if mo!=1]
    if args.verbose: print(mexp)
    (v,i,b)=polymake.check_polytope(pm,mexp)
    if args.verbose: print(v,i,b)

    print('Number of monomials: ',len(mons_list))
    print('Number of vertices: ',v)
    print('Number of interior lattice points: ',i)
    print('Number of boundary points: ',b)
    if i==0 or (v==1 and i==1 and len(mons_list)==1):
         print('Newton polytope is OK - no interior points')
         if (i==0 and b>len(mons_list)) :
             print('Newton polytope has additional vertices on boundary')
         if v<b:
             print('Some monomials do not lie on vertices')
    else:
         print('Newton polytope is NOT OK')

         
for wx in braid_orbit_method(W,w,args.num,args.inversions):
     print("\n==== Case ",qt," ====\n")
     print("Long element decomposition: ",wx,"\n")
     qt+=1
     pt=[]
     pg=[]
     for js in IN:
         print("==== Begin ====")
         print("Simple root: ",js)
         #jsbad = is_badroot(C,js)
         start_i=[p for p in [[i+1,x] for i,x in enumerate(crystal_op.roots_from_reduced_word(wx,type)) if len(x.monomials())==1] if p[1]==al[js]][0][0]
         start_wl = wx[start_i-1]
         start=t[start_i]
         print("Start ",start)
         if args.verbose: print(crystal_op.roots_from_reduced_word(wx,type))
         b_start = crystal_op.bn(C,wx,js,t,start)
         if args.verbose: print("Start b_n: ",start,b_start)
         St,Gs,sinks = crystal_generate.generate_crystal(I,js,C,wx,start,mutate_method,is_badroot,t,True)[:3]
         if args.verbose: print(str(St))
         if args.verbose: print(*St, sep='\n')
         print("Number of monomials: ",sum([len(z) for z in St]))
         if args.verbose: print(str(Gs))
         print("Sinks (",len(sinks),"):")
         print(str(sinks))


         mons=list(chain.from_iterable(St))
         ans = sum(mons)
         pt+=mons
         print("Half potential: ",ans)
         if args.verbose: print(mons)
         if args.polymake: check_poly(mons,t) 


         stop=[x[0][0] for x in Gs if x[1]==1][0]
         b_stop = crystal_op.bn(C,wx,js,t,stop)
         if args.verbose: print("Stop b_n: ",stop, b_stop)
         bb = [x-y for x,y in zip(b_stop[1:],b_start[1:]) ]
         if args.verbose: print(bb)
         ee = [ GV.basis()[y[0]-1] - GV.basis()[y[1]-1]  for y in [ y for y in [(x,crystal_op.plus(wx,wx[x-1],x)) for x in range(1,len(t))] if y[1]<+Infinity] ]
         if args.verbose: print(ee)
         GVV=GV.subspace_with_basis(ee)

#         print(GVV)

         if contains_vector(GVV,GV(bb)):
             print("Gross-Hacking-Keel-Kontsevich potential support: ", GVV.coordinates(bb))
             if len([j for j in GVV.coordinates(bb) if j<0])>0:
                print("found negative GHKK support coordinate")
             if args.verbose: print(Y[max(loc for loc, val in enumerate(wx) if val == js)+1])
             Y_start=product([Y[list(ee[i]).index(1)+1]^v for i,v in enumerate(GVV.coordinates(bb)) if v!=0])*Y[max(loc for loc, val in enumerate(wx) if val == js)+1]
             if args.verbose: print(Y_start)
             YS={}
             YS[St[0][0]]=Y_start
             YP=[Y_start]
             for edge in Gs:
                if edge[1] != 1:
                    ny=YS[edge[0][0]]
                    ny*=Y[edge[0][2]]^-1
                    YS[edge[1]]=ny
                    YP.append(ny)
             YPP=sum(list(set(YP)))
             print("GHKK potential:", YPP)
             YPI=expand((sum(YP)-Y_start)/Y_start)
             YPII=YPI.subs(dict([(y,1/y) for y in Y])) 
             if YPI==YPII :
                print("GHKK is symmetric")
             print(str(YS))
            
             Ymons=list(set(YP))
             pg+=Ymons
             if args.polymake: check_poly(Ymons,Y) 

         else:
             print("ERROR: b_start-b_stop is not contained in span of e_l-e_{l+} !!!")
         print("==== End ====")

     print('Full half-potential: ',sum(pt))
     if args.polymake: check_poly(set(pt),t) 
     print('Full GHKK potential: ',sum(pg))
     if args.polymake: check_poly(set(pg),Y) 

if args.polymake: polymake.polymake_stop(pm)


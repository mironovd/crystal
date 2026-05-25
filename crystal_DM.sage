from numpy import prod
import os
import argparse
import glob
import re
import itertools
from itertools import chain
from collections import Counter
import copy

from sage.symbolic.expression_conversions import laurent_polynomial

from sage_import import sage_import
sys.modules['sage_import'].__dict__['extendpath_rel']('modules')

sage_import('braid_orbit')
sage_import('crystal_op')
sage_import('crystal_step')
sage_import('crystal_generate')
sage_import('braid_orbit_classes_walk')
sage_import('polymake')

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def flatten(nested):
    result = []
    for item in nested:
        if isinstance(item, list):
            result.extend(flatten(item))
        else:
            result.append(item)
    return result



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
required.add_argument('-l', '--minor-level', metavar='dm_level', dest='dm_level', type=int, default=-1,
        help='root of algebra')
optional.add_argument('-i', '--inversions', metavar='inversions', dest='inversions',type=int, default=-1,
        help='process only decompositions with set number of inversions')
optional.add_argument('-n', '--num', metavar='N', dest='num',action='store', type=int, default=0,
        help='process only N random decompositions')
#optional.add_argument('-p', '--polymake',  dest='polymake',action='store_true', default=true,
#        help='use Polymake to check Newton polytopes')
#optional.add_argument('-P', '--no-polymake',  dest='polymake',action='store_false',  default=true,
#        help='do not use Polymake to check Newton polytopes')
optional.add_argument('-v', '--verbose',  dest='verbose',action='store_true', default=false,
        help='verbose info')



parser.add_argument('word',type=int,metavar='word',default=[],nargs='*',
        help='start calculations from word (needs to be reduced word corresponding to longest element of Weyl group)')

args = parser.parse_args()


type=[args.type,args.rank]

pm=None

#if args.polymake:
#    pm=polymake.polymake_start()

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
print(N,w)

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
#if args.type in ['A','D','E']:
#    braid_orbit_method = braid_orbit_classes_walk.braid_orbit_classes_walk_generator

def wminus(wu,ii):
    return crystal_op.minus(wu,wu[ii-1],ii)


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


def grsplit(gr, path, w, i, acc, pacc, tacc):
    print(i,path, wminus(w,i))
#    acc.append((i,path,gr.vertices()))
    if wminus(w,i)>=1:
        gri=gr.copy()
        gri.delete_edges([e for e in gri.edges() if e[2] == wminus(w,i)])
        comps=gri.connected_components_subgraphs()
        print(gri.vertices(), gri.edges())
        if gri.vertices()==[1]:
            return
        cs=[ [comp, set([m.degree(t[i]) for m in comp.vertices() ])] for comp in comps ]
        print([ [ [[m,m.degree(t[i]),i] for m in comp.vertices() ]] for comp in comps ])
        powers = [x[1] for x in cs]
        assert all(len(s) == 1 for s in powers), f"Different powers in components"
        powers = [list(x)[0] for x in powers]
        cs = [ [c[0], list(c[1])[0] ] for c in cs]
        print([c[0].vertices() for c in cs])
        print(powers)
        pacc.append(sorted(powers))
        def si(path):
            j=N
            r=[]
            for q in path:
                if q=='+':
                    r.append('s_'+str(w[j-1]))
                j-=1
            return " ".join(r)

        def addcomp(grr,pw, label):
#            assert len(pw)==1, f"Different powers in component {comp[1]},{comp[0].vertices()}"
#            pw = comp[1]
#            grr=comp[0]
#            print(pw, grr.vertices())
            gr_new=Graph([
                    [x.subs({t[i]:1}) for x in grr.vertices()],
                    [(x[0].subs({t[i]:1}), x[1].subs({t[i]:1}), x[2]) for x in grr.edges()]
                ])
            newpath=path+[label]
            acc.append((i,newpath,si(newpath),grr.vertices()))
            grsplit(gr_new, newpath, w, i-1, acc, pacc, tacc)
            
        


        if len(powers)==1:
            addcomp(cs[0][0], powers[0], '0')
        if len(powers)==2:
#            assert sorted(powers)==[-1,0], "Found two components without [-1,0]"
            if sorted(powers)==[-1,0]:
                grn,p = [item for item in cs if item[1] == -1][0]
                addcomp(grn, p, '-')
                grn, p = [item for item in cs if item[1] == 0][0]
                addcomp(grn, p, '+')
            if min(powers)>=0:
                p=max(powers)
                grn,p = [item for item in cs if item[1] == p][0]
                addcomp(grn, p, '-')
                p=min(powers)
                tacc.append((i, path+['X'], si(path), [item for item in cs if item[1] == p][0][0].vertices()))
# ^^^^^
# 0 1 => 1 ('-')
# 1 2 => 2 ('-')

        if len(powers)>2:
            assert min(powers)==-2 and max(powers)==0, "Found splitting of powers not in form -2 ... 0 "
            grn, p = [item for item in cs if item[1] == -2][0]
            addcomp(grn, p, '-')
            grn, p = [item for item in cs if item[1] == 0][0]
            addcomp(grn, p, '+')
            for x in [item for item in cs if item[1] == -1]:
                tacc.append((i,path+['*'],si(path),x[0].vertices()))
            
            
print(w)
         
for wx in braid_orbit_method(W,w,args.num,args.inversions):
     print("\n==== Case ",qt," ====\n")
     print("Long element decomposition: ",wx,"\n")
     qt+=1
     pt=0
     pg=0
     ptm=[]
     pgm=[]
     vm=vmat(wx,args.rank,t)
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
#         mons.remove(1)

         Gsm = {}
         for mon in mons:
              Gsm[mon]={
                  "in":[],
                  "inmarks":[],
                  "indata":[],
                  "out":[],
                  "outmarks":[],
                  "outdata":[],
              }    
         
         for edge in Gs:
             Gsm[edge[1]]['indata'].append((edge[0][0],edge[0][2]))
             Gsm[edge[1]]['in'].append(edge[0][0])
             Gsm[edge[1]]['inmarks'].append(edge[0][2])
             Gsm[edge[0][0]]['outdata'].append((edge[1],edge[0][2]))
             Gsm[edge[0][0]]['outmarks'].append(edge[0][2])
             Gsm[edge[0][0]]['out'].append(edge[1])
             

#         print(str(Gsm))

#         print(Gs, St)

         i_mons = dict(enumerate([mon for mon in mons if mon!=1]))
         g_mons = [mon for mon in mons if mon!=1]
         mons_i = {mon: i for i, mon in enumerate(mons)}
         mons_ii = [i for i, mon in enumerate(mons)]
#         print(i_mons,mons_i)
         edges = [(mons_i[x[0][0]],mons_i[x[-1]],x[0][-1]) for x in Gs if not (x[0][0]==1 or x[-1]==1)]
         xedges = [(x[0][0],x[-1],x[0][-1]) for x in Gs if not (x[0][0]==1 or x[-1]==1)]
#         print(edges)
         
         G = Graph([g_mons, xedges])

#         print(G.edges(),N)
         gracc=[]
         pacc=[]
         addacc=[]
         grsplit(G,[''],w,N,gracc,pacc,addacc)

         filtered_fr = [item for item in gracc if item[0] == args.dm_level] if args.dm_level!=-1 else gracc

         for item in sorted(filtered_fr, key=lambda t: t[0], reverse=True):
             print(item)

         print("Splitting types: ", Counter(tuple(item) for item in pacc))
         print("Orphaned components found on levels: ", set([x[0] for x in addacc]))

         if args.dm_level==-1:
            print(addacc)
         else:
            print("Orphaned components on level ",args.dm_level,": ",[item for item in addacc if item[0]==args.dm_level])

#     print('Full half-potential: ',pt)
#     if args.polymake: check_poly(ptm,t) 

#if args.polymake: polymake.polymake_stop(pm)





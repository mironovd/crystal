import itertools
from sage.combinat.gray_codes import product as gray_product

def braid_orbit_classes_generator(weyl_group,start_word,num,irrelevant):
    gen = braid_orbit_classes(RootSystem(weyl_group.cartan_type()),start_word)
    n=0
    for el in gen:
        n += 1
        if num>0 and n>num:
            break
        yield el
    return

def braid_orbit_classes(root_system,start_word):
    R=root_system
    L=R.root_lattice()
    X=[x.element for x in R.root_poset()]
    XX=[(x[0],x[0]+x[1],x[1]) for x in itertools.combinations(X,2) if (x[0]+x[1]) in X ]
    start =  root_ordered_set_to_marks(XX,word_to_root_ordered_set(L,start_word))
    mm = start

    graph = make_graph(X,XX,start)    

    n=0

    for m,i in gray_product([2]*len(start)) :
        
        #print('considering candidate', mm)
        if graph.is_directed_acyclic():
            p=graph.topological_sort()
            yield root_ordered_set_to_word(L,p)

        try: graph.reverse_edges([ [XX[m][0],XX[m][1]], [XX[m][0],XX[m][2]] , [XX[m][1],XX[m][2]] ],inplace=True)
        except: graph.reverse_edges([ [XX[m][1],XX[m][0]], [XX[m][2],XX[m][0]] , [XX[m][2],XX[m][1]] ],inplace=True)
        mm[m]=1-mm[m]
    return


def make_graph(vertices,triangles,marks):
    edges = sum([
            [ [triangles[i][0],triangles[i][1]], [triangles[i][0],triangles[i][2]], [triangles[i][1],triangles[i][2]] ] if m==0
            else [ [triangles[i][1],triangles[i][0]], [triangles[i][2],triangles[i][0]], [triangles[i][2],triangles[i][1]]  ] for
            i,m in enumerate(marks)
        ],[])
    return DiGraph([vertices,edges],format='vertices_and_edges')

def root_ordered_set_to_word(root_lattice,root_ordered_set):
    simple_roots = list(root_lattice.simple_roots())
    rank = root_lattice.rank()
    N = len(root_ordered_set)
    w = [None for iiii in range(len(root_ordered_set))]
    w[N-1] = [i+1 for i in range(0,rank) if simple_roots[i]==root_ordered_set[N-1] ][0]
    for i in reversed(range(0,N-1)):
        rt = root_ordered_set[i]
        for s in reversed(w[i+1:N]):
            rt = rt.simple_reflection(s)
        w[i] = [j+1 for j in range(0,rank) if simple_roots[j]==rt ][0]
    return w


def word_to_root_ordered_set(root_lattice,word):
    simple_roots = root_lattice.simple_roots()
    rank = root_lattice.rank()
    N = len(word)
    ro = [simple_roots[word[-1]]]
    for i in reversed(range(0,N-1)):
        rt = simple_roots[word[i]]
        for s in word[i+1:N]:
            rt = rt.simple_reflection(s)
        ro = [rt]+ro
    return ro


def root_ordered_set_to_marks(triangles,root_ordered_set):
    return [ 0 if root_ordered_set.index(k[0])<root_ordered_set.index(k[2]) else 1  for k in triangles]



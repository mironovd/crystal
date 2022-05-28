import itertools
from sage.combinat.gray_codes import product as gray_product

def braid_orbit_classes_walk_generator(weyl_group,start_word,num,irrelevant):
    gen = braid_orbit_classes_walk(RootSystem(weyl_group.cartan_type()),start_word)
    n=0
    for el in gen:
        n += 1
        if num>0 and n>=num:
            break
        yield el
    return

def braid_orbit_classes_walk(root_system,start_word):
    def graph_to_word(L,g):
        if g.is_directed_acyclic():
            p=g.topological_sort()
            return root_ordered_set_to_word(L,p)
        else:
            return None

    R=root_system
    L=R.root_lattice()
    X=[x.element for x in R.root_poset()]
    XX=[(x[0],x[0]+x[1],x[1]) for x in itertools.combinations(X,2) if (x[0]+x[1]) in X ]
    N = len(XX)
    start =  root_ordered_set_to_marks(XX,word_to_root_ordered_set(L,start_word))
    mm = start

    graph = make_graph(X,XX,start)    

    graphs = [ graph ]
    test_mms = [ mm ]
    mms = set([tuple(mm)])

    loop_ind = 0
    list_len = 1

    ans = graph_to_word(L,graph)
    if not ans is None:
#        print("add :",mm,ans,graph)
        yield ans    
    
    n=0

    while loop_ind < list_len:
        test_mm = test_mms[loop_ind]
        graph = graphs[loop_ind]
        loop_ind += 1
        for m in range(0,N):
            new_mm = test_mm.copy()
            new_mm[m] = 1 - new_mm[m]
            new_graph = graph.to_directed()
#            print(new_graph)
            try: new_graph.reverse_edges([ [XX[m][0],XX[m][1]], [XX[m][0],XX[m][2]] , [XX[m][1],XX[m][2]] ],inplace=True)
            except: new_graph.reverse_edges([ [XX[m][1],XX[m][0]], [XX[m][2],XX[m][0]] , [XX[m][2],XX[m][1]] ],inplace=True)
#            print("     ",new_mm,new_graph)
            ans = graph_to_word(L,new_graph)
            if not ans is None:
#                print("ans :",new_mm, ans)
                if tuple(new_mm) not in mms:
#                    print("add :",new_mm,ans,new_graph)
                    test_mms.append(new_mm)
                    graphs.append(new_graph)
                    mms.add(tuple(new_mm))
                    list_len+=1
                    yield ans
            n+=1
    print(n)
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



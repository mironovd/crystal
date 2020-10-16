from itertools import chain

from sage_import import sage_import
sage_import('crystal_op')

class Found(Exception): pass

def Step(I,C,n,wc,mutate_node_func,St,Gs,badroot,t):
    St.append([])
    nodes=St[n]
    all_nodes = list(chain.from_iterable(St))
    res=[]
    new_nodes = [[node,mutate_node_func(I,C,wc,node,badroot,t)] for node in nodes]
    for r in new_nodes:
        node=r[0]
        for new_node in r[1]:
            if (not (new_node[1] in res)) and (not(new_node[1] in all_nodes)):
                res.append(new_node[1])
#                print(new_node)
            Gs.append(((node,new_node[0],new_node[2],new_node[3]),new_node[1]))
    St[n+1]=res
    if len(res)==0 and St[n]!=[1]:
        for old_node in St[n]:
#            try:
                for l in I:
                    for k in crystal_op.RR(wc,node,l,t):# range(1,len(wc)+1):
                        if (node*crystal_op.Ainv(C,wc,l,k,t))==1:
                            res=[1]
                            Gs.append(((old_node,l,k,k),1))
                            find=True
#                            raise Found
#            except Found:
#                1
                        
    St[n+1]=res
    return St, Gs

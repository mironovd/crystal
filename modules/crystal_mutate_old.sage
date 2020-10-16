from sage_import import sage_import
sage_import('crystal_op')


def mutate_node(I,C,w,node,badroot,t):
    poss= [[x[0],x[1]] for x in [ [l,[k for k in crystal_op.RR(w,node,l,t) if not (crystal_op.plus(w,l,k) in crystal_op.RR(w,node,l,t))]] for l in I] if len(x[1])>0]
    poss = [ [[x[0],y] for y in x[1]] for x in poss ]
    poss = [j for i in poss for j in i]
    return [[x[0], node*crystal_op.Ainv(C,w,x[0],x[1],t), x[1],x[1] ]  for x in poss]



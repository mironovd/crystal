from sage_import import sage_import
sage_import('crystal_step')


def generate_crystal(I,C,wc,start,mutate_node_func,badroot,t,verbose):
    xx=true
    st = 0
    St=[[start]]
    Gs=[]
    while xx:
        if verbose:
            print("===\n Step: ",str(st))
        St,Gs = crystal_step.Step(I,C,st,wc,mutate_node_func,St,Gs,badroot,t)
        if verbose:
            print(St[st+1])
        if len(St[st+1])==0 or St[st+1]==[1]:
            xx=false
        st+=1
    return St, Gs


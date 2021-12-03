from sage_import import sage_import
sage_import('crystal_step')


def generate_crystal(I,js,C,wc,start,mutate_node_func,badroot,t,verbose):
    accumulator={}
    sinks=[]
    st = 0
    xx=true
    St=[[start]]
    Gs=[]
    while xx:
        if verbose:
            print("===\n Step: ",str(st))
        St,Gs = crystal_step.Step(I,js,C,st,wc,mutate_node_func,St,Gs,badroot,t,accumulator,sinks)
        if verbose:
            print(St[st+1])
        if len(St[st+1])==0 or St[st+1]==[1]:
            xx=false
        st+=1
    return St, Gs, sinks


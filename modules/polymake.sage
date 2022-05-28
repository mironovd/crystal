import pexpect

def polymake_command(pe, command):
    pe.sendline(command)
    pe.expect ('polytope > ',timeout=300000)
    x=pe.before
    y=x.split(b';\r\n')[1].split(b'\r\n')[0:-1]
    y=[z.decode('utf-8') for z in y]
    return '\n'.join(y)

def check_polytope(pe,polytope):
    polymake_command(pe,'$p=new Polytope(POINTS=>' + str([ [1]+z for z in polytope]) + ');')
    polymake_command(pe,'print $p->VERTICES;')
    vert_n=None
    int_n=None
    bou_n=None
    try:
        vert_n=int(polymake_command(pe,'print $p->N_VERTICES;'))
    except:
        pass
    try:
        int_n=int(polymake_command(pe,'print $p->N_INTERIOR_LATTICE_POINTS;'))
    except:
        pass
    try:
        bou_n=int(polymake_command(pe,'print $p->N_BOUNDARY_LATTICE_POINTS;'))
    except:
        pass
    return (vert_n,int_n,bou_n)




def polymake_start():
    pe = pexpect.spawn('polymake')
    #pe.logfile = sys.stdout
    pe.expect ('polytope > ')
    check_polytope(pe,[[0,0]])
    return(pe)

def polymake_stop(pe):
    pe.sendline('exit;')
    pe.close()


from pymol import stored

def ppos2abegos(ppos):
    abegos=[]
    for resi in range(0,len(ppos[:,0])):
        phi   = float(ppos[resi, 0])
        psi   = float(ppos[resi, 1])
        omega = float(ppos[resi, 2])
        abegos.append(dihd2abego(phi,psi,omega))
    if (len(ppos[:,0] == len(abegos))):
        return abegos
    else:
        return None

def dihd2abego(phi=-60.0,psi=-45.0,omega=180.0,cisbin=30.0,A_upper=50.0,A_lower=-75.0,G_upper=100,G_lower=-100.0):
    if (omega <= cisbin and omega >= -cisbin):
        return "O"

    if (phi <= 0):
        if (A_lower <= psi and psi <= A_upper):
            return "A"
        else:
            return "B"
    else:
        if (G_lower <= psi and psi <= G_upper):
            return "G"
        else:
            return "E"

import numpy as np

def getomega(target):
    stored.resi = []
    cmd.select("tmpsel", target)
    cmd.iterate("tmpsel and name CA", "stored.resi.append(int(resi))")
    resis = np.array(stored.resi)
    #print(resis)
    omegas=[]
    for resi in resis[1:]:
        #print(resi)    
        om=cmd.get_dihedral("tmpsel and name C and resi " + str(resi-1),
                         "tmpsel and name O and resi " + str(resi-1),
                         "tmpsel and name N and resi " + str(resi),
                         "tmpsel and name CA and resi " + str(resi))
        omegas.append(om)
    omegas.append(180)
#    omegas.append(180)
    return omegas

def labegO(target="all"):
    objs=cmd.get_object_list(target)
    #print(objs)
    for tgt in objs:
        tgt = tgt + " and " + target
        phipsi=cmd.phi_psi(tgt)
        keys=list(phipsi)
        ppos=[]
        omegas=getomega(tgt)
        count=0
        for key in keys:                        
            #print(omega)
            phi = phipsi[key][0]
            psi = phipsi[key][1]       
            omega = omegas[count]
            count = count + 1
            ppos.append([phi,psi,omega])
        
        ppos=np.array(ppos)
        myabego = ppos2abegos(ppos)
        ii = 0
        for key in keys:
            label=str(myabego[ii])
            # Be careful of label should be double quoted when evaluated
            qlabel='\"' + label + '\"'
            cmd.label(tgt +" and index "+str(key[1]),expression=qlabel)
            ii+=1
            cmd.set("label_size","30")

pymol.cmd.extend("labegO", labegO)
cmd.auto_arg[0]['labegO'] = cmd.auto_arg[0]['align']


def labego(target="all", option=""):
    objs=cmd.get_object_list(target)
    for tgt in objs:
        tgt = tgt + " and " + target
        phipsi=cmd.phi_psi(tgt)
        keys=list(phipsi)
        ppos=[]

        for key in keys:
            phi = phipsi[key][0]
            psi = phipsi[key][1]
            omega = 180.0 # temporary
            ppos.append([phi,psi,omega])
        ppos=np.array(ppos)
        myabego = ppos2abegos(ppos)
        ii = 0
        for key in keys:
            label=str(myabego[ii])
            # Be careful of label should be double quoted when evaluated
            qlabel='\"' + label + '\"'
            cmd.label(tgt +" and index "+str(key[1]),expression=qlabel)
            ii+=1
            cmd.set("label_size","30")

pymol.cmd.extend("labego", labego)
cmd.auto_arg[0]['labego'] = cmd.auto_arg[0]['align']

from pymol import stored

def ppos2abegos(ppos, divideB):
    abegos = []    
    if ppos.ndim == 2:
        for resi in range(0, len(ppos[:, 0])):
            phi = float(ppos[resi, 0])
            psi = float(ppos[resi, 1])
            omega = float(ppos[resi, 2])
            abegos.append(dihd2abego(phi, psi, omega, divideB=divideB))
        if (len(ppos[:, 0] == len(abegos))):
            return abegos
        else:
            return None
    else:
        return None

def dihd2abego(phi=-60.0,psi=-45.0,omega=180.0,
               cisbin=30.0,A_upper=50.0,A_lower=-75.0,
               G_upper=100,G_lower=-100.0,
                P_thre = -90, divideB=False):
    if (omega <= cisbin and omega >= -cisbin):
        return "O"

    if (phi <= 0):
        if (A_lower <= psi and psi <= A_upper):
            return "A"
        else:
            if (divideB):
                if ( phi > P_thre):
                    return "P"
                else:
                    return "S"
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

def labegO(target="polymer.protein"):
    objs=cmd.get_object_list(target)
    print(objs)
    for tgt_tmp in objs:
        tgt = tgt_tmp + " and " + target
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
        myabego = ppos2abegos(ppos,False)
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

def lapsegO(target="polymer.protein"):
    objs=cmd.get_object_list(target)
    print(objs)
    for tgt_tmp in objs:
        tgt = tgt_tmp + " and " + target
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
        myabego = ppos2abegos(ppos,True)
        ii = 0
        for key in keys:
            label=str(myabego[ii])
            # Be careful of label should be double quoted when evaluated
            qlabel='\"' + label + '\"'
            cmd.label(tgt +" and index "+str(key[1]),expression=qlabel)
            ii+=1
            cmd.set("label_size","30")

pymol.cmd.extend("lapsegO", lapsegO)
cmd.auto_arg[0]['lapsegO'] = cmd.auto_arg[0]['delete']


def labego(target="polymer.protein"):
    objs=cmd.get_object_list(target)
    for tgt_tmp in objs:
        tgt = tgt_tmp + " and " + target
        phipsi=cmd.phi_psi(tgt)
        keys=list(phipsi)
        ppos=[]

        for key in keys:
            phi = phipsi[key][0]
            psi = phipsi[key][1]
            omega = 180.0 # temporary
            ppos.append([phi,psi,omega])
        ppos=np.array(ppos)
        myabego = ppos2abegos(ppos,divideB=False)
        ii = 0
        for key in keys:
            label=str(myabego[ii])
            # Be careful of label should be double quoted when evaluated
            qlabel='\"' + label + '\"'
            cmd.label(tgt +" and index "+str(key[1]),expression=qlabel)
            ii+=1
            cmd.set("label_size","30")

def ramaval(target="polymer.protein"):
    objs=cmd.get_object_list(target)
    for tgt_tmp in objs:
        tgt = tgt_tmp + " and " + target
        phipsi=cmd.phi_psi(tgt)
        keys=list(phipsi)
        ppos=[]

        for key in keys:
            phi = phipsi[key][0]
            psi = phipsi[key][1]
            omega = 180.0 # temporary
            ppos.append([phi,psi,omega])
        ppos=np.array(ppos)
        myabego = ppos2abegos(ppos,divideB=False)
        ii = 0
        for key in keys:
            label=str(round(ppos[ii,0]))+"_"+str(int(ppos[ii,1]))
            # Be careful of label should be double quoted when evaluated
            qlabel='\"' + label + '\"'
            cmd.label(tgt +" and index "+str(key[1]),expression=qlabel)
            ii+=1
            cmd.set("label_size","30")

def ramabego(target="polymer.protein"):
    objs=cmd.get_object_list(target)
    for tgt_tmp in objs:
        tgt = tgt_tmp + " and " + target
        phipsi=cmd.phi_psi(tgt)
        keys=list(phipsi)
        ppos=[]

        for key in keys:
            phi = phipsi[key][0]
            psi = phipsi[key][1]
            omega = 180.0 # temporary
            ppos.append([phi,psi,omega])
        ppos=np.array(ppos)
        myabego = ppos2abegos(ppos,divideB=False)
        ii = 0
        for key in keys:
            label=str(myabego[ii])+"_"+str(round(ppos[ii,0]))+"_"+str(int(ppos[ii,1]))
            # Be careful of label should be double quoted when evaluated
            qlabel='\"' + label + '\"'
            cmd.label(tgt +" and index "+str(key[1]),expression=qlabel)
            ii+=1
            cmd.set("label_size","30")
            cmd.set("label_position",[0,0,10])    
            

def lapsego(target="polymer.protein"):
    objs=cmd.get_object_list(target)
    for tgt_tmp in objs:
        tgt = tgt_tmp + " and " + target
        phipsi=cmd.phi_psi(tgt)
        keys=list(phipsi)
        ppos=[]

        for key in keys:
            phi = phipsi[key][0]
            psi = phipsi[key][1]
            omega = 180.0 # temporary
            ppos.append([phi,psi,omega])
        ppos=np.array(ppos)
        myabego = ppos2abegos(ppos,divideB=True)
        ii = 0
        for key in keys:
            label=str(myabego[ii])
            # Be careful of label should be double quoted when evaluated
            qlabel='\"' + label + '\"'
            cmd.label(tgt +" and index "+str(key[1]),expression=qlabel)
            ii+=1
            cmd.set("label_size","30")
            cmd.set("label_position",[0,0,10])    

def ramapsego(target="polymer.protein"):
    objs=cmd.get_object_list(target)
    for tgt_tmp in objs:
        tgt = tgt_tmp + " and " + target
        phipsi=cmd.phi_psi(tgt)
        keys=list(phipsi)
        ppos=[]

        for key in keys:
            phi = phipsi[key][0]
            psi = phipsi[key][1]
            omega = 180.0 # temporary
            ppos.append([phi,psi,omega])
        ppos=np.array(ppos)
        myabego = ppos2abegos(ppos,divideB=True)
        ii = 0
        for key in keys:
            label=str(myabego[ii])+"_"+str(round(ppos[ii,0]))+"_"+str(int(ppos[ii,1]))
            # Be careful of label should be double quoted when evaluated
            qlabel='\"' + label + '\"'
            cmd.label(tgt +" and index "+str(key[1]),expression=qlabel)
            ii+=1
            cmd.set("label_size","30")            

pymol.cmd.extend("labego", labego)
pymol.cmd.extend("ramaval", ramaval)
pymol.cmd.extend("ramabego", ramabego)
pymol.cmd.extend("ramapsego", ramapsego)
pymol.cmd.extend("lapsego", lapsego)

cmd.auto_arg[0]['ramaval'] = cmd.auto_arg[0]['delete']
cmd.auto_arg[0]['labego'] = cmd.auto_arg[0]['delete']
cmd.auto_arg[0]['lapsego'] = cmd.auto_arg[0]['delete']
cmd.auto_arg[0]['ramapsego'] = cmd.auto_arg[0]['delete']
cmd.auto_arg[0]['ramabego'] = cmd.auto_arg[0]['delete']

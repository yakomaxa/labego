import pymol
import re
from pymol import stored
import numpy as np


def ppos2abegos(ppos, divideB):
    abegos = []
    for resi in range(0, len(ppos[:, 0])):
        phi = float(ppos[resi, 0])
        psi = float(ppos[resi, 1])
        omega = float(ppos[resi, 2])
        abegos.append(dihd2abego(phi, psi, omega, divideB=divideB))
    if (len(ppos[:, 0] == len(abegos))):
        return abegos
    else:
        return None


def dihd2abego(phi=-60.0, psi=-45.0, omega=180.0,
               cisbin=30.0, A_upper=50.0, A_lower=-75.0,
               G_upper=100, G_lower=-100.0,
               P_thre=-90, divideB=False):
    if (omega <= cisbin and omega >= -cisbin):
        return "O"

    if (phi <= 0):
        if (A_lower <= psi and psi <= A_upper):
            return "A"
        else:
            if (divideB):
                if (phi > P_thre):
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


def getomega(target):
    stored.resi = []
    cmd.select("tmpsel", target)
    cmd.iterate("tmpsel and name CA", "stored.resi.append(int(resi))")
    resis = np.array(stored.resi)
    # print(resis)
    omegas = []
    for resi in resis[1:]:
        # print(resi)
        om = cmd.get_dihedral("tmpsel and name C and resi " + str(resi - 1),
                              "tmpsel and name O and resi " + str(resi - 1),
                              "tmpsel and name N and resi " + str(resi),
                              "tmpsel and name CA and resi " + str(resi))
        omegas.append(om)
    omegas.append(180)
    #    omegas.append(180)
    return omegas


def labegO(target="polymer.protein"):
    objs = cmd.get_object_list(target)
    print(objs)
    for tgt_tmp in objs:
        tgt = tgt_tmp + " and " + target
        phipsi = cmd.phi_psi(tgt)
        keys = list(phipsi)
        ppos = []
        omegas = getomega(tgt)
        count = 0
        for key in keys:
            # print(omega)
            phi = phipsi[key][0]
            psi = phipsi[key][1]
            omega = omegas[count]
            count = count + 1
            ppos.append([phi, psi, omega])

        ppos = np.array(ppos)
        myabego = ppos2abegos(ppos, False)
        ii = 0
        for key in keys:
            label = str(myabego[ii])
            # Be careful of label should be double quoted when evaluated
            qlabel = '\"' + label + '\"'
            cmd.label(tgt + " and index " + str(key[1]), expression=qlabel)
            ii += 1
            cmd.set("label_size", "30")


pymol.cmd.extend("labegO", labegO)
cmd.auto_arg[0]['labegO'] = cmd.auto_arg[0]['align']


def lapsegO(target="polymer.protein"):
    objs = cmd.get_object_list(target)
    print(objs)
    for tgt_tmp in objs:
        tgt = tgt_tmp + " and " + target
        phipsi = cmd.phi_psi(tgt)
        keys = list(phipsi)
        ppos = []
        omegas = getomega(tgt)
        count = 0
        for key in keys:
            # print(omega)
            phi = phipsi[key][0]
            psi = phipsi[key][1]
            omega = omegas[count]
            count = count + 1
            ppos.append([phi, psi, omega])

        ppos = np.array(ppos)
        myabego = ppos2abegos(ppos, True)
        ii = 0
        for key in keys:
            label = str(myabego[ii])
            # Be careful of label should be double quoted when evaluated
            qlabel = '\"' + label + '\"'
            cmd.label(tgt + " and index " + str(key[1]), expression=qlabel)
            ii += 1
            cmd.set("label_size", "30")


pymol.cmd.extend("lapsegO", lapsegO)
cmd.auto_arg[0]['lapsegO'] = cmd.auto_arg[0]['delete']


def labego(target="polymer.protein"):
    objs = cmd.get_object_list(target)
    for tgt_tmp in objs:
        tgt = tgt_tmp + " and " + target
        phipsi = cmd.phi_psi(tgt)
        keys = list(phipsi)
        ppos = []

        for key in keys:
            phi = phipsi[key][0]
            psi = phipsi[key][1]
            omega = 180.0  # temporary
            ppos.append([phi, psi, omega])
        ppos = np.array(ppos)
        myabego = ppos2abegos(ppos, divideB=False)
        ii = 0
        for key in keys:
            label = str(myabego[ii])
            # Be careful of label should be double quoted when evaluated
            qlabel = '\"' + label + '\"'
            cmd.label(tgt + " and index " + str(key[1]), expression=qlabel)
            ii += 1
            cmd.set("label_size", "30")


def ramaval(target="polymer.protein"):
    objs = cmd.get_object_list(target)
    for tgt_tmp in objs:
        tgt = tgt_tmp + " and " + target
        phipsi = cmd.phi_psi(tgt)
        keys = list(phipsi)
        ppos = []

        for key in keys:
            phi = phipsi[key][0]
            psi = phipsi[key][1]
            omega = 180.0  # temporary
            ppos.append([phi, psi, omega])
        ppos = np.array(ppos)
        myabego = ppos2abegos(ppos, divideB=False)
        ii = 0
        for key in keys:
            label = str(round(ppos[ii, 0])) + "_" + str(int(ppos[ii, 1]))
            # Be careful of label should be double quoted when evaluated
            qlabel = '\"' + label + '\"'
            cmd.label(tgt + " and index " + str(key[1]), expression=qlabel)
            ii += 1
            cmd.set("label_size", "30")


def ramabego(target="polymer.protein"):
    objs = cmd.get_object_list(target)
    for tgt_tmp in objs:
        tgt = tgt_tmp + " and " + target
        phipsi = cmd.phi_psi(tgt)
        keys = list(phipsi)
        ppos = []

        for key in keys:
            phi = phipsi[key][0]
            psi = phipsi[key][1]
            omega = 180.0  # temporary
            ppos.append([phi, psi, omega])
        ppos = np.array(ppos)
        myabego = ppos2abegos(ppos, divideB=False)
        ii = 0
        for key in keys:
            label = str(myabego[ii]) + "_" + str(round(ppos[ii, 0])) + "_" + str(int(ppos[ii, 1]))
            # Be careful of label should be double quoted when evaluated
            qlabel = '\"' + label + '\"'
            cmd.label(tgt + " and index " + str(key[1]), expression=qlabel)
            ii += 1
            cmd.set("label_size", "30")
            cmd.set("label_position", [0, 0, 10])


def lapsego(target="polymer.protein"):
    objs = cmd.get_object_list(target)
    for tgt_tmp in objs:
        tgt = tgt_tmp + " and " + target
        phipsi = cmd.phi_psi(tgt)
        keys = list(phipsi)
        ppos = []

        for key in keys:
            phi = phipsi[key][0]
            psi = phipsi[key][1]
            omega = 180.0  # temporary
            ppos.append([phi, psi, omega])
        ppos = np.array(ppos)
        myabego = ppos2abegos(ppos, divideB=True)
        ii = 0
        for key in keys:
            label = str(myabego[ii])
            # Be careful of label should be double quoted when evaluated
            qlabel = '\"' + label + '\"'
            cmd.label(tgt + " and index " + str(key[1]), expression=qlabel)
            ii += 1
            cmd.set("label_size", "30")
            cmd.set("label_position", [0, 0, 10])


def ramapsego(target="polymer.protein"):
    objs = cmd.get_object_list(target)
    for tgt_tmp in objs:
        tgt = tgt_tmp + " and " + target
        phipsi = cmd.phi_psi(tgt)
        keys = list(phipsi)
        ppos = []

        for key in keys:
            phi = phipsi[key][0]
            psi = phipsi[key][1]
            omega = 180.0  # temporary
            ppos.append([phi, psi, omega])
        ppos = np.array(ppos)
        myabego = ppos2abegos(ppos, divideB=True)
        ii = 0
        for key in keys:
            label = str(myabego[ii]) + "_" + str(round(ppos[ii, 0])) + "_" + str(int(ppos[ii, 1]))
            # Be careful of label should be double quoted when evaluated
            qlabel = '\"' + label + '\"'
            cmd.label(tgt + " and index " + str(key[1]), expression=qlabel)
            ii += 1
            cmd.set("label_size", "30")


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


####################################
def get_abego_gap_filled(target="all", divideB=False):
    aa1 = list("ACDEFGHIKLMNPQRSTVWY")
    aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
    aa3to1 = dict(zip(aa3, aa1))
    _target = target
    chains = pymol.cmd.get_chains(target)

    fastas = []
    abegos = []
    bfactors = []
    min_indexes = []

    for chain in chains:
        target = _target + " and chain " + chain

        keep = "A"
        objct = target
        remStr = "(" + "%s and not (alt ''+%s)" % (objct, keep) + ")"

        selection = target + " and polymer.protein" + " and not " + remStr

        stored.bfactor = []
        stored.ca_index = []
        stored.ca_resname = []
        stored.ca_resi = []

        pymol.cmd.iterate(selection + " and name ca", "stored.bfactor.append(b)")
        pymol.cmd.iterate(selection + " and name ca", "stored.ca_index.append(index)")
        pymol.cmd.iterate(selection + " and name ca", "stored.ca_resname.append(resn)")
        pymol.cmd.iterate(selection + " and name ca", "stored.ca_resi.append(resi)")

        if 0 in stored.bfactor:
            fasta_visible = []
            resi_visible = []
            index_visible = []
            b_visible = []
            print("zero b-factor!")
            index = -1
            for b in stored.bfactor:
                index = index + 1
                if b > 0.0:
                    fasta_visible.append(stored.ca_resname[index])
                    resi_visible.append(stored.ca_resi[index])
                    index_visible.append(stored.ca_index[index])
                    b_visible.append(b)

            ca_resname = fasta_visible
            ca_index = index_visible
            ca_resi = resi_visible
            ca_bfactor = b_visible
        else:
            ca_resname = stored.ca_resname
            ca_index = stored.ca_index
            ca_resi = stored.ca_resi
            ca_bfactor = stored.bfactor

        aa_1 = []
        for aa_3 in ca_resname:
            aa_1.append(aa3to1[aa_3])
        myfasta = aa_1

        phipsi = pymol.cmd.phi_psi(selection)
        keys = list(phipsi)

        ca_index = np.array(ca_index)
        ca_resi = np.array(ca_resi)

        ppos = []
        resis = []
        omega = 180.0  # temporary
        for key in keys:
            resis.append(int(ca_resi[(key[1] == ca_index)][0]))
            phi = phipsi[key][0]
            psi = phipsi[key][1]
            ppos.append([phi, psi, omega])

        ppos = np.array(ppos)
        myabego = ppos2abegos(ppos, divideB=divideB)

        resis = np.array(resis)
        gap = resis[1:] - resis[0:-1] - 1
        index = np.arange(len(gap))
        index_missing = index[gap > 0]
        gap_sizes = gap[gap > 0]
        for site, length in zip(list(reversed(index_missing)), list(reversed(gap_sizes))):
            ins = "X" * (length)
            myabego.insert(site + 1, ins)

        intersect = ca_resi.astype(int)
        min_index = min(intersect)
        followed_by_gaps = intersect[1:] - intersect[0:-1] > 1
        gap_sizes = (intersect[1:] - intersect[0:-1] - 1)[followed_by_gaps]
        index = np.arange(stop=len(followed_by_gaps))
        index = index[followed_by_gaps]
        for site, length in zip(list(reversed(index)), list(reversed(gap_sizes))):
            ins = "B" * (length)
            myfasta.insert(site + 1, ins)
            for l in range(length, 0, -1):
                ca_bfactor.insert(site + l, 999)

        myabego = list("".join(myabego))
        myfasta = list("X" + "".join(myfasta) + "X")

        myfasta = "".join(myfasta)
        myabego = "X" + "".join(myabego) + "X"

        fastas.append(myfasta)
        abegos.append(myabego)
        bfactors.append(ca_bfactor)
        min_indexes.append(min_index)

    return abegos, fastas, bfactors, min_indexes, chains


def iterative_search(sequence=None, query=None):
    result = re.search(query, sequence)
    if (result == None):
        return None, None
    else:
        cumm = 0
        inits = []
        ends = []
        sequence_tmp = sequence
        while (len(sequence_tmp) > 0):
            result = re.search(query, sequence_tmp)
            if (result == None):
                break
            else:
                i = result.span()[0]
                e = result.span()[1]
                abs_i = i + cumm
                abs_e = e + cumm - 1

                inits.append(abs_i)
                ends.append(abs_e)

                sequence_tmp = sequence_tmp[e:]
                cumm += e
    return inits, ends


def _abego(target, query="GBB", mode="labego"):
    for obj in pymol.cmd.get_object_list(target):
        abegos, _, _, min_resis, chains = get_abego_gap_filled(target=obj)
        for abego, min_resi, chain in zip(abegos, min_resis, chains):
            inits, ends = iterative_search(sequence=abego, query=query)
            if inits is None:
                print(query + " is not found in chain " + chain)
            else:
                diff = min_resi - 1
                for ii, ee in zip(inits, ends):
                    i = diff + ii + 1
                    e = diff + ee + 1
                    if mode == "labego":
                        labego(obj + " and resi " + str(i) + "-" + str(e) + " and chain " + chain)
                    elif mode == "lapsego":
                        lapsego(obj + " and resi " + str(i) + "-" + str(e) + " and chain " + chain)
                    elif mode == "labegO":
                        labegO(obj + " and resi " + str(i) + "-" + str(e) + " and chain " + chain)
                    elif mode == "lapsegO":
                        lapsegO(obj + " and resi " + str(i) + "-" + str(e) + " and chain " + chain)
                    elif mode == "ramabego":
                        ramabego(obj + " and resi " + str(i) + "-" + str(e) + " and chain " + chain)
                    elif mode == "ramapsego":
                        ramapsego(obj + " and resi " + str(i) + "-" + str(e) + " and chain " + chain)
    return 0


def abego_labego(target=None, query="GBB"):
    if target is None:
        target = "all"
    _abego(target=target, query=query, mode="labego")
    return 0


def abego_lapsego(target=None, query="GBB"):
    if target is None:
        target = "all"
    _abego(target=target, query=query, mode="lapsego")
    return 0


def abego_labegO(target=None, query="GBB"):
    if target is None:
        target = "all"
    _abego(target=target, query=query, mode="labegO")
    return 0


def abego_lapsegO(target=None, query="GBB"):
    if target is None:
        target = "all"
    _abego(target=target, query=query, mode="lapsegO")
    return 0


def abego_ramabego(target=None, query="GBB"):
    if target is None:
        target = "all"
    _abego(target=target, query=query, mode="ramabego")
    return 0


def abego_ramapsego(target=None, query="GBB"):
    if target is None:
        target = "all"
    _abego(target=target, query=query, mode="ramapsego")
    return 0


pymol.cmd.extend("abego_labego", abego_labego)
pymol.cmd.extend("abego_labegO", abego_labegO)
pymol.cmd.extend("abego_lapsego", abego_lapsego)
pymol.cmd.extend("abego_lapsegO", abego_lapsegO)
pymol.cmd.extend("abego_ramabego", abego_ramabego)
pymol.cmd.extend("abego_ramapsego", abego_ramapsego)

cmd.auto_arg[0]['abego_labego'] = cmd.auto_arg[0]['delete']
cmd.auto_arg[0]['abego_labegO'] = cmd.auto_arg[0]['delete']
cmd.auto_arg[0]['abego_lapsego'] = cmd.auto_arg[0]['delete']
cmd.auto_arg[0]['abego_lapsegO'] = cmd.auto_arg[0]['delete']
cmd.auto_arg[0]['abego_ramabego'] = cmd.auto_arg[0]['delete']
cmd.auto_arg[0]['abego_ramapsego'] = cmd.auto_arg[0]['delete']


def abego_show(visual="line", target=None, query="GBB", divideB=False):
    if target is None:
        target = "all"

    for obj in pymol.cmd.get_object_list(target):
        abegos, _, _, min_resis, chains = get_abego_gap_filled(target=obj, divideB=divideB)
        for abego, min_resi, chain in zip(abegos, min_resis, chains):
            inits, ends = iterative_search(sequence=abego, query=query)
            if inits is None:
                print(query + " is not found in chain " + chain)
            else:
                diff = min_resi - 1
                for ii, ee in zip(inits, ends):
                    i = diff + ii + 1
                    e = diff + ee + 1
                    pymol.cmd.show(visual, obj + " and resi " + str(i) + "-" + str(e) + " and chain " + chain)
    return 0


def abego_hide(visual="everything", target=None, query="GBB", divideB=False):
    if target is None:
        target = "all"

    for obj in pymol.cmd.get_object_list(target):
        abegos, _, _, min_resis, chains = get_abego_gap_filled(target=obj, divideB=divideB)
        for abego, min_resi, chain in zip(abegos, min_resis, chains):
            inits, ends = iterative_search(sequence=abego, query=query)
            if inits is None:
                print(query + " is not found in chain " + chain)
            else:
                diff = min_resi - 1
                for ii, ee in zip(inits, ends):
                    i = diff + ii + 1
                    e = diff + ee + 1
                    pymol.cmd.hide(visual, obj + " and resi " + str(i) + "-" + str(e) + " and chain " + chain)
    return 0


def abego_color(color="white", target=None, query="GBB", divideB=False):
    if target is None:
        target = "all"

    for obj in pymol.cmd.get_object_list(target):
        abegos, _, _, min_resis, chains = get_abego_gap_filled(target=obj, divideB=divideB)
        for abego, min_resi, chain in zip(abegos, min_resis, chains):
            inits, ends = iterative_search(sequence=abego, query=query)
            if inits is None:
                print(query + " is not found in chain " + chain)
            else:
                diff = min_resi - 1
                for ii, ee in zip(inits, ends):
                    i = diff + ii + 1
                    e = diff + ee + 1
                    pymol.cmd.color(color, obj + " and resi " + str(i) + "-" + str(e) + " and chain " + chain)
    return 0


def abego_select(name=None, target=None, query="GBB", divideB=False):
    if name is None:
        name = query

    if target is None:
        target = "all"

    for obj in pymol.cmd.get_object_list(target):
        abegos, _, _, min_resis, chains = get_abego_gap_filled(target=obj, divideB=divideB)
        for abego, min_resi, chain in zip(abegos, min_resis, chains):
            inits, ends = iterative_search(sequence=abego, query=query)
            if inits is None:
                print(query + " is not found in chain " + chain)
            else:
                diff = min_resi - 1
                index = 0
                for ii, ee in zip(inits, ends):
                    index += 1
                    i = diff + ii + 1
                    e = diff + ee + 1
                    pymol.cmd.select(obj + "_" + name + "_" + str(index).zfill(4),
                                     target + " and resi " + str(i) + "-" + str(e) + " and chain " + chain)
    return 0


def abego_fit(mobile=None, target=None, query="GBB", index_t=0, index_m=0, mode="pair_fit", divideB=False):
    abegos_t, _, _, min_resis_t, chains_t = get_abego_gap_filled(target=target, divideB=divideB)
    abegos_m, _, _, min_resis_m, chains_m = get_abego_gap_filled(target=mobile, divideB=divideB)

    i_ts = []
    e_ts = []
    c_ts = []
    for abego_t, min_resi_t, chain_t in zip(abegos_t, min_resis_t, chains_t):
        inits_t, ends_t = iterative_search(sequence=abego_t, query=query)
        diff_t = min_resi_t - 1
        if inits_t is not None:
            for i, e in zip(inits_t, ends_t):
                i_ts.append(diff_t + i)
                e_ts.append(diff_t + e)
                c_ts.append(chain_t)

    i_ms = []
    e_ms = []
    c_ms = []
    for abego_m, min_resi_m, chain_m in zip(abegos_m, min_resis_m, chains_m):
        inits_m, ends_m = iterative_search(sequence=abego_m, query=query)
        diff_m = min_resi_m - 1
        if inits_m is not None:
            for i, e in zip(inits_m, ends_m):
                i_ms.append(diff_m + i)
                e_ms.append(diff_m + e)
                c_ms.append(chain_m)

    if (len(i_ms) == 0):
        print("target does not have " + query)
        return 0

    if (len(i_ms) == 0):
        print("mobile does not have " + query)
        return 0

    if len(i_ms) > int(index_m) and len(i_ts) > int(index_t):
        i_m = i_ms[int(index_m)]
        e_m = e_ms[int(index_m)]
        c_m = c_ms[int(index_m)]

        i_t = i_ts[int(index_t)]
        e_t = e_ts[int(index_t)]
        c_t = c_ts[int(index_t)]

        pymol.cmd.do(mode + " " + mobile + " and resi " + str(i_m) + "-" + str(e_m) + " and name ca and chain " + c_m + " , "
                     + target + " and resi " + str(i_t) + "-" + str(e_t) + " and name ca and chain " + c_t)
    else:
        print("You don't have enough motifs")

    return 0

def abego_fitto(target=None, query="GBB", index_t=0, index_m=0, mode="pair_fit", divideB=False):
    for obj in pymol.cmd.get_object_list("all"):
        abego_fit(obj,target,query,index_t,index_m)
    return

pymol.cmd.extend("abego_fitto", abego_fitto)
pymol.cmd.auto_arg[0]['abego_fitto'] = pymol.cmd.auto_arg[0]['align']

def abego_create(name=None, target=None, query="AAAGBBAAA", divideB=False):

    if name is None:
        name = query

    if target is None:
        target = "all"

    for obj in pymol.cmd.get_object_list(target):
        abego, _, _, min_resi, chain = get_abego_gap_filled(target=obj, divideB=divideB)
        inits, ends = iterative_search(sequence=abego, query=query)
        diff = min_resi - 1
        index = 0
        for ii, ee in zip(inits, ends):
            index += 1
            i = diff + ii + 1
            e = diff + ee + 1
            pymol.cmd.create(target + "_" + name + "_" + str(index).zfill(4),
                             target + " and resi " + str(i) + "-" + str(e) + " and chain " + chain)
        return 0


pymol.cmd.extend("abego_show", abego_show)
pymol.cmd.auto_arg[0]['abego_show'] = pymol.cmd.auto_arg[0]['show']
pymol.cmd.auto_arg[1]['abego_show'] = pymol.cmd.auto_arg[1]['show']

pymol.cmd.extend("abego_hide", abego_hide)
pymol.cmd.auto_arg[0]['abego_hide'] = pymol.cmd.auto_arg[0]['hide']
pymol.cmd.auto_arg[1]['abego_hide'] = pymol.cmd.auto_arg[1]['hide']

pymol.cmd.extend("abego_color", abego_color)
pymol.cmd.auto_arg[0]['abego_color'] = pymol.cmd.auto_arg[0]['color']
pymol.cmd.auto_arg[1]['abego_color'] = pymol.cmd.auto_arg[1]['color']

pymol.cmd.extend("abego_select", abego_select)
# pymol.cmd.auto_arg[0]['abego_select'] = pymol.cmd.auto_arg[0]['select']
pymol.cmd.auto_arg[1]['abego_select'] = pymol.cmd.auto_arg[1]['select']

pymol.cmd.extend("abego_create", abego_create)
pymol.cmd.auto_arg[0]['abego_create'] = pymol.cmd.auto_arg[0]['create']
pymol.cmd.auto_arg[1]['abego_create'] = pymol.cmd.auto_arg[1]['create']

pymol.cmd.extend("abego_fit", abego_fit)
pymol.cmd.auto_arg[0]['abego_fit'] = pymol.cmd.auto_arg[0]['align']
pymol.cmd.auto_arg[1]['abego_fit'] = pymol.cmd.auto_arg[1]['align']
# pymol.cmd.auto_arg[4]['abego_fit'] = [lambda: pymol.cmd.Shortcut(['super','align','mican']), '4th argument', ', ']

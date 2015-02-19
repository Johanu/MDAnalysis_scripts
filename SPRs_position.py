from __future__ import division
import matplotlib.pyplot as plt
# import numpy as np
import MDAnalysis as md
# from MDAnalysis import *
# from MDAnalysis.analysis.distances import dist


def calculate_dists(gro_file, xtc_file, ligand_name):
    u = md.Universe(gro_file, xtc_file)
    select_ligand = u.selectAtoms("resname " + ligand_name)
    select_res31 = u.selectAtoms("backbone and resnum 31")
    COM_distance = []
    print gro_file[:-4], select_ligand

    for i in u.trajectory:
        ligand_COM = select_ligand.centerOfMass()
        res31_COM = select_res31.centerOfMass()
        #COM_distance.append(res31_COM[2] - ligand_COM[2])
        COM_distance.append(ligand_COM[2] - res31_COM[2])
    print max(COM_distance), COM_distance.index(max(COM_distance))
    print min(COM_distance), COM_distance.index(min(COM_distance))
    return COM_distance

S31N = calculate_dists('./mutants/2HSP/2HSP_S31N/ligands/2HSP_S31N_SPRP/minim_new.gro', './mutants/2HSP/2HSP_S31N/ligands/2HSP_S31N_SPRP/equil_new.xtc', 'SPR')
V27A = calculate_dists('./mutants/2HSP/2HSP_V27A/ligands/2HSP_V27A_SPRP/minim_new.gro', './mutants/2HSP/2HSP_V27A/ligands/2HSP_V27A_SPRP/equil_new.xtc', 'SPR')
SPRP = calculate_dists('./WT_2HSP/ligands/WT_2HSP_SPRP/minim_new.gro', './WT_2HSP/ligands/WT_2HSP_SPRP/equil_new.xtc', 'SPR')
SPRC = calculate_dists('./WT_2HSP/ligands/WT_2HSP_SPRC/minim.gro', './WT_2HSP/ligands/WT_2HSP_SPRC/equil.xtc', 'SPRC')
SPRN = calculate_dists('./WT_2HSP/ligands/WT_2HSP_SPRN/minim.gro', './WT_2HSP/ligands/WT_2HSP_SPRN/equil.xtc', 'SPRN')
x_vals = [x / 100 for x in range(0, len(SPRP))]
x_S31N = [x / 100 for x in range(0, len(S31N))]
plt.plot(x_S31N, S31N, label='SPRP S31N')
plt.plot(x_vals, V27A, label='SPRP V27A')
plt.plot(x_vals, SPRP, label='SPRP WT')
plt.plot(x_vals, SPRC, label='SPRC WT')
plt.plot(x_vals, SPRN, label='SPRN WT')
leg = plt.legend(ncol=3, loc=9, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.xlabel('Time / ns')
plt.ylabel(ur'Z difference to resnum 31 COM / $\AA$')
# plt.show()
plt.savefig('SPRs_COM.png', dpi=300)
plt.close()

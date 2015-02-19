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
V27A_1 = calculate_dists('./mutants/2HSP/2HSP_V27A/ligands/2HSP_V27A_SPRP/minim_new.gro', './mutants/2HSP/2HSP_V27A/ligands/2HSP_V27A_SPRP/equil_new.xtc', 'SPR')
V27A_2 = calculate_dists('/nike-storage/passagrilli/MD_inhibitors/SPRPPRDN/SIMS/V27A/equil_new_ext.gro', '/nike-storage/passagrilli/MD_inhibitors/SPRPPRDN/SIMS/V27A/equil_new_ext.xtc', 'SPR')
V27A = V27A_1 + V27A_2
WT2HSP = calculate_dists('./WT_2HSP/ligands/WT_2HSP_SPRP/minim_new.gro', './WT_2HSP/ligands/WT_2HSP_SPRP/equil_new.xtc', 'SPR')
#SPRP = calculate_dists('./WT_2HSP/ligands/WT_2HSP_SPRP/minim_new.gro', './WT_2HSP/ligands/WT_2HSP_SPRP/equil_new.xtc', 'SPR')
#SPRC = calculate_dists('./WT_2HSP/ligands/WT_2HSP_SPRC/minim.gro', './WT_2HSP/ligands/WT_2HSP_SPRC/equil.xtc', 'SPRC')
#SPRN = calculate_dists('./WT_2HSP/ligands/WT_2HSP_SPRN/minim.gro', './WT_2HSP/ligands/WT_2HSP_SPRN/equil.xtc', 'SPRN')
x_vals = [x / 100 for x in range(0, len(WT2HSP))]
x_V27A = [x / 100 for x in range(0, len(V27A))]
#x_vals = [x / 100 for x in range(0, len(SPRP))]
#x_SPRN = [x / 100 for x in range(0, len(SPRN))]
plt.plot(x_vals, S31N, label='S31N')
plt.plot(x_V27A, V27A, label='V27A')
plt.plot(x_vals, WT2HSP, label='WT')
#plt.plot(x_vals, SPRP, label='SPRP')
#plt.plot(x_vals, SPRC, label='SPRC')
#plt.plot(x_SPRN, SPRN, label='SPRN')
leg = plt.legend(ncol=3, loc=9, fancybox=True)
leg.get_frame().set_alpha(0.5)
plt.xlabel('Time / ns')
plt.ylabel(ur'Z difference to resnum 31 COM / $\AA$')
# plt.show()
plt.savefig('SPRP_COM_longerV27A.png', dpi=300)
plt.close()

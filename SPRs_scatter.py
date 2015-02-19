from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as md
from matplotlib.mlab import griddata
from matplotlib.colors import LogNorm
import pylab

def calculate_dists(gro_file, xtc_file, ligand_name):
    u = md.Universe(gro_file, xtc_file)
    select_ligand = u.selectAtoms("resname " + ligand_name)
    select_res31 = u.selectAtoms("backbone and resnum 31")
    COM_distance_x = []
    COM_distance_y = []
    COM_distance_z = []
    print gro_file[:-4], select_ligand

    for i in u.trajectory:
        ligand_COM = select_ligand.centerOfMass()
        res31_COM = select_res31.centerOfMass()
        #COM_distance.append(res31_COM[2] - ligand_COM[2])
        COM_distance_x.append(ligand_COM[0] - res31_COM[0])
        COM_distance_y.append(ligand_COM[1] - res31_COM[1])
        COM_distance_z.append(ligand_COM[2] - res31_COM[2])
    return [COM_distance_x, COM_distance_y, COM_distance_z]

S31N = calculate_dists('./mutants/2HSP/2HSP_S31N/ligands/2HSP_S31N_SPRP/minim_new.gro', './mutants/2HSP/2HSP_S31N/ligands/2HSP_S31N_SPRP/equil_new.xtc', 'SPR')
V27A = calculate_dists('./mutants/2HSP/2HSP_V27A/ligands/2HSP_V27A_SPRP/minim_new.gro', './mutants/2HSP/2HSP_V27A/ligands/2HSP_V27A_SPRP/equil_new.xtc', 'SPR')
SPRP = calculate_dists('./WT_2HSP/ligands/WT_2HSP_SPRP/minim_new.gro', './WT_2HSP/ligands/WT_2HSP_SPRP/equil_new.xtc', 'SPR')
SPRC = calculate_dists('./WT_2HSP/ligands/WT_2HSP_SPRC/minim.gro', './WT_2HSP/ligands/WT_2HSP_SPRC/equil.xtc', 'SPRC')
SPRN = calculate_dists('./WT_2HSP/ligands/WT_2HSP_SPRN/minim.gro', './WT_2HSP/ligands/WT_2HSP_SPRN/equil.xtc', 'SPRN')
#x_vals = [x / 100 for x in range(0, len(SPRP))]
#x_S31N = [x / 100 for x in range(0, len(S31N))]
#plt.plot(x_vals, SPRN, label='SPRN WT')
data = [[], [], []]
for i, val in enumerate(SPRP[0]):
    data[0].append(SPRP[0][i])
    data[0].append(S31N[0][i])
    data[0].append(V27A[0][i])
    data[0].append(SPRC[0][i])
    data[0].append(SPRN[0][i])
    data[1].append(SPRP[1][i])
    data[1].append(S31N[1][i])
    data[1].append(V27A[1][i])
    data[1].append(SPRC[1][i])
    data[1].append(SPRN[1][i])

plt.scatter(data[0], data[1], marker='o', c='b', s=5, zorder=10, norm=True)
#leg.get_frame().set_alpha(0.5)
plt.xlabel('Time / ns')
plt.ylabel(ur'Z difference to resnum 31 COM / $\AA$')
# plt.show()
plt.savefig('SPRs_scatter.png', dpi=300)
plt.close()

#H, xedges, yedges, img = plt.hist2d(data[0], data[1], norm=LogNorm())
#extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#im = ax.imshow(H, cmap=plt.cm.jet, interpolation='spline36', norm=LogNorm(), extent=extent)
#fig.colorbar(im, ax=ax)
#plt.savefig('SPRs_hist.png', dpi=300)
#plt.close()
plt.hist2d(data[0], data[1], bins=40)
plt.colorbar()
plt.savefig('SPRs_hist2.png', dpi=300)
plt.close()

from __future__ import division
import matplotlib.pyplot as plt
import MDAnalysis as md
import numpy as np

def calculate_dists(gro_file, xtc_file):
    u = md.Universe(gro_file, xtc_file)
    select_group1 = u.selectAtoms("backbone and (resnum 50 or resnum 51)")
    select_group2 = u.selectAtoms("backbone and (resnum 149 or resnum 150)")
    select_group3 = u.selectAtoms("backbone and (resnum 50 or resnum 51 or resnum 149 or resnum 150)")
    select_group4 = u.selectAtoms("backbone and (resnum 25 or resnum 124)")
    for i in select_group1:
        print "Loop1 ", i
    for i in select_group2:
        print "Loop2 ", i
    for i in select_group4:
        print "ASP ", i
    COM_distance = []
    COM_distance_ASP = []
    COM_distance_ASP1 = []
    COM_distance_ASP2 = []
    max_dist = 0
    index = 0
    min_dist = 100
    index_min = 0
    max_dist_1 = 0
    index_1 = 0
    min_dist_1 = 100
    index_min_1 = 0
    max_dist_2 = 0
    index_2 = 0
    min_dist_2 = 100
    index_min_2 = 0
    max_dist_3 = 0
    index_3 = 0
    min_dist_3 = 100
    index_min_3 = 0
    #group1_COM = select_group1.centerOfMass()
    #group2_COM = select_group2.centerOfMass()
    #print group1_COM
    #print group2_COM
    #print np.sqrt(np.dot(group1_COM-group2_COM, group1_COM-group2_COM))
    #print np.linalg.norm(group1_COM - group2_COM)
    for i in u.trajectory:
        group1_COM = select_group1.centerOfMass()
        group2_COM = select_group2.centerOfMass()
        dist = np.linalg.norm(group1_COM - group2_COM)
        COM_distance.append(dist)
        if dist > max_dist:
            max_dist = dist
            index = i.frame
        if dist < min_dist:
            min_dist = dist
            index_min = i.frame
        group3_COM = select_group3.centerOfMass()
        group4_COM = select_group4.centerOfMass()
        dist1 = np.linalg.norm(group3_COM - group4_COM)
        COM_distance_ASP.append(dist1)
        if dist1 > max_dist_1:
            max_dist_1 = dist1
            index_1 = i.frame
        if dist1 < min_dist_1:
            min_dist_1 = dist1
            index_min_1 = i.frame
        dist2 = np.linalg.norm(group1_COM - group4_COM)
        dist3 = np.linalg.norm(group2_COM - group4_COM)
        COM_distance_ASP1.append(dist2)
        COM_distance_ASP2.append(dist3)
        if dist2 > max_dist_2:
            max_dist_2 = dist2
            index_2 = i.frame
        if dist2 < min_dist_2:
            min_dist_2 = dist2
            index_min_2 = i.frame
        if dist3 > max_dist_3:
            max_dist_3 = dist3
            index_3 = i.frame
        if dist3 < min_dist_3:
            min_dist_3 = dist3
            index_min_3 = i.frame
    print 'Max interloop distance: ', max_dist, index
    print 'Min interloop distance: ', min_dist, index_min
    print 'Max loops-ASP distance: ', max_dist_1, index_1
    print 'Min loops-ASP distance: ', min_dist_1, index_min_1
    print 'Max loop1-ASP distance: ', max_dist_2, index_2
    print 'Min loop1-ASP distance: ', min_dist_2, index_min_2
    print 'Max loop2-ASP distance: ', max_dist_3, index_3
    print 'Min loop2-ASP distance: ', min_dist_3, index_min_3
    return COM_distance, COM_distance_ASP, COM_distance_ASP1, COM_distance_ASP2

coil_distance, ASP_distance, ASP_distance1, ASP_distance2 = calculate_dists('structure.pdb', 'equ.dcd')
x_vals = [x / 10 for x in range(0, len(coil_distance))]
plt.plot(x_vals, coil_distance, linewidth=0.5)
#leg = plt.legend(ncol=3, loc=9, fancybox=True)
#leg.get_frame().set_alpha(0.5)
plt.xlabel('Time / ns')
plt.ylabel(ur'Loop COM distance / $\AA$')
plt.axhline(y=9.84, linewidth=1, color = 'red')
plt.axhline(y=11.11, linewidth=1, color = 'green')
plt.savefig('coil_COMdistance.png', dpi=300)
plt.close()
plt.plot(x_vals, ASP_distance, linewidth=0.5)
plt.plot(x_vals, ASP_distance1, linewidth=0.5)
plt.plot(x_vals, ASP_distance2, linewidth=0.5)
print 'Loop1 average: ', np.average(ASP_distance1[500:]), np.std(ASP_distance1[500:])
print 'Loop2 average: ', np.average(ASP_distance2[500:]), np.std(ASP_distance2[500:])
plt.xlabel('Time / ns')
plt.ylabel(ur'Loop COM distance / $\AA$')
plt.axhline(y=21.29, linewidth=1, color = '#C45AEC', label='PR20')
plt.axhline(y=15.18, linewidth=1, color = '#C45AEC')
plt.axhline(y=20.36, linewidth=1, color = '#EAC117', label='PR')
plt.axhline(y=15.11, linewidth=1, color = '#EAC117')
plt.axhline(y=np.average(ASP_distance1), linewidth=1, color = 'green', label='Loop1 average')
plt.axhline(y=np.average(ASP_distance2), linewidth=1, color = 'red', label='Loop2 average')
leg = plt.legend(fancybox=True, loc=2, framealpha=0.5)
#leg.get_frame().set_alpha(0.5)
plt.savefig('ASP_COMdistance.png', dpi=300)
plt.close()


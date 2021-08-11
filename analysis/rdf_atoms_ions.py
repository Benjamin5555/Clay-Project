#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 15:39:04 2021

@author: sarah
"""

# Average radial distribution function for two groups of atoms

import MDAnalysis as mda
from MDAnalysis.analysis import lineardensity as lin
from MDAnalysis.analysis import density
from MDAnalysis.analysis import rdf
from MDAnalysis.topology import tpr
from MDAnalysis.topology import TPRParser
from definitions import *
import sys
sys.path.append("../../../")
from clayAnalysis import ClayAnalysis 


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import sys


tpr = sys.argv[1] 
trr = sys.argv[2]

u=mda.Universe(tpr, trr,tpr_resid_from_one=False)

cal = ClayAnalysis(u)

surfaces = cal.generate_surface_group("waters")
w_surf_grp = cal.combine_atomgroups(surfaces[0])+cal.combine_atomgroups(surfaces[1])
print(np.unique(w_surf_grp.types))
w_surf_grp= w_surf_grp.select_atoms("type o*")
ATs = u.atoms.select_atoms("type at*")
STs = u.atoms.select_atoms("type st*")
print(w_surf_grp)
print(ATs)
print(STs)

for ions in ["resname Cs", "resname K", "resname Na"]:
    ion_grp = u.atoms.select_atoms(ions)
    print(np.unique(u.atoms.types))
    for surf_atm in [w_surf_grp,ATs,STs]: 
        average_rdf(ion_grp, surf_atm)
    
    plt.savefig(str(ions)+".png") 
    plt.show()

#for ions in ["resname Cs", "resname K", "resname Na"]:
#    ion_grp = u.atoms.select_atoms(ions)
#    print(np.unique(u.atoms.types))
#    surf_grp = u.atoms.select_atoms("resname NON*")
#    print(ion_grp.n_atoms)
#    average_rdf(ion_grp, surf_grp)
#    plt.show()
#
#clay=mda.AtomGroup(u.select_atoms('resname NON*'))
#atomgroup_coords(clay)
#position_density(tpr, trr, 'resname NON*', False)
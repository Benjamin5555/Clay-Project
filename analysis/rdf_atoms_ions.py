#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 15:39:04 2021

@author: sarah

Creates a plot of the radial distribution of clay surface atoms (ATs, STs and oxygens) to a set of 
cations (Cs, K, Na edit as necessary).

Args:
    tpr: Topology file name
    trr: Trajectory file name
    sys: describes name of system
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
sys.path.append("ClayAnalysis")

from clayAnalysis import ClayAnalysis 
import clayAnalysis

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


import sys


tpr = sys.argv[1] 
trr = sys.argv[2]
sys = sys.argv[3]

u=mda.Universe(tpr, trr,tpr_resid_from_one=False)

cal = ClayAnalysis(u)
	
surfaces = cal.generate_surface_group("waters")
w_surf_grp = cal.combine_atomgroups(surfaces[0])+cal.combine_atomgroups(surfaces[1])
print(np.unique(w_surf_grp.types))

w_surf_grp= w_surf_grp.select_atoms("type o*")

ATs = u.atoms.select_atoms("type at*")
STs = u.atoms.select_atoms("type st*")
clayAnalysis.plot_group(w_surf_grp,label="Waters")
clayAnalysis.plot_group(ATs)
clayAnalysis.plot_group(STs)
plt.show()
	
for ions in ["resname Cs", "resname K", "resname Na"]:
    ion_grp = u.atoms.select_atoms(ions)
    for surf_atm in [w_surf_grp,ATs,STs]: 
        average_rdf(ion_grp, surf_atm)
    
    plt.savefig(str(sys)+"_"+str(ions)+".png") 


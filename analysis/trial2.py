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


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import sys

print("THIS ONE")

tpr = sys.argv[1] 
trr = sys.argv[2]

u=mda.Universe(tpr, trr,tpr_resid_from_one=False)

SOL_upper=u.select_atoms('resname SOL and prop z > 89')
Cs_upper=u.select_atoms('resname Cs and prop z > 89 and prop z < 95')
clay_oxy=u.select_atoms('resname NON2 and name O* and prop z > 89')

ags=[[clay_oxy, SOL_upper], [Cs_upper, SOL_upper]]
average_rdf(Cs_upper, SOL_upper)
plt.show()


#clay=mda.AtomGroup(u.select_atoms('resname NON*'))
#atomgroup_coords(clay)
#position_density(tpr, trr, 'resname NON*', False)

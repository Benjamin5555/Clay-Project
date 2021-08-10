#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 11:11:04 2021

@author: sarah
"""

import MDAnalysis as mda
from MDAnalysis.analysis import lineardensity as lin
from MDAnalysis.analysis import density
from MDAnalysis.topology import tpr
from MDAnalysis.topology import TPRParser

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def gen_density_dataframe(atom_selection):
    # assign bulk ions AtomGroup
    print("Selection resnames found:")
    print(np.unique(atom_selection.resnames))
    # ionic z-density profiles
    denisties=lin.LinearDensity(atom_selection, binsize=0.25).run()
    
    # z-values
    z_max = u.dimensions[3]
    
    # create dataframe from profiles
    densities_df=pd.DataFrame(denisties.results['z'], columns = ['pos', 'char'])
    
    num_entries = len(densities_df['pos'])
    
    index=pd.Index(np.linspace(0, z_max, num_entries))
    
    
    densities_df=densities_df.set_index(index)
    densities_df.reset_index(inplace=True)
    densities_df=densities_df.rename(columns = {'index':'z-coordinate'})
    
    #Normalise
    densities_df['posNorm']=(densities_df['pos']-densities_df['pos'].mean())/densities_df['pos'].std()
    densities_df['charNorm']=(densities_df['char']-densities_df['char'].mean())/densities_df['char'].std()
    
    return densities_df 


####PLOTTING
def plot_density_profile(densities_df,label="",fig=None,cat_dens=None):
    if fig == None or cat_dens== None:
        fig, cat_dens=plt.subplots()
    #cat_dens.plot(densities_df['z-coordinate'], densities_df['char'], '--', color='black', label=str(label)+" charge density")
    cat_dens.plot(densities_df['z-coordinate'], densities_df['pos'], ':', label=str(label) +" mass density")
    cat_dens.set_xlabel("z-coordinate (Å)")
    cat_dens.set_ylabel("density")
    cat_dens.legend()
    return fig, cat_dens


def plot_norm_density_profile(densities_df,label="",fig=None,cat_norm=None):
    if fig ==None or cat_norm == None:
        fig, cat_norm=plt.subplots()
    cat_norm.plot(densities_df['z-coordinate'], densities_df['charNorm'], '--', color='black', label=str(label)+ " normalised charge density")
    cat_norm.plot(densities_df['z-coordinate'], densities_df['posNorm'], ':', label=str(label)+" normalised mass density")
    cat_norm.set_xlabel("z-coordinate (Å)")
    cat_norm.set_ylabel("normalised densities")
    cat_norm.legend()
    return fig, cat_norm



# create Universe
u=mda.Universe("../../TestFiles/control.tpr", "../../TestFiles/control.trr")

# assign clay AtomGroup
clay=mda.AtomGroup(u.select_atoms('resname NON*'))


dynamic_ions=mda.AtomGroup(u.select_atoms('not resname NON* and  not resname Cl and not resname SOL and not resname iSL and not resname Ca'))
ionic_density_df= gen_density_dataframe(dynamic_ions)
fig,cat_dens = plot_density_profile(ionic_density_df,"Cs")
figNorm,cat_norm = plot_norm_density_profile(ionic_density_df,"Cs")

clay_density_df= gen_density_dataframe(clay)
plot_density_profile(clay_density_df,"clay",fig,cat_dens)
plot_norm_density_profile(ionic_density_df,"clay", figNorm, cat_norm)
plt.show()

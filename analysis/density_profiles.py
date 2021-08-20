#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 11:11:04 2021

@author: sarah

Produces a graph of the mass densiy for the clay and a selection (MDAnalysis selection syntax) as 
well as a normalised version of this also.

Example (Command Line) 
    >>> python3 density_profiles top_path trr_path selection_string

"""

import MDAnalysis as mda
from MDAnalysis.analysis import lineardensity as lin
from MDAnalysis.analysis import density
from MDAnalysis.topology import tpr
from MDAnalysis.topology import TPRParser
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def gen_density_dataframe(atom_selection):
    """
        Generates mass and charge density dataframe for a given input selection

        Args:
            atom_selecion:  A collection of atoms that are to have their densities calculated

        Returns:
            Pandas dataframe containing position ('pos') and charge ('char') densities along the z
            axis of the system. Also gives normalised versions of these (posNorm and charNorm)

        Examples:
            >>> u=mda.Universe(top_path, trj_path)
            # assign clay AtomGroup
            >>> clay=mda.AtomGroup(u.select_atoms('resname NON*'))
            >>> ions=mda.AtomGroup(u.select_atoms(sel_str))
            >>> ionic_density_df= gen_density_dataframe(ions)
            >>> norm_charges = ionic_density_df["charNorm"]

    """
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
    """
        Creates a plot of the mass density profile of a selecion, given its dataframe

        Args:
            densities_df: Pandas dataframe of densities as produced by gen_densities_dataframe
            label:   Label for current selection contained within the dataframe
            fig:     (Optional) Figure object to add the plot to
            cat_norm (Optional) axes object to add the plot to
        Returns:
            fig, ax objects that can be plotted via plt.show() or added to
    """

    if fig == None or cat_dens== None:
        fig, cat_dens=plt.subplots()
    #cat_dens.plot(densities_df['z-coordinate'], densities_df['char'], '--', color='black', label=str(label)+" charge density")
    cat_dens.plot(densities_df['z-coordinate'], densities_df['pos'], ':', label=str(label) +" mass density")
    cat_dens.set_xlabel("z-coordinate (Å)")
    cat_dens.set_ylabel("density")
    cat_dens.legend()
    return fig, cat_dens


def plot_norm_density_profile(densities_df,label="",fig=None,cat_norm=None):
    """
        Creates a plot of the normalised mass density profile of a selecion, given its dataframe

        Args:
            densities_df: Pandas dataframe of densities as produced by gen_densities_dataframe
            label:   Label for current selection contained within the dataframe
            fig:     (Optional) Figure object to add the plot to
            cat_norm (Optional) axes object to add the plot to
        Returns:
            fig, ax objects that can be plotted via plt.show() or added to


    """
    if fig ==None or cat_norm == None:
        fig, cat_norm=plt.subplots()
    cat_norm.plot(densities_df['z-coordinate'], densities_df['charNorm'], '--', color='black', label=str(label)+ " normalised charge density")
    cat_norm.plot(densities_df['z-coordinate'], densities_df['posNorm'], ':', label=str(label)+" normalised mass density")
    cat_norm.set_xlabel("z-coordinate (Å)")
    cat_norm.set_ylabel("normalised densities")
    cat_norm.legend()
    return fig, cat_norm



if __name__ == "__main__":
    if(len(sys.argv)<2):
        print("USAGE: python3 density_profiles top_path trr_path selection_string")
    top_path=  sys.argv[1]
    trj_path = sys.argv[2]
    sel_str = sys.argv[3] 
    # create Universe
    u=mda.Universe(top_path, trj_path)
    
    # assign clay AtomGroup
    clay=mda.AtomGroup(u.select_atoms('resname NON*'))
    
    
    dynamic_ions=mda.AtomGroup(u.select_atoms(sel_str))
    ionic_density_df= gen_density_dataframe(dynamic_ions)
    label = ",".join(np.unique(sel_str.names))
    fig,cat_dens = plot_density_profile(ionic_density_df,label)
    figNorm,cat_norm = plot_norm_density_profile(ionic_density_df,label)
    
    clay_density_df= gen_density_dataframe(clay)
    plot_density_profile(clay_density_df,"clay",fig,cat_dens)
    plot_norm_density_profile(ionic_density_df,"clay", figNorm, cat_norm)
    plt.show()

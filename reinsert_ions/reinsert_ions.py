#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 10:59:09 2021

@author: sarah


Removes ions that are very close to the surface of the clay in the initial time step and reinserts 
them away from the surface so that there are no cations that start the simulation adsorbed (which 
may mess up adsorption times calculation)

Also acts to reorders the .gro file such that it matches that of a topology file, required for 
gromacs

Args:
    topology_path: path to topology of system
    gro_path: path to input .gro coordinates file of system where ions can be `at' surface
    out_gro_path: output .gro path with ions reinserted away from surface
    ions to move in order appear in topology: Ordering of ions in the topology file (unique i.e. Cs \n Cs \n Ca  is just Cs Ca)

Examples:
    >>> $ python3 reinsert_ions.py topo.top in.gro out.gro Cs Cl

"""

import pathlib as pl
import re
import sys
import typing as tp
import subprocess as sp
import os
import sys
import MDAnalysis as mda
import numpy as np

from numpy import random

#%%

def get_clay_sel_str(*uc_stem: str) -> mda.AtomGroup:
    """
    Generate string for selecting clay atoms in an mda.Universe
    :param *uc_stem: Residue name of clay, e.g. 'NON'
    :type *uc_stem: str
    :return: Selection str to be used in select_atoms function
    :rtype: mda.AtomGroup

    """
    if uc_stem:
        clay_sel_strs = ('', f'{uc_stem[0]}')
        print(clay_sel_strs)
    else:
        clay_sel_strs = ('not ', ' '.join(['SOL', 'iSL', 'GLY',
                                           'Ca', 'Mg', 'Na',
                                           'K', 'Cl','Cs']))
        print(clay_sel_strs)
    clay_sel_str = ('{}resname {}*'.format(*clay_sel_strs))
    return clay_sel_str

def find_clay_minmax_z_positions(clay_atomgroup: mda.AtomGroup) -> dict:
    """
    Find minium and maximum z-positions of clay atoms
    :param clay_atomgroup: Clay atomgroup selection
    :type clay_atomgroup: mda.AtomGroup
    :return: Dictionary with maximum and minimum z-positions in A
    :rtype: dict

    """
    clay_z = clay_atomgroup.atoms.positions[:, 2]
    z_dict = {'max': np.max(clay_z),
              'min': np.min(clay_z)}
    return z_dict

def get_atoms_in_outside_clay_sel_str(seltype,
                                      selname: list,
                                      z_minmax_dict: dict,
                                      where,
                                      sel_dist: float = 0,
                                      ) -> mda.AtomGroup:
    """
    Select atoms, within sel_dist in A, outside or inside of clay
    :param seltype: selection by name or resname
    :type seltype: ['name', 'resname']
    :param selname: list of str of ions to be selected
    :type selname: str, list
    :param z_minmax_dict: Dictionary with minimum and maximum clay z-positions
    :type z_minmax_dict: dict
    :param where: select inside or outside clay
    :type where: ['inside', 'outside']
    :param sel_dist: select within a certain distance from clay, defaults to 0
    :type sel_dist: float
    :return: Atom selection within a certain distance inside or outside of clay
    :rtype: mda.AtomGroup

    """
    inout_dict = {'inside': 'not',
                  'outside': ''
                  }
    inoutsel = inout_dict[where]
    if type(selname) == list:
        selname = ' '.join(selname)
    atom_selection_str = (f'{seltype} {selname} and {inoutsel} '
                          f'((prop z > {z_minmax_dict["max"]} '
                          f' and prop z < {z_minmax_dict["max"] + sel_dist}) or '
                          f'prop z < {z_minmax_dict["min"]} '
                          f' and prop z > {z_minmax_dict["min"] - sel_dist})')
    return atom_selection_str

def select_atoms_in_outside_clay(u: mda.Universe,
                                 seltype,
                                 selname: str,
                                 where,
                                 clay_dist: float=0,
                                 *uc_stem) -> mda.AtomGroup:
    """
    Select atoms inside or outside of clay (within a distance specified by clay-dist)
    :param u: Universe to insert into
    :type u: mda.Universe
    :param seltype: Select by name or resname
    :type seltype: ['name', 'resname']
    :param selname: list or str for atom selection
    :type selname: str, list
    :param where: Select inside or outside of clay boundaries
    :type where: ['inside', 'outside']
    :param clay_dist: Select within distance in A from clay (if outside)
    :type clay_dist: float
    :param *uc_stem: Residue name of clay, e.g. 'NON', optional
    :type *uc_stem: str
    :return: Selected atoms
    :rtype: mda.AtomGroup

    """
    if where=='inside' and clay_dist > 0:
        raise ValueError('clay_dist cannot be > 0 for selecting atoms inside clay.')
    clay_sel_str = get_clay_sel_str(*uc_stem)
    clay = u.select_atoms(clay_sel_str)
    clay_z_minmax_dict = find_clay_minmax_z_positions(clay)
    atom_sel_str = get_atoms_in_outside_clay_sel_str(seltype,
                                                     selname,
                                                     clay_z_minmax_dict,
                                                     where,
                                                     clay_dist)
    atom_selection = u.select_atoms(atom_sel_str)
    return atom_selection

#%%

def get_atom_type_numbers_dict(ag: mda.AtomGroup) -> dict:
    """
    Generate dictionary with atom types and numbers in Atom Group.
    :param ag: AtomGroup for dictionary construction
    :type ag: mda.AtomGroup
    :return: Dictionary mapping atom name to number
    :rtype: dict

    """
    atom_type_names_array = np.unique(ag.atoms.names)
    atom_type_numbers_list = list(map(lambda i: ag.select_atoms(f'name {i}').n_atoms,
                                      atom_type_names_array))
    atom_types_dict = dict(zip(atom_type_names_array, atom_type_numbers_list))
    
    return atom_types_dict


#%%

def remove_atoms_from_grofile(u: mda.Universe, ag: mda.AtomGroup,
                              outgro: str) -> None:
    """
    Remove atoms in atom selection from a universe and write to .gro file.
    :param u: Universe with atoms in atom group
    :type u: mda.Universe
    :param ag: atom group to remove
    :type ag: mda.AtomGroup
    :param outgro: output filename
    :type outgro: str
    :return: None
    :rtype: None

    """
    new_u = u.atoms - ag
    new_u.write(outgro)


#%%

def write_atom_pos_datfile(atom_group: mda.AtomGroup,
                           *single_atom_type: str) -> None:
    """
    Writes positions of mda.AtomGroup
    :param atom_group: Atom group for positions
    :type atom_group: mda.AtomGroup
    :param single_atom_type: Atom type inside atom group
    :type single_atom_type: str
    :return: None
    :rtype: None

    """
    if single_atom_type:
        atom_group = atom_group.select_atoms(f'name {single_atom_type}')
        outname = f'{single_atom_type}.dat'
    else:
        outname = 'atom_pos.dat'
    pos_in_nm = np.round(atom_group.positions / 10, 4)
    np.savetxt(outname, pos_in_nm, delimiter=' ', fmt='%.4f')


#%%

def get_file_numlines(filepath: str) -> int:
    """
    Count number of lines in a file, equals number of atoms in position .dat file
    :param filepath: Path to input file
    :type filepath: str
    :return: number of lines in file
    :rtype: int

    """
    with open(filepath, 'r') as file:
        lines = sum(1 for _ in file)
    return lines


def execute_bash_alias(command, **outputargs):
    output = sp.run(['/bin/bash', '-i', '-c', *command], **outputargs)
    return output


def run_gmx_insert_molecules(confgro: str, insertgro: str, outgro: str,
                             nmols: int, pos: str='', scale: float='',
                             radius: float='', dr: float='') -> None:
    """
    
    :param confgro: configuration wuhere molecule is inserted
    :type confgro: str
    :param insertgro: configuration to be inserted
    :type insertgro: str
    :param outgro: configuration with inserted molecules
    :type outgro: str
    :param nmols: number of times insertgro is added
    :type nmols: int
    :param pos: path to .dat file with insertion positions, defaults to ''
    :type pos: str, optional
    :param scale: scaling factor for vdW radius, defaults to ''
    :type scale: float, optional
    :param radius: overwrites vdW radius, defaults to ''
    :type radius: float, optional
    :param dr: accepted deviation from specified insert positions, defaults to ''
    :type dr: float, optional
    :return: None
    :rtype: None

    """
    if pos != '':
        pos = f' -ip {pos}'
    if scale != '':
        scale = f' -scale {scale}'
    if radius != '':
        radius = f' -radius {radius}'
    if dr != '':
        radius = f' -dr {dr}'
    print(f'gmx insert-molecules -f {confgro} -ci {insertgro}'
                       f'{pos} -o {outgro} -nmol {nmols}{scale}{radius}')
    execute_bash_alias([f'gmx insert-molecules -f {confgro} -ci {insertgro}'
                       f'{pos} -o {outgro} -nmol {nmols}{scale}{radius}'])


def add_ion_replacement_waters(grofile: str, outgro: str, positions: str,
                               nmols: int, scale: float=0.35,
                               radius: float =0.105) -> None:
    """
    Insert SPC water in specified positions
    :param grofile: configuration where water is inserted
    :type grofile: str
    :param outgro: configuration with inserted waters
    :type outgro: str
    :param positions: path to .dat file with insertion positions
    :type positions: str
    :param nmols: number of waters added
    :type nmols: int
    :param scale: scaling factor for vdW radius, defaults to 0.05
    :type scale: float, optional
    :param radius: overwrites vdW radius, defaults to ''
    :type radius: float, optional
    :return: None
    :rtype: None

    """
    run_gmx_insert_molecules(grofile, 'scripts/reinsert_ions/spc_mol.gro', outgro, nmols,
                             pos=positions, scale=scale, radius=radius)


def insert_ion_waters(ion_datfile: str, ingro: str, outgro: str) -> None:
    """
    Call gmx insert-molecules to insert waters in positions specified in .dat file
    :param ion_datfile: insertion positions file
    :type ion_datfile: str
    :param ingro: configuration where waters should be inserted
    :type ingro: str
    :param outgro: configuration with inserted waters
    :type outgro: str
    :return: None
    :rtype: None

    """
    nmols = get_file_numlines(ion_datfile)
    add_ion_replacement_waters(ingro, outgro, ion_datfile,
                                   nmols)

#%%

def insert_molecules_centered(ingro: str, outgro: str, atom_dict: dict,
                              clay_dist=0) -> None:
    """
    Replace SOL by nmols ions_type within a layer at clay_dist from clay.
    :param u: Universe for insertion
    :type u: mda.Universe
    :param outgro: output gro file
    :type outgro: str
    :param clay_dist: Distance in A from clay for insertion, defaults to 0
    :type clay_dist: TYPE, optional
    :param ion_dict: dictionary with atom types and numbers
    :type ion_dict: dict
    :return: None
    :rtype: None

    """
    
    u=mda.Universe(ingro)
    
    i=0
    for pair in atom_dict:
        nmols=list(atom_dict.values())[i]
        ion_type=list(atom_dict.keys())[i]
        
        close_sol=select_atoms_in_outside_clay(u, 'resname', 'SOL', 'outside',
                                                30).residues
        all_sol=u.select_atoms('resname SOL').residues
        
        far_sol_layer=all_sol - close_sol
        
        sol_indices = far_sol_layer.resids
        sub_index_list = random.choice(sol_indices, nmols,
                                   replace=False).astype(str).tolist()
        print(sub_index_list)
        
        sol_sub_atoms = far_sol_layer.atoms.select_atoms(f'resid {" ".join(sub_index_list)}')
        sol_h_atoms = sol_sub_atoms.select_atoms('name HW*')
        print(len(sol_h_atoms))
        
        sol_o_atoms = sol_sub_atoms - sol_h_atoms
        print(sol_o_atoms)
        print(sol_o_atoms.names)
        
        sol_o_atoms.names = np.full(nmols, ion_type)
        print(sol_o_atoms.names)
        sol_o_atoms.residues.resnames = np.full(nmols, ion_type)
        new_u = u.atoms - sol_h_atoms
        u=new_u
        i += 1
  

    u.write(outgro, reindex=True)

#%%

def substitute_bulk_ions_gro(u: mda.Universe, ion_list: list, outgro: str,
                         clay_dist: float=30) -> None:
    """
    Reinsert bulk ions from within a dist_around clay in A in the centre of SOL.
    :type u: mda.Universe
    :param ion_list: ion name(s)
    :type ion_list: list, str
    :param clay_dist: Distance within which ions will be reinserted
    :type clay_dist: float
    :param outgro: output conformation file name
    :type outgro: str
    :return: None
    :rtype: None

    """
    if type(ion_list) == str:
        ion_list = [ion_list]
        
    #select ions within clay_dist of clay   
    reposition_ions = select_atoms_in_outside_clay(u, 'name', ion_list,
                                                   'outside', clay_dist)
    
    #make dictionary containing atom type and number of ions in reposition_ions
    ions_dict = get_atom_type_numbers_dict(reposition_ions)
    
    #remove ions in reposition_ions from universe and .gro file
    remove_atoms_from_grofile(u, reposition_ions, outgro)
    
    #writes .dat file containing the positions of the ions in reposition_ions
    write_atom_pos_datfile(reposition_ions)
    
    #calls gmx insert_molecules to insert SPC water into position ions were removed from
    insert_ion_waters('atom_pos.dat', outgro, outgro)

    #replace SOL with ions
    insert_molecules_centered(outgro, outgro, ions_dict)
    

#%%

def run_gmx_editconf(confgro: str) -> None:
    """
    
    :param confgro: configuration where mresidues are reordered
    :type confgro: str
    """
    
    print(f'gmx editconf -f {confgro} -o {confgro} -resnr 1')
    execute_bash_alias([f'gmx editconf -f {confgro} -o {confgro} -resnr 1'])


def reindex_resnr(grofile: str) -> None:
    """
    Reindex the residue numbers in grofile
    :param grofile: grofile to reorder
    :type grofile: str
    :return: None
    :rtype: None

    """
    run_gmx_editconf(grofile)


def reorder_gro(ingro, reorder_list):
    '''
    Selects atoms in reorder_list above and below clay, removes atom selection from overall gro file and adds atom selection to bottom of gro file 
    
    Parameters
    -------------
    ingro: .gro file you wish to reorder
    
    reorder_list: residues you wish to reorder, in the order you want them to appear in the .gro
    '''
    
    i=0
    for entry in reorder_list:
        u=mda.Universe(ingro)
        
        ion_type=reorder_list[i]
        ag=select_atoms_in_outside_clay(u, 'resname', ion_type, 'outside', 100)
        
        
        new_u=u.atoms - ag.atoms
        new_u2=new_u.atoms + ag.atoms
        
        u=new_u2
        u.write(ingro)
        i += 1
   
    reindex_resnr(ingro)       



if __name__ == "__main__":
	if(len(sys.argv)<2):
		print("Usage:\n topology_path, gro_path, out_gro_path, ions to move in order appear in topology")
		print("e.g. $ python3 reinsert_ions.py topo.top in.gro out.gro Cs Cl")
	
	topo_path = sys.argv[1]
	gro_path = sys.argv[2]
	out_gro_path =sys.argv[3]
	ions_to_move = sys.argv[4:]	

	print(topo_path,gro_path)
	u = mda.Universe(topo_path,gro_path,topology_format='ITP')

	substitute_bulk_ions_gro(u,ions_to_move,out_gro_path)
	
	reorder_gro(out_gro_path,ions_to_move)

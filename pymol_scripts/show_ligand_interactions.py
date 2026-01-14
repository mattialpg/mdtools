#show_ligand_interactions v.1.0
# author: Thomas Evangelidis, 2019
# License: BSD-2-Clause

from pymol import cmd, util

def show_ligand_interactions(recsel='not hetatm', ligsel='org', cutoff=4):
	"""
DESCRIPTION

	Visualize interactions between receptor and ligand.

ARGUMENTS

	recsel = string: atom selection of the receptor
	ligsel = string: atom selections of the ligand
	cutoff = float: show as sticks all receptor residues within this distance from the ligand {default: 4.0}
	"""

	# cmd.set('cartoon_transparency', 0.2)
	cmd.set('dash_radius', 0.1)
	# cmd.show("sticks", "(hydro)");

	cmd.select('ligand', ligsel)
	cmd.select('receptor', recsel)
	cmd.select('pocket', 'byres (receptor within %s of ligand)' % cutoff);
	cmd.show('sticks', 'pocket')
	# cmd.hide('(h. and (e. c extend 1))')
	
	# Find H-bonds
	cmd.set('h_bond_max_angle', 30)
	cmd.set('h_bond_cutoff_center', 3.6)
	cmd.set('h_bond_cutoff_edge', 3.2)
	cmd.dist('ligand_Hbonds', 'ligand', 'receptor', 3.8, mode=2)

	# now set the label options
	cmd.set('label_size', 20)
	cmd.set('label_position', [0,0,10])

cmd.extend('show_ligand_interactions', show_ligand_interactions)

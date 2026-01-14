# -*- coding: utf-8 -*-
from pymol import cmd, stored, CmdException
from pymol.cgo import *
from random import randint  
from glob import glob
from os.path import sep, basename

def del_noligand(selection='sele'):
    for obj in cmd.get_object_list(selection):
        cmd.select('lig', obj + ' & org')
        if cmd.count_atoms('lig') == 0:
            cmd.delete(obj)
            print('*** ' + obj + ' deleted ***')
    cmd.delete('lig')
cmd.extend('del_noligand', del_noligand)

def group_noligand(selection='sele'):
    for obj in cmd.get_object_list(selection):
        cmd.select('lig', obj + ' & org')
        if cmd.count_atoms('lig') == 0:
            cmd.group('apo', obj)
            print('*** ' + obj + ' moved to apo group ***')
cmd.extend('group_noligand', group_noligand)

def export_as(extension, selection='sele'):
    for obj in cmd.get_object_list(selection):
        filename = str(obj + '.' + extension)
        cmd.save(filename, obj)
cmd.extend('export_as', export_as)

def to_group(group_name, selection='sele'):
    for l in cmd.get_object_list(selection):
        cmd.group(group_name, l)
cmd.extend('to_group', to_group)

def prepare_medoid():
    biolip_list = ['ACE', 'HEX', 'TMA', 'SOH', 'P25', 'CCN', 'PR', 'PTN', 'NO3', 'TCN', 'BU1', 'BCN', 'CB3', 'HCS', 'NBN',\
                   'SO2', 'MO6', 'MOH', 'CAC', 'MLT', 'KR', '6PH', 'MOS', 'UNL', 'MO3', 'SR', 'CD3', 'PB', 'ACM', 'LUT',\
                   'PMS', 'OF3', 'SCN', 'DHB', 'E4N', '13P', '3PG', 'CYC', 'NC', 'BEN', 'NAO', 'PHQ', 'EPE', 'BME', 'TB',\
                   'ETE', 'EU', 'OES', 'EAP', 'ETX', 'BEZ', '5AD', 'OC2', 'OLA', 'GD3', 'CIT', 'DVT', 'OC6', 'MW1', 'OC3',\
                   'SRT', 'LCO', 'BNZ', 'PPV', 'STE', 'PEG', 'RU', 'PGE', 'MPO', 'B3P', 'OGA', 'IPA', 'LU', 'EDO', 'MAC',\
                   '9PE', 'IPH', 'MBN', 'C1O', '1PE', 'YF3', 'PEF', 'GD', '8PE', 'DKA', 'RB', 'YB', 'GGD', 'SE4', 'LHG',\
                   'SMO', 'DGD', 'CMO', 'MLI', 'MW2', 'DTT', 'DOD', '7PH', 'PBM', 'AU', 'FOR', 'PSC', 'TG1', 'KAI', '1PG',\
                   'DGA', 'IR', 'PE4', 'VO4', 'ACN', 'AG', 'MO4', 'OCL', '6UL', 'CHT', 'RHD', 'CPS', 'IR3', 'OC4', 'MTE',\
                   'HGC', 'CR', 'PC1', 'HC4', 'TEA', 'BOG', 'PEO', 'PE5', '144', 'IUM', 'LMG', 'SQU', 'MMC', 'GOL', 'NVP',\
                   'AU3', '3PH', 'PT4', 'PGO', 'ICT', 'OCM', 'BCR', 'PG4', 'L4P', 'OPC', 'OXM', 'SQD', 'PQ9', 'BAM', 'PI',\
                   'PL9', 'IRI', '15P', 'MAE', 'MBO', 'FMT', 'L1P', 'DUD', 'PGV', 'CD1', 'P33', 'DTU', 'XAT', 'CD',\
                   'THE', 'U1', 'NA', 'MW3', 'BHG', 'Y1', 'OCT', 'BET', 'MPD', 'HTO', 'IBM', 'D01', 'HAI', 'HED', 'CAD',\
                   'CUZ', 'TLA', 'SO4', 'OC5', 'ETF', 'MRD', 'PT', 'PHB', 'URE', 'MLA', 'TGL', 'PLM', 'NET', 'LAC', 'AUC',\
                   'UNX', 'GA', 'DMS', 'MO2', 'LA', 'NI', 'TE', 'THJ', 'NHE', 'HAE', 'MO1', 'DAO', '3PE', 'LMU', 'DHJ',\
                   'FLC', 'SAL', 'GAI', 'ORO', 'HEZ', 'TAM', 'TRA', 'NEX', 'CXS', 'LCP', 'OCN', 'PER', 'ACY', 'MH2',\
                   'ARS', '12P', 'L3P', 'PUT', 'IN', 'CS', 'NAW', 'SB', 'GUN', 'SX', 'CON', 'C2O', 'EMC', 'BO4', 'BNG',\
                   'MN5', '__O', 'K', 'CYN', 'H2S', 'MH3', 'YT3', 'P22', 'KO4', '1AG', 'CE', 'IPL', 'PG6', 'MO5', 'F09',\
                   'HO', 'AL', 'TRS', 'EOH', 'GCP', 'MSE', 'AKR', 'NCO', 'PO4', 'L2P', 'LDA', 'SIN', 'DMI', 'SM', 'DTD',\
                   'SGM', 'DIO', 'PPI', 'DDQ', 'DPO', 'HCA', 'CO5', 'PD', 'OS', 'OH', 'NA6', 'NAG', 'W', 'ENC', 'NA5',\
                   'LI1', 'P4C', 'GLV', 'DMF', 'ACT', 'BTB', '6PL', 'BGL', 'OF1', 'N8E', 'LMT', 'THM', 'EU3', 'PGR', 'NA2',\
                   'FOL', '543', '_CP', 'PEK', 'NSP', 'PEE', 'OCO', 'CHD', 'CO2', 'TBU', 'UMQ', 'MES', 'NH4', 'CD5', 'HTG',\
                   'DEP', 'OC1', 'KDO', '2PE', 'PE3', 'IOD', 'NDG', 'CL', 'HG', 'F', 'XE', 'TL', 'BA', 'LI', 'BR', 'TAU',\
                   'TCA', 'SPD', 'SPM', 'SAR', 'SUC', 'PAM', 'SPH', 'BE7', 'P4G', 'OLC', 'OLB', 'LFA', 'D10', 'D12', 'DD9',\
                   'HP6', 'R16', 'PX4', 'TRD', 'UND', 'FTT', 'MYR', 'RG1', 'IMD', 'DMN', 'KEN', 'C14', 'UPL', 'CMJ', 'ULI',\
                   'MYS', 'TWT', 'M2M', 'P15', 'PG0', 'PEU', 'AE3', 'TOE', 'ME2', 'PE8', '6JZ', '7PE', 'P3G', '7PG', 'PG5',\
                   '16P', 'XPE', 'PGF', 'AE4', '7E8', '7E9', 'MVC', 'TAR', 'DMR', 'LMR', 'NER', '02U', 'NGZ', 'LXB', 'A2G',\
                   'BM3', 'NAA', 'NGA', 'LXZ', 'PX6', 'PA8', 'LPP', 'PX2', 'MYY', 'PX8', 'PD7', 'XP4', 'XPA', 'PEV', '6PE',\
                   'PEX', 'PEH', 'PTY', 'YB2', 'PGT', 'CN3', 'AGA', 'DGG', 'CD4', 'CN6', 'CDL', 'PG8', 'MGE', 'DTV', 'L44',\
                   'L2C', '4AG', 'B3H', '1EM', 'DDR', 'I42', 'CNS', 'PC7', 'HGP', 'PC8', 'HGX', 'LIO', 'PLD', 'PC2', 'PCF',\
                   'MC3', 'P1O', 'PLC', 'PC6', 'HSH', 'BXC', 'HSG', 'DPG', '2DP', 'POV', 'PCW', 'GVT', 'CE9', 'CXE', 'C10',\
                   'CE1', 'SPJ', 'SPZ', 'SPK', 'SPW', 'HT3', 'HTH', '2OP', '3NI', 'BO3', 'DET', 'D1D', 'SWE', 'SOG']

    # Remove extra chains and non-LOI molecules 
    cmd.remove('solvent + inorganic')
    # cmd.remove('visible and not chain A')
    for b in biolip_list: cmd.remove('resn ' + b)
    
    # Group chains with no bound ligand
    # group_noligand()
    
    # Align to the first object
    obj_list = cmd.get_object_list()
    for obj in obj_list[1:]:
        cmd.align(obj, obj_list[0])
cmd.extend ('prepare_medoid', prepare_medoid)

def show_ligand_interactions(recsel='not hetatm', ligsel='org', cutoff=4):
    """
DESCRIPTION
    Visualize interactions between receptor and ligand.

ARGUMENTS
    recsel = string: atom selection of the receptor
    ligsel = string: atom selections of the ligand
    cutoff = float: select residues in the pocket within this distance from the ligand
    """

    # cmd.set('cartoon_transparency', 0.2)
    cmd.set('dash_radius', 0.15)
    # cmd.show("sticks", "(hydro)")

    cmd.select('ligand', ligsel)
    cmd.select('receptor', recsel)
    cmd.select('pocket', 'byres (receptor within %s of ligand)' % cutoff)
    # cmd.show('sticks', 'pocket')
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

### EDIT MOLECULAR/OBJ PROPERTIES ###
def renumber(selection, atom='CA', offset=1):
    stored.resi1 = []
    cmd.iterate(selection + ' and name ' + atom, 'stored.resi1.append(resi)')
    for i,r in enumerate(stored.resi1):
        cmd.alter("resi %s" %r, "resi=%d" %(i+offset))
    cmd.sort()
cmd.extend("renumber", renumber)

def count_molecules(selection='sele'):
    """
    Returns the number of distinct molecules (residues)
    in a given selection.
    """
    stored.resi1 = []
    cmd.iterate(selection, 'stored.resi1.append(resi)')
    print(len(set(stored.resi1)))
cmd.extend('count_molecules', count_molecules)

def move_to(GroupName):
    for obj in cmd.get_object_list('sele'):
        cmd.group(GroupName, members=obj)
cmd.extend('move_to', move_to)

def exchange_atoms(atom1, atom2, selection=None):
    """
    Exchange atom1 with atom2.
    If no selection is given, the function works on all the residues 
    """
    # Get residue numbers
    stored.resi1 = []
    if selection is None: selection = atom1
    cmd.iterate(selection, 'stored.resi1.append(resi)')
    stored.resi1 = list(set(stored.resi1)) # Keep unique values
    
    for r in stored.resi1:
        #Get atom ids
        id1 = cmd.id_atom("resi %s and %s" %(r, atom1))
        id2 = cmd.id_atom("resi %s and %s" %(r, atom2))
        # Get coordinates
        xyz_1 = cmd.get_coords("id %s" %id1, 1)
        xyz_2 = cmd.get_coords("id %s" %id2, 1)
        # Exchange coordinates
        cmd.alter_state(1, "id %s" %id1, "(x,y,z)=(%s,%s,%s)" %(xyz_2[0][0],xyz_2[0][1],xyz_2[0][2]))
        cmd.alter_state(1, "id %s" %id2, "(x,y,z)=(%s,%s,%s)" %(xyz_1[0][0],xyz_1[0][1],xyz_1[0][2]))
        
        pairs = []
        # Get list of all bonds
        model = cmd.get_model("resi %s" %r)
        for bond in model.bond: pairs.append([model.atom[bond.index[0]].id, model.atom[bond.index[1]].id])
        # Extract bonds of atom1 and atom2
        bonds_1 = [x for x in pairs if id1 in x]
        bonds_2 = [x for x in pairs if id2 in x]
        # Unbond exchanged atoms
        for b in bonds_1: cmd.unbond('id ' + str(b[0]), 'id ' + str(b[1]))
        for b in bonds_2: cmd.unbond('id ' + str(b[0]), 'id ' + str(b[1]))
        # Rebond exchanged atoms in the new position
        for b in bonds_1:
            b.remove(id1)
            cmd.bond('id ' + str(id2), 'id ' + str(b[0]))
        for b in bonds_2:
            b.remove(id2)
            cmd.bond('id ' + str(id1), 'id ' + str(b[0]))
        print("Atom %s exchanged with atom %s (resi %s)" %(id1,id2,r))
    print("Exchanged %i atom pairs" %(len(stored.resi1)))
cmd.extend('exchange_atoms', exchange_atoms)

def group_by_auth(dir='.'):
    file_list = glob(dir + '/*pdb')
    for file in file_list:
        cmd.load(file)
        name = file.replace('/', '\\').split('\\')[-1].replace('.pdb','')
        x = name.split('_')
        cmd.group(x[0], members=name)
        # cmd.set_name(file, x[1])
cmd.extend('group_by_auth', group_by_auth)

    
    
    
    

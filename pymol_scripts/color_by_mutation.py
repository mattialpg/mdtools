'''
Created by Christoph Malisi
Improved by Mattia Lo Piccolo

Creates an alignment of two proteins and superimposes them. 
Aligned residues that are different in the two (i.e. mutations) are highlighted and 
colored according to their difference in the BLOSUM90 matrix. 
Is meant to be used for similar proteins, e.g. close homologs or point mutants, 
to visualize their differences.

'''

from pymol import cmd

# Amino acids
standard_amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
special_amino_acids = ['B', 'Z', 'X', '*']  # B=Asx, Z=Glx, X=unknown, *=stop
aa_dict = {aa: idx for idx, aa in enumerate(standard_amino_acids + special_amino_acids)}

#             A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
blosum90 = [[ 5, -2, -2, -3, -1, -1, -1,  0, -2, -2, -2, -1, -2, -3, -1,  1,  0, -4, -3, -1, -2, -1, -1, -6], 
            [-2,  6, -1, -3, -5,  1, -1, -3,  0, -4, -3,  2, -2, -4, -3, -1, -2, -4, -3, -3, -2,  0, -2, -6], 
            [-2, -1,  7,  1, -4,  0, -1, -1,  0, -4, -4,  0, -3, -4, -3,  0,  0, -5, -3, -4,  4, -1, -2, -6], 
            [-3, -3,  1,  7, -5, -1,  1, -2, -2, -5, -5, -1, -4, -5, -3, -1, -2, -6, -4, -5,  4,  0, -2, -6], 
            [-1, -5, -4, -5,  9, -4, -6, -4, -5, -2, -2, -4, -2, -3, -4, -2, -2, -4, -4, -2, -4, -5, -3, -6], 
            [-1,  1,  0, -1, -4,  7,  2, -3,  1, -4, -3,  1,  0, -4, -2, -1, -1, -3, -3, -3, -1,  4, -1, -6], 
            [-1, -1, -1,  1, -6,  2,  6, -3, -1, -4, -4,  0, -3, -5, -2, -1, -1, -5, -4, -3,  0,  4, -2, -6], 
            [ 0, -3, -1, -2, -4, -3, -3,  6, -3, -5, -5, -2, -4, -5, -3, -1, -3, -4, -5, -5, -2, -3, -2, -6], 
            [-2,  0,  0, -2, -5,  1, -1, -3,  8, -4, -4, -1, -3, -2, -3, -2, -2, -3,  1, -4, -1,  0, -2, -6], 
            [-2, -4, -4, -5, -2, -4, -4, -5, -4,  5,  1, -4,  1, -1, -4, -3, -1, -4, -2,  3, -5, -4, -2, -6], 
            [-2, -3, -4, -5, -2, -3, -4, -5, -4,  1,  5, -3,  2,  0, -4, -3, -2, -3, -2,  0, -5, -4, -2, -6], 
            [-1,  2,  0, -1, -4,  1,  0, -2, -1, -4, -3,  6, -2, -4, -2, -1, -1, -5, -3, -3, -1,  1, -1, -6], 
            [-2, -2, -3, -4, -2,  0, -3, -4, -3,  1,  2, -2,  7, -1, -3, -2, -1, -2, -2,  0, -4, -2, -1, -6], 
            [-3, -4, -4, -5, -3, -4, -5, -5, -2, -1,  0, -4, -1,  7, -4, -3, -3,  0,  3, -2, -4, -4, -2, -6], 
            [-1, -3, -3, -3, -4, -2, -2, -3, -3, -4, -4, -2, -3, -4,  8, -2, -2, -5, -4, -3, -3, -2, -2, -6], 
            [ 1, -1,  0, -1, -2, -1, -1, -1, -2, -3, -3, -1, -2, -3, -2,  5,  1, -4, -3, -2,  0, -1, -1, -6], 
            [ 0, -2,  0, -2, -2, -1, -1, -3, -2, -1, -2, -1, -1, -3, -2,  1,  6, -4, -2, -1, -1, -1, -1, -6], 
            [-4, -4, -5, -6, -4, -3, -5, -4, -3, -4, -3, -5, -2,  0, -5, -4, -4, 11,  2, -3, -6, -4, -3, -6], 
            [-3, -3, -3, -4, -4, -3, -4, -5,  1, -2, -2, -3, -2,  3, -4, -3, -2,  2,  8, -3, -4, -3, -2, -6], 
            [-1, -3, -4, -5, -2, -3, -3, -5, -4,  3,  0, -3,  0, -2, -3, -2, -1, -3, -3,  5, -4, -3, -2, -6], 
            [-2, -2,  4,  4, -4, -1,  0, -2, -1, -5, -5, -1, -4, -4, -3,  0, -1, -6, -4, -4,  4,  0, -2, -6], 
            [-1,  0, -1,  0, -5,  4,  4, -3,  0, -4, -4,  1, -2, -4, -2, -1, -1, -4, -3, -3,  0,  4, -1, -6], 
            [-1, -2, -2, -2, -3, -1, -2, -2, -2, -2, -2, -1, -1, -2, -2, -1, -1, -3, -2, -2, -2, -1, -2, -6], 
            [-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,  1]] 

def getBlosum90ColorName(aa1, aa2):
    '''returns a rgb color name of a color that represents the similarity of the two residues according to
       the BLOSUM90 matrix. the color is on a spectrum from blue to red, where blue is very similar, and 
       red very disimilar.'''
    # return red for residues that are not part of the 20 amino acids
    if aa1 not in aa_dict or aa2 not in aa_dict:
        return 'red'
    
    # If the two are the same, return green
    if aa1 == aa2:
        return 'green'
    i1 = aa_dict[aa1]
    i2 = aa_dict[aa2]
    b = blosum90[i1][i2]

    # 3 is the highest score for non-identical substitutions, so substract 4 to get into range [-10, -1]
    b = abs(b - 4)

    # map to (0, 1]:
    b = 1. - (b / 10.0)

    bcolor = (1.0, b, 0.0)
    col_name = '0x%02x%02x%02x' % tuple(int(b * 0xFF) for b in bcolor)
    return col_name


import networkx as nx
import string

def assign_chains(obj_name):
    model = cmd.get_model(obj_name)

    # Create the graph
    G = nx.Graph()

    # Add atoms as nodes
    for i, atom in enumerate(model.atom):
        G.add_node(i, name=atom.name, resn=atom.resn, resi=atom.resi, id=atom.id)

    # Use the bond matrix to add edges
    for bond in model.bond:
        idx1 = bond.index[0]
        idx2 = bond.index[1]
        G.add_edge(idx1, idx2, order=bond.order)

    # Connected components
    components = list(nx.connected_components(G))
    # print(f"\nTotal disconnected subgraphs: {len(components)}\n")

    chain_ids = string.ascii_uppercase[:26]
    for i, comp in enumerate(components):
        id_set = {G.nodes[node]['id'] for node in comp}
        id_selection = "id " + "+".join(str(x) for x in id_set)
        cmd.alter(f"{obj_name} and {id_selection}", f"chain='{chain_ids[i]}'")
    cmd.rebuild(obj_name)

def color_by_mutation(obj1, obj2, waters=0, labels=0):
    '''  
    USAGE: color_by_mutation selection1, selection2 [,waters [,labels ]]
    '''
    from pymol import stored, CmdException

    if cmd.count_atoms(obj1) == 0:
        print('%s is empty'%obj1)
        return
    if cmd.count_atoms(obj2) == 0:
        print('%s is empty'%obj2)
        return
    waters = int(waters)
    labels = int(labels)

    assign_chains(obj1)
    assign_chains(obj2)
    
    # Align the proteins
    cmd.super(obj1, obj2, object='aln', cycles=0)
    cmd.super(obj1, obj2)

    stored.resn1, stored.resn2 = [], []
    stored.resi1, stored.resi2 = [], []
    stored.chain1, stored.chain2 = [], []
    
    # Store residue ids, residue names and chains of aligned residues
    cmd.iterate(f'{obj1} and name CA and aln', 'stored.resn1.append(resn)')
    cmd.iterate(f'{obj1} and name CA and aln', 'stored.resi1.append(resi)')
    cmd.iterate(f'{obj2} and name CA and aln', 'stored.resi2.append(resi)')
    cmd.iterate(f'{obj2} and name CA and aln', 'stored.resn2.append(resn)')
    cmd.iterate(f'{obj1} and name CA and aln', 'stored.chain1.append(chain)')
    cmd.iterate(f'{obj2} and name CA and aln', 'stored.chain2.append(chain)')

    mutant_selection = []
    non_mutant_selection = []
    colors = []
    for n1, n2, i1, i2, c1, c2 in zip(stored.resn1, stored.resn2,
        stored.resi1, stored.resi2, stored.chain1, stored.chain2):

        sel1 = f'{obj1} and resi {i1} and chain {c1}'
        sel2 = f'{obj2} and resi {i2} and chain {c2}'

        if n1 == n2:
            non_mutant_selection += [sel1, sel2]
            color = 'limegreen'
            # print(n1, n2, sel1, sel2)
        else:
            mutant_selection += [sel1, sel2]
            color = getBlosum90ColorName(n1, n2)
            # color = 'blue'
        colors.append((color, sel1))
        colors.append((color, sel2))

    if not mutant_selection:
        raise CmdException('No mutations found')

    cmd.select('mutations', ' or '.join(mutant_selection))
    cmd.select('non_mutations', ' or '.join(non_mutant_selection))
    cmd.select('not_aligned', f'({obj1} or {obj2}) and not mutations and not non_mutations')
    
    # Create the view and coloring
    selection = f"{obj1} or {obj2}"
    cmd.hide('everything', selection)
    cmd.show('cartoon', selection)
    cmd.color('grey', 'not_aligned')
    for (col, sel) in colors:
        cmd.color(col, sel)

    cmd.hide('everything', '(hydro) and (%s or %s)'%(obj1, obj2))
    cmd.center('%s or %s'%(obj1, obj2))
    if labels:
        cmd.label('mutations and name CA','"(%s-%s-%s)"%(chain, resi, resn)')
    if waters:
        cmd.set('sphere_scale', '0.1')
        cmd.show('spheres', 'resn HOH and (%s or %s)'%(obj1, obj2))
        cmd.color('red', 'resn HOH and %s'%obj1)
        cmd.color('salmon', 'resn HOH and %s'%obj2)
    cmd.hide('cgo')
    cmd.deselect()

    print(f'''Mutations are shown in yellow-to-red colormap, based on BLOSUM90 similarity.
    Identical residues are shown in green whereas unaligned regions are represented in grey.''')
    
cmd.extend("color_by_mutation", color_by_mutation)

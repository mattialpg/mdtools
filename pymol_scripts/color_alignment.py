from pymol import cmd

# Colors for match types
match_color = "green"
conserved_color = "yellow"
nonconserved_color = "red"

# Define groups of similar amino acids for conservative substitutions
similar_groups = [
    set("STA"),   # small polar
    set("NQ"),    # amide
    set("DE"),    # acidic
    set("KR"),    # basic
    set("MILV"),  # hydrophobic
    set("FYWH")   # aromatic
]

def are_similar(res1, res2):
    if res1 == res2:
        return True
    for group in similar_groups:
        if res1 in group and res2 in group:
            return True
    return False

def color_aligned_residues(obj1, obj2, aligned_seq1, aligned_seq2):
    """
    Color residues in two objects based on aligned sequences.
    Assumes residue numbering corresponds to alignment position (excluding gaps).
    """

    # Track residue indices in each protein (exclude gaps)
    index1 = 0
    index2 = 0

    for i, (r1, r2) in enumerate(zip(aligned_seq1, aligned_seq2)):
        # Increase residue index counters only if residue is not gap
        if r1 != '-':
            index1 += 1
        if r2 != '-':
            index2 += 1

        # Skip alignment positions with gaps in either seq (no comparison)
        if r1 == '-' or r2 == '-':
            continue

        # Determine residue color based on conservation
        if r1 == r2:
            color = match_color
        elif are_similar(r1, r2):
            color = conserved_color
        else:
            color = nonconserved_color

        # Color residue in prot1 (by residue number index1)
        cmd.color(color, f"{obj1} and resi {index1}")
        # Color residue in prot2 (by residue number index2)
        cmd.color(color, f"{obj2} and resi {index2}")

# Run coloring on your objects and aligned sequences
color_aligned_residues(prot1, prot2, seq1, seq2)

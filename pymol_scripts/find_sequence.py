from pymol import cmd


def find_sequence(arg1, arg2=None, chain=None, sel_name="motif"):
    if arg2 is None:
        motif = str(arg1).strip().upper()
        obj_name = "sele" if "sele" in cmd.get_names("selections") else "enabled"
    else:
        obj_name = str(arg1).strip()
        motif = str(arg2).strip().upper()

    aa3_to_1 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }

    atom_sel = f"({obj_name}) and polymer.protein and name CA"
    if chain:
        atom_sel += f" and chain {chain}"

    model = cmd.get_model(atom_sel)

    by_chain = {}
    seen = set()
    for atom in model.atom:
        key = (atom.model, atom.chain, atom.resi)
        if key in seen:
            continue
        seen.add(key)
        one = aa3_to_1.get(atom.resn.upper())
        if not one:
            continue
        chain_key = (atom.model, atom.chain)
        by_chain.setdefault(chain_key, []).append((one, atom.resi))

    selections = []
    total_matches = 0
    for (model_id, chain_id), residues in by_chain.items():
        seq = "".join(one for one, _ in residues)
        for i in range(len(seq) - len(motif) + 1):
            if seq[i : i + len(motif)] != motif:
                continue
            total_matches += 1
            for j in range(len(motif)):
                resi_id = residues[i + j][1]
                if chain_id:
                    selections.append(
                        f"(model {model_id} and chain {chain_id} and resi {resi_id})"
                    )
                else:
                    selections.append(f"(model {model_id} and resi {resi_id})")

    if selections:
        selection_string = " or ".join(selections)
        cmd.select(sel_name, selection_string)
        print(f"Found {total_matches} match(es), selection: {sel_name}")
    else:
        print("No matches found")


cmd.extend("find_sequence", find_sequence)

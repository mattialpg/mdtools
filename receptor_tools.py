from pathlib import Path
from pymol2 import PyMOL
from pymol import stored
import subprocess
import requests
# from modeller import Environ, Model, Alignment, automodel


def extract_receptor(pdb_infile, pdb_outfile, chain='A'):
    """Extract a polymer chain from a PDB/CIF using PyMOL."""
    pdb_path = Path(pdb_infile).resolve()
    out = Path(pdb_outfile).resolve()

    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    with PyMOL() as pm:
        pm.cmd.load(str(pdb_path), 'prot')
        pm.cmd.select('sel', f"polymer and chain {chain}")

        # Check secondary structure 
        pm.cmd.dss('sel')
        stored.missing_ss = 0
        pm.cmd.iterate('sel and name CA', 'stored.missing_ss += 1 if ss in ("", "C") else 0')
        if stored.missing_ss:
            print(f"Warning: {stored.missing_ss} residues lack secondary structure in chain {chain}")

        pm.cmd.save(str(out), 'sel')


def reconstruct_loops(pdb_id, chain='A'):
    """Simplified MODELLER-based loop reconstruction."""
    env = Environ()
    model_code = f"{pdb_id}_{chain}"
    pdb_path = Path(f"{pdb_id}.pdb")
    aln = Alignment(env)
    mdl = Model(env, file=pdb_path.stem, model_segment=(f"FIRST:{chain}", f"LAST:{chain}"))
    aln.append_model(mdl, align_codes='template', atom_files=f"{pdb_path}")
    aln.append(file=f"{pdb_id}.ali", align_codes=model_code)
    aln.align2d()
    a = automodel.AutoModel(env, alnfile=f"template-{model_code}.ali",
                            knowns='template', sequence=model_code)
    a.starting_model = 1
    a.ending_model = 1
    a.make()
    return Path(f"bestmodel_{model_code}.pdb")


def fix_receptor(pdb_id, chain='A', run_loops=True):
    if run_loops:
        reconstruct_loops(pdb_id, chain)
    return extract_receptor(f"{pdb_id}.pdb", chain)

from pathlib import Path
from pymol2 import PyMOL
from itertools import pairwise
import subprocess

# import logging
# logger = logging.getLogger(__name__)


def extract_receptor(pdb_infile, pdb_outfile, chain='A'):
    """Extract a polymer chain from a PDB/CIF using PyMOL."""
    #TODO: adjust residue numbering offset in protein
    
    pdb_path = Path(pdb_infile).resolve()
    out = Path(pdb_outfile).resolve()

    if not pdb_path.exists():
        logger.error(f"PDB file not found: {pdb_path}")

    with PyMOL() as pm:
        pm.cmd.load(str(pdb_path), 'prot')
        pm.cmd.select('sel', f"polymer and chain {chain}")
        pm.cmd.save(str(out), 'sel')

        # Check for residue gaps
        residues = []
        pm.cmd.iterate('sel and name CA', 'residues.append(resv)',
            space={'residues': residues})
        gaps = [(a + 1, b - 1) for a, b in pairwise(sorted(set(residues))) if b - a > 1]
        if gaps:
            gap_text = ', '.join(f"{s}-{e}" if s != e else f"{s}" for s, e in gaps)
            print(f"\nMissing residues in chain {chain}: {gap_text}")
            return True


def reconstruct_loops(pdb_id, chain='A'):
    """Simplified MODELLER-based loop reconstruction."""
    import os
    from pathlib import Path
    from urllib.request import urlopen
    from Bio import SeqIO
    from modeller import Environ, Model, Alignment
    from modeller.automodel import automodel, assess

    pdb_path = Path(f"{pdb_id}.pdb")
    model_code = f"{pdb_id}.{chain}"
    model_folder = Path("protein.modeller")
    model_folder.mkdir(exist_ok=True)

    # Dowload full fasta sequence
    url = f"https://www.rcsb.org/fasta/chain/{pdb_id}.{chain}/download"
    full_fasta = model_folder / f"{pdb_id}.fasta"
    with urlopen(url) as r:
        fasta_text = r.read().decode('utf-8')
    Path(full_fasta).write_text(fasta_text)

    # Convert to PIR format (.ali)
    rec = next(SeqIO.parse(full_fasta, 'fasta'))
    seq = str(rec.seq).replace('*', '').replace('X', 'X')
    ali_path = model_folder / f"{pdb_id}.ali"
    ali_text = (
        f">P1;{model_code}\n"
        f"sequence:{model_code}:::::::0.00:0.00\n"
        f"{seq}*\n")
    ali_path.write_text(ali_text)
    subprocess.run(['cp', str(pdb_path), str(model_folder)], check=True)

    # Load modeller
    env = Environ()
    aln = Alignment(env)
    mdl = Model(env, file=pdb_path.stem, model_segment=(f"FIRST:{chain}", f"LAST:{chain}"))
    aln.append_model(mdl, align_codes='template', atom_files=str(pdb_path))
    aln.append(file=str(ali_path), align_codes=model_code)

    # Align template to target sequence
    aln.align2d()
    alnfile = model_folder / f"template-{model_code}.ali"
    aln.write(file=str(alnfile), alignment_format='PIR')

    # Reconstruct loops
    cwd = os.getcwd()
    os.chdir(model_folder)
    a = automodel(env, alnfile=alnfile.name, knowns='template',
        sequence=model_code, assess_methods=(assess.DOPE,))
    a.starting_model = 1
    a.ending_model = 1
    a.make()
    os.chdir(cwd)

    # Select best model
    best = min(a.outputs, key=lambda x: x['DOPE score'])
    best_model = best['name']
    subprocess.run(['cp', str(model_folder / best_model), "protein.pdb"], check=True)
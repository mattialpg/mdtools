from pathlib import Path
from pymol2 import PyMOL
from itertools import pairwise

import logging
logger = logging.getLogger(__name__)


def extract_receptor(pdb_infile, pdb_outfile, chain='A'):
    """Extract a polymer chain from a PDB/CIF using PyMOL."""
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
            logger.warning(f"Missing residues in chain {chain}: {gap_text}")
            return True


def reconstruct_loops(pdb_id, chain='A'):
    """Simplified MODELLER-based loop reconstruction."""
    from modeller import Environ, Model, Alignment, automodel

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

    logger.info(Path(f"bestmodel_{model_code}.pdb"))
from pathlib import Path
from itertools import pairwise
import os
import re
import shutil
import subprocess
from urllib.request import urlopen

from Bio import SeqIO
from modeller import Environ, Model, Alignment
from modeller.automodel import LoopModel, assess, refine
from pymol2 import PyMOL


ARCHIVE_DIR = Path("/home/mattia/modeller_archive")
MODEL_DIR_NAME = "protein.modeller"


def _copy_archived_model(archived_model_dir, model_dir, pdb_outfile):
    if not archived_model_dir.is_dir():
        return False
    if model_dir.exists():
        shutil.rmtree(model_dir)
    shutil.copytree(archived_model_dir, model_dir)
    print(f"\nCopied archived model to {model_dir}\n")

    archived_protein = model_dir / "protein.pdb"
    if archived_protein.exists():
        shutil.copy2(archived_protein, pdb_outfile)
        print(f"Using archived reconstructed receptor: {archived_protein}")
    return True


def _archive_new_model(model_dir, archived_model_dir):
    if not model_dir.exists():
        return
    archived_model_dir.parent.mkdir(parents=True, exist_ok=True)
    if archived_model_dir.exists():
        shutil.rmtree(archived_model_dir)
    shutil.copytree(model_dir, archived_model_dir)
    print(f"Archived reconstructed model to: {archived_model_dir}")


def extract_receptor(pdb_infile, pdb_outfile, chain="A", domain=None):
    """Extract polymer selection and report whether sequence gaps are present."""
    #TODO: check enginerised residues

    pdb_path = Path(pdb_infile).resolve()
    out = Path(pdb_outfile).resolve()
    if not pdb_path.exists():
        print(f"PDB file not found: {pdb_path}")

    with PyMOL() as pm:
        pm.cmd.load(str(pdb_path), "prot")
        metal_resn = "+".join(["MG", "ZN", "MN", "FE", "CU", "CO", "NI", "CA"])
        sel_string = f"(polymer and chain {chain}) or (inorganic and chain {chain} and resn {metal_resn})"
        if domain:
            sel_string = (
                f"((polymer and chain {chain} and resi {domain}) or "
                f"(inorganic and chain {chain} and resn {metal_resn}))"
            )
        pm.cmd.select("sel", sel_string)
        pm.cmd.save(str(out), "sel")

        residues = []
        pm.cmd.iterate("sel and name CA", "residues.append(resv)", space={"residues": residues})
        gaps = [(a + 1, b - 1) for a, b in pairwise(sorted(set(residues))) if b - a > 1]
        if gaps:
            gap_text = ", ".join(f"{s}-{e}" if s != e else f"{s}" for s, e in gaps)
            print(f"\nMissing residues in chain {chain}: {gap_text}")
            return True
    return False


def reconstruct_loops(pdb_id, chain="A", domain=None):
    """Run MODELLER loop reconstruction and write best model as protein.pdb in CWD."""
    pdb_path = Path(f"{pdb_id}.pdb")
    model_dir = Path(MODEL_DIR_NAME)
    model_dir.mkdir(exist_ok=True)

    full_fasta = model_dir / f"{pdb_id}.fasta"
    with urlopen(f"https://www.rcsb.org/fasta/chain/{pdb_id}.{chain}/download") as response:
        full_fasta.write_text(response.read().decode("utf-8"))

    start = end = None
    if domain:
        match = re.fullmatch(R"\s*(\d+)\s*-\s*(\d+)\s*", str(domain))
        start, end = map(int, match.groups())
    rec = next(SeqIO.parse(full_fasta, "fasta"))
    seq = str(rec.seq).replace("*", "")
    if start is not None:
        seq = seq[start - 1 : end]

    ali_path = model_dir / f"{pdb_id}.ali"
    ali_path.write_text(">P1;target\n"
        "sequence:target:::::::0.00:0.00\n"
        f"{seq}*\n")

    shutil.copy(str(pdb_path), str(model_dir / pdb_path.name))

    env = Environ()
    aln = Alignment(env)
    if start is not None:
        model_segment = (f"{start}:{chain}", f"{end}:{chain}")
    else:
        model_segment = (f"FIRST:{chain}", f"LAST:{chain}")
    mdl = Model(env, file=pdb_path.stem, model_segment=model_segment)
    aln.append_model(mdl, align_codes="template", atom_files=str(pdb_path))
    aln.append(file=str(ali_path), align_codes="target")
    aln.align2d()

    alnfile = model_dir / "template-target.ali"
    aln.write(file=str(alnfile), alignment_format="PIR")

    cwd = os.getcwd()
    os.chdir(model_dir)
    try:
        modeler = LoopModel(env, alnfile=alnfile.name,
            knowns="template", sequence="target",
            assess_methods=(assess.DOPE,))
        modeler.starting_model = 1
        modeler.ending_model = 1
        modeler.loop.starting_model = 1
        modeler.loop.ending_model = 1
        modeler.loop.md_level = refine.slow
        modeler.make()
    finally:
        os.chdir(cwd)

    valid_models = [item for item in modeler.outputs if item.get("failure") is None]
    if not valid_models:
        raise RuntimeError("MODELLER did not produce any successful models.")
    best = min(valid_models, key=lambda item: item["DOPE score"])
    best_model = model_dir / best["name"]

    # Re-align best model to original structure
    with PyMOL() as pm:
        pm.cmd.load("protein.pdb", "ref")
        pm.cmd.load(str(best_model), "best")
        pm.cmd.align("best", "ref")
        pm.cmd.create("merged", "best or (ref and inorganic)")
        pm.cmd.save("protein.pdb", "merged")


def prepare_receptor(configs, archive_dir=ARCHIVE_DIR):
    """Prepare receptor PDB/PDBQT using archive-or-reconstruct logic."""
    workdir = Path(configs["workdir"])
    protein_cfg = configs.get("protein", {}) if isinstance(configs.get("protein"), dict) else {}
    receptor_id = configs.get("receptor_id") or protein_cfg.get("id")
    receptor_name = configs.get("receptor_name") or protein_cfg.get("name") or "protein"
    print(f"Preparing receptor for {receptor_id}")

    parts = receptor_id.split(":")
    pdb_id, chain = parts[0], parts[1]
    domain = parts[2] if len(parts) > 2 else None

    pdb_infile = workdir / f"{pdb_id}.pdb"
    pdb_outfile = workdir / f"{receptor_name}.pdb"
    receptor_pdbqt = pdb_outfile.with_suffix(".pdbqt")

    has_missing_loops = extract_receptor(pdb_infile, pdb_outfile, chain, domain)
    if has_missing_loops:
        modeller_dir_name = f"{receptor_id.replace(':', '_')}.modeller"
        archived_model_dir = Path(archive_dir) / modeller_dir_name
        model_dir = workdir / MODEL_DIR_NAME

        used_archive = _copy_archived_model(archived_model_dir, model_dir, pdb_outfile)
        if not used_archive:
            answer = input("Do you want to reconstruct loops? [Y/n]: ").strip().lower()
            if answer in ("", "y", "yes"):
                reconstruct_loops(pdb_id, chain, domain)
                _archive_new_model(model_dir, archived_model_dir)

    subprocess.run([configs["obabel"], str(pdb_outfile), "-xr", "-O", str(receptor_pdbqt)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return receptor_pdbqt

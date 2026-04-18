import sys
from pathlib import Path
import subprocess
import json
from pymol2 import PyMOL
from rdkit import Chem
from rdkit.Chem import AllChem

def check_ionisation(ligand_smiles, ph=7.0, report_file="unipka_report.txt"):
    script_path = Path(__file__).resolve().parent.parent / "docking_pipeline" / "unipka_ops.py"
    try:
        proc = subprocess.run([str(script_path), ligand_smiles, str(ph), str(report_file)], 
            check=True, capture_output=True, text=True)
        payload = json.loads(proc.stdout.strip() or "{}")
        chosen = payload.get("smiles", ligand_smiles)
    except Exception as exc:
        print(f"WARNING: Uni-pKa ionisation failed, using input SMILES ({exc})")
        return ligand_smiles
    return chosen


def prepare_ligand(configs):
    """Prepare ligand from SMILES and return the ligand PDBQT path."""
    workdir = Path(configs["workdir"])
    ligand_cfg = {}
    ligands = configs.get("ligands")
    if isinstance(ligands, list) and ligands and isinstance(ligands[0], dict):
        ligand_cfg = ligands[0]
    elif isinstance(ligands, dict):
        ligand_cfg = ligands

    ligand_id = configs.get("ligand_id") or ligand_cfg.get("id")
    ligand_name = configs.get("ligand_name") or ligand_cfg.get("name") or "ligand"
    ligand_smiles = configs.get("ligand_smiles") or ligand_cfg.get("smiles")
    obabel = configs["obabel"]

    print(f"Preparing ligand {ligand_id}")
    sdf_free = workdir / f"{ligand_name}_free.sdf"
    pdbqt_free = workdir / f"{ligand_name}_free.pdbqt"

    subprocess.run([obabel, f"-:{ligand_smiles}", "-O", str(sdf_free), "--gen3d"],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run([obabel, str(sdf_free), "-O", str(pdbqt_free)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    print(f"Ligand prepared: {pdbqt_free}")
    return pdbqt_free


def extract_ligand(ligand_id, pdb_id, sdf_outfile):
    """Extract ligand from receptor PDB and save as SDF."""
    pdb_path = Path(f"{pdb_id}.pdb").resolve()
    out = Path(sdf_outfile).resolve()
    if not pdb_path.exists():
        print(f"PDB file not found: {pdb_path}")

    lig_resname, lig_chain, lig_resid = ligand_id.split(":")

    with PyMOL() as pm:
        pm.cmd.load(str(pdb_path), "prot")
        sel_string = f"r. {lig_resname} and chain {lig_chain} and resi {lig_resid}"
        pm.cmd.select("sel", sel_string)
        pm.cmd.create("ligand_obj", "sel")
        pm.cmd.h_add("ligand_obj")
        pm.cmd.save(str(out), "ligand_obj", format="sdf")

    subprocess.run(["sed", "-i", "1s/.*//", str(out)], check=True)
    subprocess.run(["sed", "-i", f"2s/.*/{lig_resname}/", str(out)], check=True)
    subprocess.run(["sed", "-i", "3s/.*//", str(out)], check=True)

from pathlib import Path
import json
import re

from pymol2 import PyMOL

import sys
sys.path.append('/home/mattia/aiddtools')
import interaction_tools as inttools
from .tools import run_in_terminal, run_capture


def check_ionisation(ligand_smiles, ph=7.0, report_file="unipka_report.txt"):
    script_path = Path(__file__).resolve().with_name("unipka_ops.py")
    try:
        proc = run_capture([str(script_path), ligand_smiles, str(ph), str(report_file)])
        payload = json.loads(proc.stdout.strip() or "{}")
        chosen = payload.get("smiles", ligand_smiles)
    except Exception as exc:
        print(f"WARNING: Uni-pKa ionisation failed, using input SMILES ({exc})")
        return ligand_smiles
    return chosen


def prepare_ligand(ctx):
    """Prepare ligand from SMILES and return the ligand PDBQT path."""
    workdir = Path(ctx.workdir)
    ligand_id = ctx.query_ligand["id"]
    ligand_name = ctx.ligand_name
    ligand_smiles = ctx.query_ligand["smiles"]
    obabel = ctx.obabel

    print(f"Preparing ligand {ligand_id}")
    sdf_free = workdir / f"{ligand_name}_free.sdf"
    pdbqt_free = workdir / f"{ligand_name}_free.pdbqt"

    run_in_terminal([obabel, f"-:{ligand_smiles}", "-O", str(sdf_free), "--gen3d"])
    run_in_terminal([obabel, str(sdf_free), "-O", str(pdbqt_free)])

    print(f"Ligand prepared: {pdbqt_free}")
    ctx.ligand_pdbqt = pdbqt_free
    return ctx


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

    run_in_terminal(["sed", "-i", "1s/.*//", str(out)])
    run_in_terminal(["sed", "-i", f"2s/.*/{lig_resname}/", str(out)])
    run_in_terminal(["sed", "-i", "3s/.*//", str(out)])


def select_cognate_ligand(ctx):
    ligands = inttools.get_ligands(f"{ctx.pdb_id}.pdb")
    if not ligands:
        print(f"No ligands detected in {ctx.pdb_id}.pdb")
        return ctx
    if ctx.chain is not None:
        ligands = {k: v for k, v in ligands.items() if k.split(':')[1] == ctx.chain}

    sorted_ligands = sorted(ligands.keys(), key=lambda key: (key.split(':')[1], 
        int(re.match(R"^\s*(-?\d+)", key.split(':')[2]).group(1)), key.split(':')[0],),)

    print("Detected ligands:")
    for i, key in enumerate(sorted_ligands, start=1):
        loi_status, _ = ligands[key]
        mark = f"  ({loi_status})" if loi_status == 'LOI' else ''
        print(f"   ({i}) Ligand {key}{mark}")

    choice = int(input("\nSelect a ligand to replace: ").strip())

    keys = list(sorted_ligands)
    ctx.cognate_lig_id, ctx.chain, _ = keys[choice - 1].split(':')

    ctx.receptor_id = f"{ctx.pdb_id}:{ctx.chain}"
    if ctx.domain is not None:
        ctx.receptor_id += f":{ctx.domain}"

    if len(keys) > 1:
        prompt = ("Select ligand(s) to keep or press Enter to continue: ")
        selected = input(prompt).strip()
        if selected:
            ctx.ligand_name = "ligand1"
            ctx.query_ligand["name"] = ctx.ligand_name
            for i, idx in enumerate(map(int, selected.replace(',', ' ').split()), start=2):
                key = keys[idx - 1]
                _, smiles = ligands[key]
                ctx.ligands.append({
                    "id": key,
                    "smiles": smiles,
                    "name": f"ligand{i}",
                    "resname": key.split(':')[0],
                    "role": "kept",})
                extract_ligand(key, ctx.pdb_id, f"ligand{i}.sdf")
    return ctx

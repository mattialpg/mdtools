#!/usr/bin/env python3
import argparse
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import quote
import requests
import tools.utils as utils


def cmd_download_pdb(args):
    outpath = utils.download_pdb(args.pdb_id, outdir=args.outdir)
    print(f"Downloaded {Path(outpath).name}")
    return 0

def _resolve_name_to_smiles_pubchem(name: str) -> str:
    encoded = quote(name.strip(), safe="")
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{encoded}/property/CanonicalSMILES/TXT")
    response = requests.get(url, timeout=20)
    if response.status_code != 200:
        raise RuntimeError(
            f"Could not resolve ligand name '{name}' on PubChem (HTTP {response.status_code}).")

    smiles = response.text.strip().splitlines()[0].strip() if response.text.strip() else ""
    if not smiles:
        raise RuntimeError(f"PubChem returned no SMILES for ligand name '{name}'.")
    return smiles


def _normalize_ligand_to_smiles(ligand_value: str) -> str:
    text = ligand_value.strip()
    if not text:
        raise RuntimeError("Ligand value is empty.")

    try:
        from rdkit import Chem
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.error")
        if Chem.MolFromSmiles(text) is not None:
            return text
    except Exception:
        # Fallback heuristic if RDKit is unavailable in the runtime.
        if re.fullmatch(R"[A-Za-z0-9@+\-\[\]\(\)=#$\\/%.]+", text):
            return text

    resolved = _resolve_name_to_smiles_pubchem(text)
    print(f"Resolved ligand '{text}' -> {resolved}")
    return resolved


def cmd_docking_session(args):
    cmd = [sys.executable, str(Path(__file__).with_name("docking_session.py"))]
    if args.pdb:
        cmd.extend(["--pdb", args.pdb])
    if args.lig_new:
        resolved_smiles = _normalize_ligand_to_smiles(args.lig_new)
        cmd.extend(["--lig_new", resolved_smiles])
    if args.lig_file:
        cmd.extend(["--lig_file", args.lig_file])
    if args.verbose:
        cmd.append("--verbose")
    return subprocess.run(cmd).returncode


def cmd_interaction_analysis(args):
    cmd = [sys.executable,
        str(Path(__file__).with_name("interaction_analysis.py")),
        *args.extra]
    return subprocess.run(cmd).returncode


def build_parser():
    parser = argparse.ArgumentParser(
        description="mdtools command line interface for agent and pipeline usage")
    sub = parser.add_subparsers(dest="command", required=True)

    p_download = sub.add_parser("download-pdb", help="Download a PDB by ID")
    p_download.add_argument("--pdb-id", required=True, dest="pdb_id")
    p_download.add_argument("--outdir", default=".")
    p_download.set_defaults(func=cmd_download_pdb)

    p_docking = sub.add_parser(
        "docking-session", help="Run mdtools docking_session.py workflow")
    p_docking.add_argument("--pdb", help="Receptor ID (e.g. 1ABC or 1ABC:A:45-210)")
    p_docking.add_argument("--lig-new", dest="lig_new", help="New ligand SMILES")
    p_docking.add_argument("--lig-file", dest="lig_file", help="CSV file of ligands")
    p_docking.add_argument("--verbose", action="store_true")
    p_docking.set_defaults(func=cmd_docking_session)

    p_int = sub.add_parser("interaction-analysis",
        help="Forward arguments to interaction_analysis.py")
    p_int.add_argument("extra", nargs=argparse.REMAINDER)
    p_int.set_defaults(func=cmd_interaction_analysis)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    code = args.func(args)
    raise SystemExit(code)


if __name__ == "__main__":
    main()

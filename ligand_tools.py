from pathlib import Path
import subprocess


def prepare_ligand(configs):
    """Prepare ligand from SMILES and return the ligand PDBQT path."""
    workdir = Path(configs["workdir"])
    ligand_cfg = configs.get("ligand", {}) if isinstance(configs.get("ligand"), dict) else {}
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

from pathlib import Path
import subprocess


def prepare_ligand(configs):
    """Prepare ligand from SMILES and return the ligand PDBQT path."""
    workdir = Path(configs["workdir"])
    ligand_id = configs["ligand_id"]
    ligand_name = configs["ligand_name"]
    ligand_smiles = configs["ligand_smiles"]
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

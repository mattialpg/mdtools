from dataclasses import dataclass, field
from pathlib import Path

TOOL_PATHS = {
    "pymol": "/home/mattia/.miniforge/envs/my-chem/bin/pymol",
    "obabel": "/home/mattia/.miniforge/envs/my-chem/bin/obabel",
    "vina": "/home/mattia/.bin/smina.static"}


@dataclass
class RunContext:
    receptor_id: str
    workdir: Path

    receptor_name: str = "protein"
    ligand_name: str = "ligand"

    pdb_id: str = ""
    chain: str | None = None
    domain: str | None = None
    cognate_lig_id: str | None = None

    ligands: list[dict] = field(default_factory=list)

    pymol: str = ""
    obabel: str = ""
    vina: str = ""

    box_center: tuple | None = None
    box_size: tuple | None = None

    receptor_pdbqt: Path | None = None
    ligand_pdbqt: Path | None = None
    docked_name: str | None = None
    docked_file: Path | None = None

    energy_range: int = 4
    exhaustiveness: int = 16
    num_modes: int = 8

    @classmethod
    def from_job(cls, job):
        receptor_id, ligand_id, ligand_smiles = job
        parts = receptor_id.split(":")
        pdb_id = parts[0]
        chain = parts[1] if len(parts) > 1 else None
        domain = parts[2] if len(parts) > 2 else None
        workdir = Path(f"{pdb_id}_{ligand_id}")

        ligands = [{
            "id": ligand_id,
            "smiles": ligand_smiles,
            "name": "ligand",
            "resname": "LIG",
            "role": "docked",}]

        return cls(
            receptor_id=receptor_id,
            workdir=workdir,
            pdb_id=pdb_id,
            chain=chain,
            domain=domain,
            ligands=ligands,
            pymol=TOOL_PATHS["pymol"],
            obabel=TOOL_PATHS["obabel"],
            vina=TOOL_PATHS["vina"],)

    @property
    def query_ligand(self):
        return self.ligands[0]

    @property
    def ligand_smiles(self):
        return self.ligands[0]["smiles"]

    @ligand_smiles.setter
    def ligand_smiles(self, value):
        self.ligands[0]["smiles"] = value

    @property
    def ligand_resname(self):
        return self.ligands[0]["resname"]

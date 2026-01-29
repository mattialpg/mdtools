import os, sys, re
import subprocess, shutil
from pathlib import Path
import logging
import pandas as pd

import yaml
def tuple_as_literal(dumper, data):
    """ Serialises tuples as plain strings """
    return dumper.represent_scalar("tag:yaml.org,2002:str", str(data))
yaml.SafeDumper.add_representer(tuple, tuple_as_literal)

# Import AIDDTools utils
sys.path.append("/home/mattia/aiddtools")
import interaction_tools as inttools
import receptor_tools
import utils


class Preprocess:
    def __init__(self, args):
        self.logger = utils.get_logger()

        self.receptor_id = args.pdb
        self.lig_new_smiles = args.lig_new
        self.lig_new_id = "LIG"
        self.lig_new_md_id = "LIG"
        
        self.pymol = '/home/mattia/.miniforge3/envs/fraglib/bin/pymol'
        self.obabel = '/home/mattia/.miniforge3/envs/fraglib/bin/obabel'
        self.vina = '/home/mattia/.bin/smina.static'


    def choose_ligand(self):
        ligands = inttools.get_ligands(f"{self.receptor_id}.pdb")

        if not ligands:
            self.logger.warning(f"No ligands detected in {self.receptor_id}.pdb")
            return

        self.logger.info("Detected ligands:")
        for i, (lig_id, loi_status) in enumerate(ligands.items(), start=1):
            mark = f"  ({loi_status})" if loi_status == "Non-LOI" else ""
            print(f"   ({i}) Ligand {lig_id}{mark}")

        choice = int(input("\nSelect a ligand to replace: ").strip())

        keys = list(ligands.keys())
        selected_key = keys[choice - 1]
        resname, chain, _ = selected_key.split(":")

        self.lig_orig_id = resname
        self.chain = chain


    def get_box(self):
        pymol_commands = (
            f"run /home/mattia/mdtools/pymol_scripts/draw_bounding_box.py;"
            f"load {self.receptor_id}.pdb, prot;"
            f"select lig, chain {self.chain} and resn {self.lig_orig_id};"
            f"draw_bounding_box lig")
        result = subprocess.run([self.pymol, "-cqd", pymol_commands],
            check=True, capture_output=True, text=True).stdout.strip()

        lines = result.splitlines()[-2:]
        for line in lines:
            if line.startswith("Box center"):
                box_center = re.sub(R"^Box center\s*", "", line).strip()
            elif line.startswith("Box dimensions"):
                box_size = re.sub(R"^Box dimensions\s*", "", line).strip()

        self.box_center = tuple(float(x) for x in eval(box_center))
        self.box_size = tuple(float(x) for x in eval(box_size))


    def write_conf_file(self, workdir):
        conf_file = workdir / "config.yaml"

        configs = {
            "receptor_id": self.receptor_id,
            "receptor_name": "protein",
            "chain": self.chain,

            "ligand_id": self.lig_new_id,
            "ligand_smiles": self.lig_new_smiles,
            "ligand_name": "ligand",
            "ligand_md_id": self.lig_new_md_id,

            "box_center": self.box_center,
            "box_size": self.box_size,

            "energy_range": 4,
            "exhaustiveness": 16,
            "num_modes": 8,

            "workdir": str(workdir),
            "pymol": self.pymol,
            "obabel": self.obabel,
            "vina": self.vina,}

        # Export to config file
        config_text = yaml.safe_dump(configs, sort_keys=False)
        for key in ("ligand_id:", "box_center:", "energy_range:", "workdir:"):
            config_text = config_text.replace(f"\n{key}", f"\n\n{key}")
        conf_file.write_text(config_text)

        self.logger.info(f"Configuration file written: {conf_file}")
        return configs
    

class DockingSession:
    def __init__(self, configs):
        self.logger = utils.get_logger()

        # Bind configs as attributes
        for key, value in configs.items():
            setattr(self, key, value)
        self.workdir = Path(self.workdir)
        self.box = {"center": tuple(self.box_center),
            "size": tuple(self.box_size)}

        # Runtime-only attributes
        self.receptor_pdbqt = None
        self.ligand_pdbqt = None


    def prepare_receptor(self, fix_loops=False):
        self.logger.info(f"Preparing receptor for {self.receptor_id} (chain {self.chain})")

        pdb_infile = self.workdir / f"{self.receptor_id}.pdb"
        pdb_outfile = self.workdir / f"{self.receptor_name}.pdb"

        if fix_loops:
            receptor_tools.fix_receptor(self.workdir / self.receptor_id, self.chain)
        
        receptor_tools.extract_receptor(pdb_infile, pdb_outfile, self.chain)
        self.receptor_pdbqt = pdb_outfile.with_suffix(".pdbqt")
        subprocess.run([self.obabel, str(pdb_outfile), "-xr", "-O", str(self.receptor_pdbqt)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


    def prepare_ligand(self):
        self.logger.info(f"Preparing ligand {self.ligand_id}")
        mol2_free = self.workdir / f"{self.ligand_name}_free.mol2"
        pdbqt_free = self.workdir / f"{self.ligand_name}_free.pdbqt"

        # Generate 3D ligand from SMILES
        subprocess.run([self.obabel, f"-:{self.ligand_smiles}", "-O", str(mol2_free), "--gen3d"],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Convert MOL2 to PDBQT
        subprocess.run([self.obabel, str(mol2_free), "-O", str(pdbqt_free)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.logger.info(f"Ligand prepared: {pdbqt_free}")
        self.ligand_pdbqt = pdbqt_free


    def dock(self, output_prefix="docked"):
        self.logger.info(f"Running docking for {self.receptor_id}")

        self.docked_name = f"{self.ligand_name}_{output_prefix}"
        self.docked_file = self.workdir / f"{self.docked_name}.pdbqt"
        vina_conf = self.workdir / f"{self.docked_name}_vina.conf"
        log_file = self.workdir / f"{self.docked_name}.log"

        cx, cy, cz = self.box_center
        sx, sy, sz = self.box_size

        # Write vina configuration file
        vina_conf.write_text(f"""
            receptor = {self.receptor_pdbqt}
            ligand = {self.ligand_pdbqt}
            center_x = {cx}
            center_y = {cy}
            center_z = {cz}
            size_x = {sx}
            size_y = {sy}
            size_z = {sz}
            exhaustiveness = {self.exhaustiveness}
            energy_range = {self.energy_range}
            num_modes = {self.num_modes}""")

        # Run vina via subprocess
        subprocess.run([self.vina, "--config", str(vina_conf), "--out", str(self.docked_file),
            "--log", str(log_file)], check=True)
        vina_conf.unlink(missing_ok=True)

        self.logger.info(f"Docking completed. Results written to {self.docked_file}")


    def postprocess(self):
        # Split docked poses into individual MOL2 files
        subprocess.run([self.obabel, str(self.docked_file), "-O",
            str(self.workdir / f"{self.docked_name}_0*.mol2"), "-h"],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Add hydrogens using PyMOL
        pymol_script = self.workdir / "add_hydrogens.pml"
        with pymol_script.open('w') as pml:
            for f in self.workdir.glob(f"{self.docked_name}_0*.mol2"):
                pml.write(f"load {f}, lig; h_add lig; save {f}, lig; delete lig;\n")
            pml.write("quit\n")
        subprocess.run([self.pymol, "-cq", "-r", str(pymol_script)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        pymol_script.unlink()

        # Fix MOL2 files
        perl_script = Path(__file__).parent / "sort_mol2_bonds.pl"
        for f in self.workdir.glob(f"{self.docked_name}_0*.mol2"):
            # Remove non-standard header lines
            subprocess.run(["sed", "-i", "1{/^#/d;}", str(f)], check=True)

            # Sort bonds with Perl script
            subprocess.run(["perl", str(perl_script), str(f), str(f)],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            # Fix ligand name
            pattern = Rf"s/\b\(UNL\|{self.ligand_md_id}\)[[:space:]]*[0-9]*/{self.ligand_md_id}/g"
            subprocess.run(["sed", "-i", pattern, str(f)], check=True)

        result_dir = self.workdir / f"{self.ligand_name}.vina"
        if result_dir.exists():
            shutil.rmtree(result_dir)
        result_dir.mkdir()

        for f in self.workdir.glob('*'):
            if f.is_file() and (f.suffix == ".pdbqt" or f.name.startswith(self.ligand_name)):
                shutil.move(str(f), result_dir / f.name)

        self.logger.info(f"Docked files moved to folder: {result_dir}")


class RedockSession(DockingSession):
    """Redocking workflow."""
    def run(self):
        self.prepare_receptor()
        self.prepare_ligand()
        self.dock()
        self.postprocess()
        return True


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run a redocking workflow")
    parser.add_argument("--pdb", required=True, help="Protein PDB ID (e.g. 1ABC)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--lig_new", help="New ligand SMILES")
    group.add_argument("--lig_file", help="CSV file containing ligand SMILES")
    args = parser.parse_args()

    prep = Preprocess(args)
    prep.choose_ligand()
    prep.get_box()

    if args.lig_new:
        ligands_new = [("LIG", args.lig_new)]
    else:
        df = pd.read_csv(args.lig_file)
        ligands_new = [(row["ID"], row["ISOSMILES"]) for _, row in df.iterrows()]
    
    for lig_new_id, lig_new_smiles in ligands_new:
        prep.lig_new_id = lig_new_id
        prep.lig_new_smiles = lig_new_smiles

        workdir = Path(f"{args.pdb}_{lig_new_id}").resolve()
        pdb_src = Path(f"{args.pdb}.pdb")
        workdir.mkdir(exist_ok=True)
        shutil.copy2(pdb_src, workdir / pdb_src.name)

        configs = prep.write_conf_file(workdir)

        redock = RedockSession(configs)
        redock.run()

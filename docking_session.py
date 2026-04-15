import os, sys, re, logging
import subprocess, shutil
from pathlib import Path
from ast import literal_eval
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

import yaml
def tuple_as_literal(dumper, data):
    """Serialize tuples as plain '(x, y, z)' strings in YAML."""
    return dumper.represent_scalar('tag:yaml.org,2002:str', str(data))
yaml.SafeDumper.add_representer(tuple, tuple_as_literal)

import receptor_tools
import ligand_tools
import utils

# Import AIDDTools utils
sys.path.append('/home/mattia/aiddtools')
import interaction_tools as inttools


class Preprocess:
    def __init__(self, receptor_id, lig_new_id, lig_new_smiles):
        self.logger = logging.getLogger(__name__)

        self.receptor_name = 'protein'
        self.receptor_id = receptor_id
        parts = receptor_id.split(':')
        self.pdb_id = parts[0]
        self.chain = parts[1] if len(parts) > 1 else None
        self.domain = parts[2] if len(parts) > 2 else None

        self.lig_new_id = lig_new_id
        self.lig_new_smiles = lig_new_smiles
        self.ligand_name = 'ligand'
        self.lig_new_md_id = 'LIG'

        self.ligands = [{"id": self.lig_new_id, "smiles": self.lig_new_smiles,
            "name": self.ligand_name, "md_id": self.lig_new_md_id, "role": "docked"}]
        
        tool_paths = utils.get_tool_paths()
        self.pymol = tool_paths['pymol']
        self.obabel = tool_paths['obabel']
        self.vina = tool_paths['vina']


    def choose_ligands(self):
        ligands = inttools.get_ligands(f"{self.pdb_id}.pdb")
        if not ligands:
            self.logger.warning(f"No ligands detected in {self.pdb_id}.pdb")
            return
        if self.chain is not None:
            ligands = {k:v for k, v in ligands.items() if k.split(':')[1] == self.chain}

        sorted_ligands = sorted(
            ligands.keys(), key=lambda key: (key.split(':')[1],
                int(re.match(R"^\s*(-?\d+)", key.split(':')[2]).group(1)),
                key.split(':')[0]))

        print("Detected ligands:")
        for i, key in enumerate(sorted_ligands, start=1):
            loi_status, _ = ligands[key]
            mark = f"  ({loi_status})" if loi_status == 'LOI' else ''
            print(f"   ({i}) Ligand {key}{mark}")

        choice = int(input("\nSelect a ligand to replace: ").strip())

        keys = list(sorted_ligands)
        self.lig_orig_id, self.chain, _ = keys[choice - 1].split(':')

        self.receptor_id = f"{self.pdb_id}:{self.chain}"
        if self.domain is not None:
            self.receptor_id += f":{self.domain}"

        # Optionally select additional ligands to keep
        if len(keys) > 1:
            prompt = ("Select ligand(s) to keep or press Enter to continue: ")
            selected = input(prompt).strip()
            if selected:
                self.ligand_name = "ligand1"
                self.ligands[0]["name"] = self.ligand_name
                for i, idx in enumerate(map(int, selected.replace(',', ' ').split()), start=2):
                    key = keys[idx - 1]
                    _, smiles = ligands[key]
                    self.ligands.append({"id": key, "smiles": smiles,
                        "name": f"ligand{i}", "md_id": key.split(':')[0], "role": "kept"})
                    ligand_tools.extract_ligand(key, self.pdb_id, f"ligand{i}.sdf")


    def get_box(self):
        pymol_commands = (
            f"run /home/mattia/mdtools/pymol_scripts/draw_bounding_box.py;"
            f"load {self.pdb_id}.pdb, prot;"
            f"select lig, chain {self.chain} and resn {self.lig_orig_id};"
            f"draw_bounding_box lig")
        result = subprocess.run([self.pymol, '-cqd', pymol_commands],
            check=True, capture_output=True, text=True).stdout.strip()

        lines = result.splitlines()[-2:]
        for line in lines:
            if line.startswith('Box center'):
                box_center = re.sub(R"^Box center\s*", '', line).strip()
            elif line.startswith('Box dimensions'):
                box_size = re.sub(R"^Box dimensions\s*", '', line).strip()

        self.box_center = tuple(float(x) for x in literal_eval(box_center))
        self.box_size = tuple(float(x) for x in literal_eval(box_size))
    

class DockingPipeline:
    def __init__(self, config_file=Path("config.yaml")):
        self.logger = logging.getLogger(__name__)

        self.configs = utils.read_config_file(config_file)
        for key, value in self.configs.items():
            setattr(self, key, value)
        self.workdir = Path(self.workdir)

    def prepare_receptor(self):
        self.receptor_pdbqt = receptor_tools.prepare_receptor(self.configs)

    def prepare_ligand(self):
        self.ligand_pdbqt = ligand_tools.prepare_ligand(self.configs)


    def dock(self, output_prefix='docked'):
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
        subprocess.run([self.vina, '--config', str(vina_conf), '--out', str(self.docked_file),
            '--log', str(log_file)], check=True)
        vina_conf.unlink(missing_ok=True)

        self.logger.info(f"Docking completed. Results written to {self.docked_file}")


    def postprocess(self):
        # Split docked poses into individual SDF files
        subprocess.run([self.obabel, str(self.docked_file), '-O',
            str(self.workdir / f"{self.docked_name}_0*.sdf")],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Add hydrogens using PyMOL
        pymol_script = self.workdir / 'add_hydrogens.pml'
        with pymol_script.open('w') as pml:
            for f in self.workdir.glob(f"{self.docked_name}_0*.sdf"):
                pml.write(f"load {f}, lig; h_add lig; save {f}, lig; delete lig;\n")
            pml.write('quit\n')
        subprocess.run([self.pymol, '-cq', '-r', str(pymol_script)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        pymol_script.unlink()

        # Fix bond order/aromaticity
        sdf_free = self.workdir / f"{self.ligand_name}_free.sdf"
        template = Chem.MolFromMolFile(str(sdf_free), sanitize=True, removeHs=False)
        for f in self.workdir.glob(f"{self.docked_name}_0*.sdf"):
            pose = Chem.MolFromMolFile(str(f), sanitize=True, removeHs=False)
            if pose is None:
                raise ValueError(f"Could not parse docked pose file: {f}")

            fixed = None
            if template is not None:
                match = pose.GetSubstructMatch(template)
                if match and len(match) == template.GetNumAtoms():
                    try:
                        reordered = Chem.RenumberAtoms(pose, list(match))
                        fixed = AllChem.AssignBondOrdersFromTemplate(template, reordered)
                    except Exception as exc:
                        self.logger.warning(
                            f"Could not transfer bond orders for {f.name}: {exc}"
                        )
                else:
                    self.logger.warning(
                        f"Template/pose atom mapping failed for {f.name}; "
                        "keeping OpenBabel bond orders."
                    )
            else:
                self.logger.warning(
                    "Template ligand_free.sdf could not be parsed; "
                    "keeping OpenBabel bond orders."
                )

            Chem.MolToMolFile(fixed if fixed is not None else pose, str(f.with_suffix('.sdf')))
            subprocess.run(["sed", "-i", f"2c\\{self.ligand_md_id}", str(f.with_suffix(".sdf"))],
                check=True)

        # # Fix MOL2 files
        # perl_script = Path(__file__).parent / 'sort_mol2_bonds.pl'
        # for f in self.workdir.glob(f"{self.docked_name}_0*.mol2"):
        #     # Remove non-standard header lines
        #     subprocess.run(['sed', '-i', '1{/^#/d;}', str(f)], check=True)

        #     # Sort bonds with Perl script
        #     subprocess.run(['perl', str(perl_script), str(f), str(f)],
        #         check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        #     # Fix ligand name
        #     pattern = Rf"s/\b\(UNL\|{self.ligand_md_id}\)[[:space:]]*[0-9]*/{self.ligand_md_id}/g"
        #     subprocess.run(['sed', '-i', pattern, str(f)], check=True)

        result_dir = self.workdir / f"{self.ligand_name}.vina"
        if result_dir.exists():
            shutil.rmtree(result_dir)
        result_dir.mkdir()

        for f in self.workdir.glob('*'):
            if f.is_file() and (f.suffix == '.pdbqt' or f.name.startswith(self.ligand_name)):
                shutil.move(str(f), result_dir / f.name)
        shutil.copy2(result_dir / f"{self.ligand_name}_docked_01.sdf",
            self.workdir / f"{self.ligand_name}.sdf")

        self.logger.info(f"Docked files moved to folder: {result_dir}")


if __name__ == '__main__':
    import argparse
    import utils
    original_cwd = os.getcwd()

    parser = argparse.ArgumentParser(description="Run a redocking workflow")
    parser.add_argument('--pdb', help="Receptor ID (e.g. 1ABC or 1ABC:A:45-210)")
    parser.add_argument('--lig_new', help="New ligand SMILES")
    parser.add_argument('--lig_id', help="New ligand name")
    parser.add_argument('--lig_file', help="CSV file containing ligand SMILES")
    parser.add_argument('--verbose', action='store_true', help="Enable verbose logging")
    args = parser.parse_args()
    logger = utils.setup_logger(level=logging.INFO if args.verbose else logging.WARNING)

    lig_id = args.lig_id if args.lig_id else 'LIG'
    if args.lig_new:
        docking_jobs = [(args.pdb, lig_id, args.lig_new)]
    elif args.lig_file:
        df = pd.read_csv(args.lig_file, sep=R'\t+', engine='python')
        docking_jobs = [(row['RECEPTOR_ID'], row['LIGAND_ID'], row['LIGAND_SMILES'])
            for _, row in df.iterrows()]

    for receptor_id, ligand_id, ligand_smiles in docking_jobs:
        pdb_id = receptor_id.split(':')[0]

        workdir = Path(f"{pdb_id}_{ligand_id}").resolve()
        workdir.mkdir(exist_ok=True)

        utils.download_pdb(pdb_id, workdir)

        os.chdir(workdir)
        try:
            prep = Preprocess(receptor_id, ligand_id, ligand_smiles)
            prep.choose_ligands()
            prep.get_box()
            config_file = Path("config.yaml")
            utils.write_config_file(dict(prep.__dict__))

            # Redock pipeline
            dock = DockingPipeline(config_file)
            dock.prepare_receptor()
            dock.prepare_ligand()
            dock.dock()
            dock.postprocess()
            print('\n\n')
        except Exception as exc:
            logger.error(f"Pipeline failed in {workdir}: {exc}")
        finally:
            os.chdir(original_cwd)

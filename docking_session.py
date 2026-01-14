from email.mime import base
import importlib
import io
from pathlib import Path
import os, sys, re
import subprocess, shutil, glob
import textwrap
import logging
import yaml

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
import receptor_tools
import utils

# Import AIDDTools utils
sys.path.append(R"..\aiddtools")
import interaction_tools as inttools

class Preprocess:
    def __init__(self, parser_args):
        self.pdb_id = parser_args.pdb
        self.lig_smiles = parser_args.lig_new
        self.box = None
        self.logger = utils.get_logger()

    def choose_ligand(self):
        ligands = inttools.get_ligands(f"{self.pdb_id}.pdb")

        if not ligands:
            self.logger.warning(f"No ligands detected in {self.pdb_id}.pdb")
            return

        self.logger.info("Detected ligands:")
        for i, (lig_id, loi_status) in enumerate(ligands.items(), start=1):
            mark = f"  ({loi_status})" if loi_status == "Non-LOI" else ""
            print(f"   ({i}) Ligand {lig_id}{mark}")

        choice = int(input("\nSelect a ligand to replace: ").strip())

        keys = list(ligands.keys())
        selected_key = keys[choice - 1]

        resname, chain, _ = selected_key.split(":")

        self.lig_orig = resname
        self.chain = chain


    def get_box(self):
        pymol_script = (
            f'run C:/Users/Idener/Documents/MEGA/pymol_scripts/draw_bounding_box.py;'
            f'load {self.pdb_id}.pdb, prot;'
            f'select lig, chain {self.chain} and resn {self.lig_orig};'
            f'draw_bounding_box lig')
        exe = R"C:\Users\Idener\AppData\Local\Schrodinger\PyMOL3\Scripts\pymol.exe"
    
        result = subprocess.run([exe, "-c", "-q", "-d", pymol_script],
            check=True, capture_output=True, text=True).stdout.strip()

        lines = result.splitlines()[-2:]
        for line in lines:
            if line.startswith("Box center"):
                box_center = re.sub(R"^Box center\s*", "", line).strip()
            elif line.startswith("Box dimensions"):
                box_dims = re.sub(R"^Box dimensions\s*", "", line).strip()

        center = [float(x) for x in eval(box_center)]
        size = [float(x) for x in eval(box_dims)]
        self.box = {"center": center, "size": size}


    def write_conf_file(self):
        outfile = Path('system.conf')

        self.lig_orig = getattr(self, "self.lig_orig", "LIG")
        self.lig_name = getattr(self, "self.lig_orig", "LIG")

        # Convert numeric lists to strings
        box_center = "({:.3f}, {:.3f}, {:.3f})".format(*self.box["center"])
        box_dims = "({:.3f}, {:.3f}, {:.3f})".format(*self.box["size"])

        content = textwrap.dedent(f"""
            receptor = {self.pdb_id}
            receptor_name = protein
            chain = {self.chain}

            ligand_smiles = {self.lig_smiles}
            ligand_name = ligand
            ligand_id = {self.lig_orig}

            box_center = {box_center}
            box_dimensions = {box_dims}

            energy_range = 4
            exhaustiveness = 16
            num_modes = 8
        """).strip() + "\n"
        
        outfile.write_text(content)
        self.logger.info(f"Configuration file written: {outfile}")
        return outfile
    

class DockingSession:

    def __init__(self, params):
        self.logger = utils.get_logger()

        # Assign attributes from configuration file 
        self.pdb_id = params.get("receptor", "unknown")
        self.receptor_name = params.get("receptor_name", "protein")
        self.chain = params.get("chain", "A")

        self.lig_smiles = params.get("ligand_smiles", "")
        self.lig_id = params.get("ligand_id", "LIG")
        self.lig_name = params.get("ligand_name", "ligand")
        self.box = params.get("box", None)

        # Runtime attributes (not stored in conf)
        self.workdir = Path(".")
        self.receptor_pdbqt = None
        self.ligand_pdbqt = None


    def prepare_receptor(self, fix_loops=False):
        self.logger.info(f"Preparing receptor for {self.pdb_id} (chain {self.chain})")
        pdb = f"{self.pdb_id}.pdb"
        if fix_loops:
            receptor_tools.fix_receptor(self.pdb_id, self.chain)
        rec = receptor_tools.extract_receptor(pdb, self.chain)
        self.receptor_pdbqt = receptor_tools.convert_to_pdbqt(rec)


    def prepare_ligand(self):
        name = 'ligand'
        self.logger.info(f"Preparing ligand {name}")
        mol2_free = Path(f"{name}_free.mol2")
        pdbqt_free = Path(f"{name}_free.pdbqt")

        # Generate 3D ligand from SMILES
        subprocess.run(["obabel", f"-:{self.lig_smiles}", "-O", str(mol2_free), "--gen3d"],
                       check=True, capture_output=True)

        # Rename residues
        text = mol2_free.read_text()
        lines = text.splitlines()
        if len(lines) > 1:
            lines[1] = name.upper()[:3]
        text = "\n".join(lines)
        text = re.sub(R"\bUNL\b", name.upper()[:3], text)
        text = re.sub(R"\bLIG1\b", name.upper()[:3], text)
        text = re.sub(R"\bLIG 1\b", name.upper()[:3], text)
        text = re.sub(R"\bLIG\b", name.upper()[:3], text)
        mol2_free.write_text(text)

        # 3) Convert MOL2 to PDBQT
        subprocess.run(["obabel", str(mol2_free), "-O", str(pdbqt_free)], check=True)
        self.logger.info(f"Ligand prepared: {pdbqt_free}")

        self.ligand_pdbqt = pdbqt_free


    def dock(self, exhaustiveness=8, energy_range=4, num_modes=8, output_prefix="docked"):
        self.logger.info(f"Running docking for {self.pdb_id}")

        receptor_pdbqt = Path(self.receptor_pdbqt)
        ligand_pdbqt = Path(self.ligand_pdbqt)
        vina_conf = self.workdir / f"{output_prefix}_vina.conf"
        docked_out = self.workdir / f"{output_prefix}.pdbqt"
        log_file = self.workdir / f"{output_prefix}.log"

        center_x, center_y, center_z = self.box["center"]
        size_x, size_y, size_z = self.box["size"]

        # Write vina configuration file
        vina_conf.write_text(f"""receptor = {receptor_pdbqt}
            ligand = {ligand_pdbqt}
            center_x = {center_x}
            center_y = {center_y}
            center_z = {center_z}
            size_x = {size_x}
            size_y = {size_y}
            size_z = {size_z}
            exhaustiveness = {exhaustiveness}
            energy_range   = {energy_range}
            num_modes      = {num_modes}
            """)

        # --- Run vina via subprocess ---
        cmd = ["vina", "--config", str(vina_conf), "--out", str(docked_out)]
        subprocess.run(cmd, check=True)
        vina_conf.unlink(missing_ok=True)

        self.logger.info(f"Docking completed. Results written to {docked_out}")
        return docked_out

    def postprocess(self):
        # Split docked poses into individual mol2 files
        base_prefix = Path(self.ligand_pdbqt).stem.replace("_free", "")
        resname = base_prefix.upper()[:3]
        subprocess.run(["obabel", "docked.pdbqt", "-O", f"{base_prefix}_docked_0*.mol2", "-h"],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Add hydrogens using PyMOL
        pml_script = Path("add_hydrogens.pml")
        with pml_script.open("w") as pml:
            for f in sorted(Path(".").glob(f"{base_prefix}_docked_0*.mol2")):
                base = f.stem
                pml.write(f"load {f}, lig; h_add lig; save {base}_h.mol2, lig; delete lig;\n")
            pml.write("quit\n")
        subprocess.run(["pymol", "-cq", str(pml_script)], check=True,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Sort bonds in mol2 files
        perl_script = Path(__file__).parent / "sort_mol2_bonds.pl"
        for f in sorted(Path(".").glob(f"{base_prefix}_docked_0*_h.mol2")):
            base = f.stem.replace("_h", "")

            # Remove nonstandard header lines
            lines = f.read_text().splitlines()
            for i, line in enumerate(lines):
                if line.strip().startswith("@<TRIPOS>MOLECULE"):
                    lines = lines[i:]
                    break
            f.write_text("\n".join(lines) + "\n")

            # Call Perl sorter
            subprocess.run(["perl", str(perl_script), str(f), f"{base}.mol2"],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            mol2_text = Path(f"{base}.mol2").read_text()
            mol2_lines = mol2_text.splitlines()
            if len(mol2_lines) > 1:
                mol2_lines[1] = resname
            mol2_text = "\n".join(mol2_lines)
            mol2_text = re.sub(R"\b(UNL[0-9A-Za-z_]*|LIG[0-9A-Za-z_ ]*)\b", resname, mol2_text)
            Path(f"{base}.mol2").write_text(mol2_text)

        for f in glob.glob(f"{base_prefix}_docked*_h.mol2"):
            os.remove(f)
        pml_script.unlink()

        result_dir = Path(f"{base_prefix}.vina")
        if result_dir.exists():
            shutil.rmtree(result_dir)
        result_dir.mkdir()

        for f in Path(".").glob("*"):
            if f.is_file() and (f.suffix == ".pdbqt" or f.name.startswith(base_prefix)):
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
    parser.add_argument("--lig_new", required=True, help="New ligand SMILES")
    args = parser.parse_args()

    prep = Preprocess(args)
    prep.choose_ligand()
    prep.get_box()
    conf_file = prep.write_conf_file()

    params = utils.read_conf_file("system.conf")

    redock = RedockSession(params)
    result = redock.run()

import shutil

from rdkit import Chem
from rdkit.Chem import AllChem
from pymol2 import PyMOL
from .tools import run_in_terminal


def postprocess(ctx):
    # Split docked poses into individual SDF files
    run_in_terminal([ctx.obabel, str(ctx.docked_file),
        '-O', str(ctx.workdir / f"{ctx.docked_name}_0*.sdf")])

    # Remove hydrogens using PyMOL (the only reliable way)
    with PyMOL() as pm:
        for f in ctx.workdir.glob(f"{ctx.docked_name}_0*.sdf"):
            pm.cmd.load(str(f), "lig")
            pm.cmd.remove("lig and hydro")
            pm.cmd.save(str(f), "lig")
            pm.cmd.delete("lig")

    template = Chem.MolFromSmiles(ctx.ligand_smiles)

    for f in ctx.workdir.glob(f"{ctx.docked_name}_0*.sdf"):
        pose = Chem.MolFromMolFile(str(f), sanitize=False)

        # Fix bond order/aromaticity
        fixed_pose = AllChem.AssignBondOrdersFromTemplate(template, pose)

        # Reorder atoms to match template
        match = fixed_pose.GetSubstructMatch(template)
        reordered_pose = Chem.RenumberAtoms(fixed_pose, list(match))

        # Build pose on template graph, then force docked coords
        checked_pose = Chem.Mol(template)
        checked_pose.AddConformer(Chem.Conformer(reordered_pose.GetConformer()))
        checked_pose = Chem.AddHs(checked_pose, addCoords=True)

        Chem.MolToMolFile(checked_pose, str(f.with_suffix('.sdf')))
        run_in_terminal(["sed", "-i", f"2c\\{ctx.ligand_resname}", str(f.with_suffix('.sdf'))])

    result_dir = ctx.workdir / f"{ctx.ligand_name}.vina"
    if result_dir.exists():
        shutil.rmtree(result_dir)
    result_dir.mkdir()

    for f in ctx.workdir.glob('*'):
        if f.is_file() and (f.suffix == '.pdbqt' or f.name.startswith(ctx.ligand_name)):
            shutil.move(str(f), result_dir / f.name)
    shutil.copy2(
        result_dir / f"{ctx.ligand_name}_docked_01.sdf",
        ctx.workdir / f"{ctx.ligand_name}.sdf",)

    print(f"Docked files moved to folder: {result_dir}")
    return ctx

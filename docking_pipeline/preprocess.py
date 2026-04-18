from pathlib import Path

from . import ligand as ligand_stage
from pymol2 import PyMOL


def check_ionisation(ctx):
    print("Calculating protonation states...")
    result_dir = ctx.workdir / f"{ctx.ligand_name}.unipka"
    result_dir.mkdir(exist_ok=True)
    smi = ligand_stage.check_ionisation(
        ctx.ligand_smiles,
        report_file=result_dir / "unipka_report.txt",)
    ctx.ligand_smiles = smi
    return ctx


def get_box(ctx):
    draw_box_script = Path(__file__).resolve().parents[1] / "pymol_scripts" / "draw_bounding_box.py"
    with PyMOL() as pm:
        pm.cmd.do(f"run {draw_box_script}")
        pm.cmd.load(f"{ctx.pdb_id}.pdb", "prot")
        pm.cmd.select("lig", f"chain {ctx.chain} and resn {ctx.cognate_lig_id}")
        extents = pm.cmd.get_extent("lig")

    min_corner, max_corner = extents
    cx = (min_corner[0] + max_corner[0]) / 2.0
    cy = (min_corner[1] + max_corner[1]) / 2.0
    cz = (min_corner[2] + max_corner[2]) / 2.0
    sx = max_corner[0] - min_corner[0]
    sy = max_corner[1] - min_corner[1]
    sz = max_corner[2] - min_corner[2]
    ctx.box_center = (float(cx), float(cy), float(cz))
    ctx.box_size = (float(sx), float(sy), float(sz))
    return ctx

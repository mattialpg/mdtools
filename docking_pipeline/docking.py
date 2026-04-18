from .tools import run_checked


def run_docking(ctx, output_prefix='docked'):
    print(f"Running docking for {ctx.receptor_id}")

    ctx.docked_name = f"{ctx.ligand_name}_{output_prefix}"
    ctx.docked_file = ctx.workdir / f"{ctx.docked_name}.pdbqt"
    vina_conf = ctx.workdir / f"{ctx.docked_name}_vina.conf"
    log_file = ctx.workdir / f"{ctx.docked_name}.log"

    cx, cy, cz = ctx.box_center
    sx, sy, sz = ctx.box_size

    vina_conf.write_text(f"""
        receptor = {ctx.receptor_pdbqt}
        ligand = {ctx.ligand_pdbqt}
        center_x = {cx}
        center_y = {cy}
        center_z = {cz}
        size_x = {sx}
        size_y = {sy}
        size_z = {sz}
        exhaustiveness = {ctx.exhaustiveness}
        energy_range = {ctx.energy_range}
        num_modes = {ctx.num_modes}""")

    run_checked([ctx.vina, '--config', str(vina_conf), '--out', str(ctx.docked_file),
        '--log', str(log_file)])
    vina_conf.unlink(missing_ok=True)

    print(f"Docking completed. Results written to {ctx.docked_file}")
    return ctx

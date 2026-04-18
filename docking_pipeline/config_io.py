from pathlib import Path
import yaml


def build_config(ctx):
    return {
        "protein": {"id": ctx.receptor_id,
            "name": ctx.receptor_name,},
        "ligands": ctx.ligands,
        "box_center": str(tuple(ctx.box_center)),
        "box_size": str(tuple(ctx.box_size)),
        "energy_range": ctx.energy_range,
        "exhaustiveness": ctx.exhaustiveness,
        "num_modes": ctx.num_modes,
        "workdir": str(ctx.workdir),
        "pymol": ctx.pymol,
        "obabel": ctx.obabel,
        "vina": ctx.vina,}


def write_config(ctx, config_file=Path("config.yaml")):
    config = build_config(ctx)
    config_text = yaml.safe_dump(config, sort_keys=False)
    for key in ("ligands:", "box_center:", "workdir:"):
        config_text = config_text.replace(f"\n{key}", f"\n\n{key}")
    config_file.write_text(config_text)


def read_config(config_file=Path("config.yaml")):
    return yaml.safe_load(Path(config_file).read_text()) or {}

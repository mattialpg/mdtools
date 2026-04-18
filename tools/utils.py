import sys
import numpy as np
from pathlib import Path
import logging, requests
import yaml
from ast import literal_eval


TOOL_PATHS = {
    "pymol": "/home/mattia/.miniforge/envs/my-chem/bin/pymol",
    "obabel": "/home/mattia/.miniforge/envs/my-chem/bin/obabel",
    "vina": "/home/mattia/.bin/smina.static",}


def setup_logger(name=None, level=logging.INFO):
    """Create a controlled logger with consistent handlers."""
    root = logging.getLogger()
    if not getattr(root, "_mdtools_logger_configured", False):
        root.handlers.clear()
        formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')

        # Console
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(formatter)
        root.addHandler(stream_handler)

        # Fatal-on-error handler
        def exit_on_error(record):
            if record.levelno >= logging.ERROR:
                sys.exit(1)
        fatal_handler = logging.Handler()
        fatal_handler.emit = exit_on_error
        root.addHandler(fatal_handler)
        root._mdtools_logger_configured = True

    root.setLevel(level)
    logger = logging.getLogger(name) if name else logging.getLogger("mdtools")
    logger.setLevel(level)

    logging.getLogger("plip").setLevel(logging.WARNING)
    logging.getLogger("rdkit").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.ERROR)

    return logger


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)
    return path

def download_pdb(pdb_id, outdir="."):
    # pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    outpath = Path(outdir) / f"{pdb_id}.pdb"

    response = requests.get(url)
    if response.status_code != 200:
        raise RuntimeError(f"Failed to download {pdb_id} (HTTP {response.status_code})")

    outpath.write_text(response.text)
    return outpath


def get_tool_paths(config_file=Path("config.yaml")):
    tools = dict(TOOL_PATHS)
    path = Path(config_file)
    if path.exists():
        configs = yaml.safe_load(path.read_text()) or {}
        for key in ("pymol", "obabel", "vina"):
            if configs.get(key):
                tools[key] = configs[key]
    return tools


def write_config_file(preprocess_dict):
    cwd = Path.cwd()
    config_file = Path('config.yaml')
    protein = preprocess_dict.get('protein')
    if not isinstance(protein, dict):
        protein = {
            'id': preprocess_dict['receptor_id'],
            'name': preprocess_dict['receptor_name'],}

    configs = {
        'protein': protein,
        
        "ligands": preprocess_dict["ligands"],

        'box_center': preprocess_dict['box_center'],
        'box_size': preprocess_dict['box_size'],
        'energy_range': 4,
        'exhaustiveness': 16,
        'num_modes': 8,

        'workdir': str(cwd),
        'pymol': preprocess_dict.get('pymol', TOOL_PATHS['pymol']),
        'obabel': preprocess_dict.get('obabel', TOOL_PATHS['obabel']),
        'vina': preprocess_dict.get('vina', TOOL_PATHS['vina'])}
        
    config_text = yaml.safe_dump(configs, sort_keys=False)
    for key in ('ligands:', 'box_center:', 'workdir:'):
        config_text = config_text.replace(f"\n{key}", f"\n\n{key}")
    config_file.write_text(config_text)

    # self.logger.info(f"Configuration file written: {config_file}")
    return configs


def read_config_file(config_file=Path("config.yaml")):
    configs = yaml.safe_load(Path(config_file).read_text()) or {}

    protein = configs.get("protein") or {}
    ligands = configs.get("ligands")
    ligand = {}
    if isinstance(ligands, list) and ligands and isinstance(ligands[0], dict):
        ligand = ligands[0]
    elif isinstance(ligands, dict):
        ligand = ligands
    if isinstance(protein, dict):
        if protein.get("id") and not configs.get("receptor_id"):
            configs["receptor_id"] = protein["id"]
        if protein.get("name") and not configs.get("receptor_name"):
            configs["receptor_name"] = protein["name"]
    if isinstance(ligand, dict):
        if ligand.get("id") and not configs.get("ligand_id"):
            configs["ligand_id"] = ligand["id"]
        if ligand.get("smiles") and not configs.get("ligand_smiles"):
            configs["ligand_smiles"] = ligand["smiles"]
        if ligand.get("name") and not configs.get("ligand_name"):
            configs["ligand_name"] = ligand["name"]
        if ligand.get("resname") and not configs.get("ligand_resname"):
            configs["ligand_resname"] = ligand["resname"]

    configs.setdefault("receptor_name", "protein")
    configs.setdefault("ligand_name", "ligand")
    configs.setdefault("ligand_resname", "LIG")

    receptor_id = configs.get("receptor_id")
    if receptor_id:
        parts = str(receptor_id).split(":")
        configs["pdb_id"] = parts[0]
        configs["chain"] = parts[1] if len(parts) > 1 else None
        configs["domain"] = parts[2] if len(parts) > 2 else None

    for key in ("box_center", "box_size"):
        value = configs.get(key)
        if isinstance(value, str):
            try:
                configs[key] = tuple(float(x) for x in literal_eval(value))
            except Exception:
                pass

    for key, value in TOOL_PATHS.items():
        configs.setdefault(key, value)

    configs["protein"] = {
        "id": configs.get("receptor_id"),
        "name": configs.get("receptor_name"),
    }
    ligand_entry = {
        "id": configs.get("ligand_id"),
        "smiles": configs.get("ligand_smiles"),
        "name": configs.get("ligand_name"),
        "resname": configs.get("ligand_resname"),
    }
    configs["ligands"] = [ligand_entry]
    # Backward compatibility for code paths that still read singular key.
    configs["ligand"] = ligand_entry

    return configs


def round_up_nice(x):
    """Round up to a 'nice' number, excluding 7 and 9 as leading digits.
    If x > 10, snap to the next multiple of 5."""
    if x == 0:
        return 0

    if x > 10:
        # Snap to the next multiple of 5
        nice = np.ceil(x / 5) * 5
        # Check the leading digit after rounding
        leading = int(str(int(nice))[0])
        if leading in (7, 9):
            nice += 5  # move to next multiple of 5
        return nice

    # For small numbers (<10), round up as before
    exponent = np.floor(np.log10(x))
    fraction = x / 10**exponent
    nice_fraction = np.ceil(fraction)
    if nice_fraction in (7, 9):
        nice_fraction += 1
        if nice_fraction >= 10:
            nice_fraction = 1
            exponent += 1
    return np.round(nice_fraction * 10**exponent, 2)


def nice_ticks(x1, x2, min_ticks=4, max_ticks=6):
    """Return 'nice' matching tick values for two axes.
    Both axes will have the same number of ticks (4–6 inclusive).
    The best pair is chosen by combined niceness score."""
    if x1 <= 0 and x2 <= 0:
        return 0, np.array([0]), 0, np.array([0])

    u1 = round_up_nice(x1)
    u2 = round_up_nice(x2)

    # Generate candidate tick lists for both axes
    results1, results2 = {}, {}
    for n_ticks in range(min_ticks, max_ticks + 1):
        for (u, results) in [(u1, results1), (u2, results2)]:
            ticks = np.round(np.linspace(0, u, n_ticks), 2)
            results[n_ticks] = ticks

    # Score each tick list
    scores1, scores2 = {}, {}
    for results, scores in [(results1, scores1), (results2, scores2)]:
        for n_ticks, ticks in results.items():
            step = np.round(ticks[1] - ticks[0], 2)
            score = 0

            # Base score for step simplicity
            if step in [1, 2, 2.5, 5, 10, 0.5, 0.25, 0.2, 0.1]:
                score += 3
            elif step in [0.05, 0.02, 0.01]:
                score += 2
            else:
                score -= 2  # penalise awkward like 0.13

            # Bonus if all tick decimals are 0 or .00/.05/.10/.20 etc. (even decimals)
            fractional_parts = np.round((ticks * 100) % 10, 0)  # last digit
            if np.all(np.isin(fractional_parts, [0, 2, 4, 6, 8])):
                score += 1

            # Penalise mixture of odd/even decimals
            elif np.any(np.isin(fractional_parts, [1, 3, 5, 7, 9])):
                score -= 1

            scores[n_ticks] = score


    # Compare tick lists with same tick count and sum their scores
    best_key, best_score = None, -np.inf
    for n_ticks in results1.keys() & results2.keys():
        total_score = scores1[n_ticks] + scores2[n_ticks]
        if total_score > best_score:
            best_key = n_ticks
            best_score = total_score
    ticks1 = results1[best_key]
    ticks2 = results2[best_key]

    return u1, ticks1, u2, ticks2


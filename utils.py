import sys
import numpy as np
from pathlib import Path
import logging, requests


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)
    return path

def download_pdb(pdb_id, outdir="."):
    pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    outpath = Path(outdir) / f"{pdb_id}.pdb"

    response = requests.get(url)
    if response.status_code != 200:
        raise RuntimeError(f"Failed to download {pdb_id} (HTTP {response.status_code})")

    outpath.write_text(response.text)
    return outpath


def setup_logger(logfile: Path):
    logfile.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=logfile,
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO,
    )
    return logging.getLogger(__name__)


def get_logger() -> logging.Logger:
    import logging, sys

    logger = logging.getLogger()
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)

        logging.getLogger("plip").setLevel(logging.WARNING)
        logging.getLogger("rdkit").setLevel(logging.WARNING)
        logging.getLogger("urllib3").setLevel(logging.ERROR)

    logger.propagate = False
    return logger

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


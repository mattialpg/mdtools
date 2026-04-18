import argparse
import logging

import pandas as pd

import tools.utils as utils
from .context import RunContext
from .runner import PipelineRunner


def parse_args():
    parser = argparse.ArgumentParser(description="Run a redocking workflow")
    parser.add_argument('--pdb', help="Receptor ID (e.g. 1ABC or 1ABC:A:45-210)")
    parser.add_argument('--lig_smiles', help="New ligand SMILES")
    parser.add_argument('--lig_id', default="LIG", help="New ligand name")
    parser.add_argument('--csv', help="CSV file containing ligand SMILES")
    parser.add_argument('--verbose', action='store_true', help="Enable verbose logging")
    return parser.parse_args()


def build_jobs(args):
    if args.lig_smiles:
        return [(args.pdb, args.lig_id, args.lig_smiles)]
    if args.csv:
        df = pd.read_csv(args.csv, sep=R'\t+', engine='python')
        return [(row['RECEPTOR_ID'], row['LIGAND_ID'], row['LIGAND_SMILES'])
            for _, row in df.iterrows()]
    return []


def main():
    args = parse_args()
    jobs = build_jobs(args)

    logger = utils.setup_logger(level=logging.INFO if args.verbose else logging.WARNING)
    runner = PipelineRunner(logger=logger)

    for job in jobs:
        ctx = RunContext.from_job(job)
        runner.run_job(ctx)


if __name__ == '__main__':
    main()

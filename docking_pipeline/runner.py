import os
from pathlib import Path

import tools.utils as utils

from . import config_io
from .ligand import select_cognate_ligand, prepare_ligand
from .preprocess import check_ionisation, get_box
from .receptor import prepare_receptor
from .docking import run_docking
from .postprocess import postprocess

class PipelineRunner:
    def __init__(self, logger):
        self.logger = logger
        self.original_cwd = Path.cwd()

    def run_job(self, ctx):
        workdir = Path(ctx.workdir).resolve()
        ctx.workdir = workdir
        workdir.mkdir(exist_ok=True)

        os.chdir(workdir)
        utils.download_pdb(ctx.pdb_id)

        try:
            for stage in [select_cognate_ligand, check_ionisation, get_box,
                prepare_receptor, prepare_ligand, run_docking, postprocess]:
                ctx = stage(ctx)
            
            config_io.write_config(ctx)
            
        except Exception as exc:
            self.logger.error(f"Pipeline failed in {workdir}: {exc}")
        finally:
            os.chdir(self.original_cwd)

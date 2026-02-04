#!/bin/bash
set -euo pipefail
trap 'echo "Error on line $LINENO"; exit 1' ERR

bash /home/mattia/mdtools/gromacs_scripts/make_topology.sh ligand
bash /home/mattia/mdtools/gromacs_scripts/make_equil.sh ligand
bash /home/mattia/mdtools/gromacs_scripts/make_prod.sh -length 50
bash /home/mattia/mdtools/gromacs_scripts/make_trj.sh -i md_50

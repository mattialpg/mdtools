#!/bin/bash
set -euo pipefail
trap 'echo "Error on line $LINENO"; exit 1' ERR

bash ~/.gromacs/make_topology.sh ligand.mol2
bash ~/.gromacs/make_equil.sh ligand
bash ~/.gromacs/make_prod.sh -length 200
bash ~/.gromacs/make_trj.sh -i md_200

#!/bin/bash
set -euo pipefail

# Usage:
#   ./run_all.sh -len NS

# Parse command-line args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -len) LENGTH="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 -len NS"
      exit 0
      ;;
    *) echo "Usage: $0 -len NS"; exit 1 ;;
  esac
done

if [[ -z "${LENGTH:-}" ]]; then
  echo "Usage: $0 -len NS"
  exit 1
fi

MDNAME="md_${LENGTH}"

bash /home/mattia/mdtools/gromacs_scripts/make_topology.sh -lig ligand
bash /home/mattia/mdtools/gromacs_scripts/make_equil.sh ligand
bash /home/mattia/mdtools/gromacs_scripts/make_prod.sh -len "${LENGTH}"
bash /home/mattia/mdtools/gromacs_scripts/make_trj.sh -i "${MDNAME}"

#!/bin/bash
set -euo pipefail
shopt -s nullglob

# Usage:
#   ./make_lig_analysis.sh [--config CONFIG_YAML]

# Parse command-line args
NET_CHARGE=0
CONFIG_YAML="config.yaml"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG_YAML="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 [--config CONFIG_YAML]"
      exit 0
      ;;
    *) echo "Usage: $0 [--config CONFIG_YAML]"; exit 1 ;;
  esac
done

MD_LENGTH=$(yq -r '.md_length // empty' "${CONFIG_YAML}")
MDNAME="md_${MD_LENGTH}"
TRJNAME="trj_${MD_LENGTH}"
LIG="LIG"

log() {
  printf "\n\033[38;2;255;255;255;48;2;15;88;157m%s\033[0m\n\n" "$1"
}

# Calculate ligand RMSD and pocket-ligand distance on full trajectory
echo -e "q" | gmx make_ndx -f "${MDNAME}.gro" -o index.ndx
printf "${LIG}\n" | gmx rms -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" \
  -n index.ndx -o trj_rmsd_lig.xvg -tu ns -fit none
gmx select -s "${MDNAME}.tpr" -n index.ndx -on pocket.ndx \
  -select "group \"Protein\" and same residue as (within 0.45 of group \"${LIG}\"); group \"${LIG}\""
sed -i '0,/^\[.*\]$/s//[ Pocket ]/' pocket.ndx
gmx pairdist -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" -n pocket.ndx \
  -o trj_dist_lig.xvg -tu ns -ref 'group "Pocket"' -sel "group \"${LIG}\"" -type min

python ~/mdtools/interaction_analysis.py --trj_name "${TRJNAME}_fit"
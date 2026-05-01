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

# Calculate ligand RMSD and pocket-ligand distance
echo -e "q" | gmx make_ndx -f "${MDNAME}.gro" -o index.ndx
printf "${LIG}\n" | gmx rms -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" \
  -n index.ndx -o trj_rmsd_lig.xvg -tu ns -fit none
gmx select -s "${MDNAME}.tpr" -n index.ndx -on pocket.ndx \
  -select "group \"Protein\" and same residue as (within 0.45 of group \"${LIG}\"); group \"${LIG}\""
sed -i '0,/^\[.*\]$/s//[ Pocket ]/' pocket.ndx
gmx pairdist -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" -n pocket.ndx \
  -o trj_dist_lig.xvg -tu ns -ref 'group "Pocket"' -sel "group \"${LIG}\"" -type min

python ~/mdtools/interaction_analysis.py --trj_name "${TRJNAME}_strip"

### ---------------------------- ###
###       Binding Affinity       ###
### ---------------------------- ###

MMPBSA_DIR="ligand.mmpbsa"
mkdir -p "${MMPBSA_DIR}"
echo -e "q" | gmx make_ndx -f "${TRJNAME}_strip.gro" -o "${MMPBSA_DIR}/mmpbsa.ndx"
gmx select -s "${TRJNAME}_strip.gro" -n "${MMPBSA_DIR}/mmpbsa.ndx" \
  -on "${MMPBSA_DIR}/complex.ndx" -select "group \"Protein\" or group \"${LIG}\""
sed -i '0,/^\[.*\]$/s//[ Protein_LIG ]/' "${MMPBSA_DIR}/complex.ndx"
printf "Protein_LIG\n" | gmx trjconv -f "${TRJNAME}_strip.gro" -s "${TRJNAME}_strip.gro" \
  -o "${MMPBSA_DIR}/${TRJNAME}_strip.pdb" -n "${MMPBSA_DIR}/complex.ndx" 

# Assign chain to receptor (A) and ligand (B)
awk 'BEGIN{OFS=""}(/^ATOM/ || /^HETATM/){
  resn = substr($0,18,3) chain = (resn == "LIG" ? "B" : "A")
  $0 = sprintf("%s%s%s", substr($0,1,21), chain, substr($0,23))}{print}
' "${MMPBSA_DIR}/${TRJNAME}_strip.pdb" > "${MMPBSA_DIR}/strip.tmp"
mv "${MMPBSA_DIR}/strip.tmp" "${MMPBSA_DIR}/${TRJNAME}_strip.pdb"

# Create MMPBSA input
( cd "${MMPBSA_DIR}" && gmx_MMPBSA --create_input gb >/dev/null )
sed -i 's/\(forcefields *= *\)\".*\"/\1\"oldff\/leaprc.ff99SBildn,leaprc.gaff\"/' "${MMPBSA_DIR}/mmpbsa.in"
sed -i 's/\(saltcon *= *\)0\.0/\10.150/' "${MMPBSA_DIR}/mmpbsa.in"

# Run MMPBSA
gmx_MMPBSA -O -nogui -i "${MMPBSA_DIR}/mmpbsa.in" -cs "${MDNAME}.tpr" \
  -cr "${MMPBSA_DIR}/${TRJNAME}_strip.pdb" \
  -ct "${TRJNAME}_strip.xtc" -ci "${MMPBSA_DIR}/mmpbsa.ndx" -cg Protein "${LIG}" \
  -cp topol.top -o "${MMPBSA_DIR}/FINAL_RESULTS_MMPBSA.dat" \
  -eo "${MMPBSA_DIR}/FINAL_RESULTS_MMPBSA.csv"
mv COMPACT_MMXSA_RESULTS.mmxsa gmx_MMPBSA.log "${MMPBSA_DIR}"

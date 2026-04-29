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

# # Calculate ligand RMSD and pocket-ligand distance
# echo -e "q" | gmx make_ndx -f "${MDNAME}.gro" -o index.ndx
# printf "${LIG}\n" | gmx rms -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" \
#   -n index.ndx -o trj_rmsd_lig.xvg -tu ns -fit none
# gmx select -s "${MDNAME}.tpr" -n index.ndx -on pocket.ndx \
#   -select "group \"Protein\" and same residue as (within 0.45 of group \"${LIG}\"); group \"${LIG}\""
# sed -i '0,/^\[.*\]$/s//[ Pocket ]/' pocket.ndx
# gmx pairdist -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" -n pocket.ndx \
#   -o trj_dist_lig.xvg -tu ns -ref 'group "Pocket"' -sel "group \"${LIG}\"" -type min

# python ~/mdtools/interaction_analysis.py --trj_name "${TRJNAME}_strip"

#
MMPBSA_DIR="ligand.mmpbsa"
mkdir -p "${MMPBSA_DIR}"
echo -e "q" | gmx make_ndx -f "${TRJNAME}_strip.gro" -o "${MMPBSA_DIR}/index_mmpbsa.ndx"
gmx select -s "${TRJNAME}_strip.gro" -n "${MMPBSA_DIR}/index_mmpbsa.ndx" \
  -on "${MMPBSA_DIR}/complex.ndx" -select "group \"Protein\" or group \"${LIG}\""
sed -i '0,/^\[.*\]$/s//[ Protein_LIG ]/' "${MMPBSA_DIR}/complex.ndx"

# Build a Protein+LIG-only complex PDB for gmx_MMPBSA and annotate chains.
printf "Protein_${LIG}\n" | gmx trjconv -f "${TRJNAME}_strip.gro" -s "${TRJNAME}_strip.gro" \
  -n "${MMPBSA_DIR}/complex.ndx" -o "${MMPBSA_DIR}/complex_prot_lig.pdb"

# Chain A: receptor atoms, Chain B: ligand atoms (resname LIG).
awk '
BEGIN{OFS=""}
(/^ATOM/ || /^HETATM/){
  resn = substr($0,18,3)
  chain = (resn == "LIG" ? "B" : "A")
  $0 = sprintf("%s%s%s", substr($0,1,21), chain, substr($0,23))
}
{print}
' "${MMPBSA_DIR}/complex_prot_lig.pdb" > "${MMPBSA_DIR}/complex_prot_lig_chain.pdb"

cat > "${MMPBSA_DIR}/mmpbsa.in" << EOF
&general
  startframe=1, endframe=999999, interval=10,
  verbose=1, keep_files=0,
/
&gb
  igb=5, saltcon=0.150
/
EOF

gmx_MMPBSA -O -nogui -i "${MMPBSA_DIR}/mmpbsa.in" -cs "${MDNAME}.tpr" \
  -cr "${MMPBSA_DIR}/complex_prot_lig_chain.pdb" \
  -ct "${TRJNAME}_strip.xtc" -ci "${MMPBSA_DIR}/index_mmpbsa.ndx" -cg Protein "${LIG}" \
  -cp topol.top -o "${MMPBSA_DIR}/FINAL_RESULTS_MMPBSA.dat" \
  -eo "${MMPBSA_DIR}/FINAL_RESULTS_MMPBSA.csv"

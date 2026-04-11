#!/bin/bash
set -euo pipefail
shopt -s nullglob

# Usage:
#   ./make_topology.sh [--config CONFIG_YAML]

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

log() {
  printf "\n\033[38;2;255;255;255;48;2;15;88;157m%s\033[0m\n\n" "$1"
}

prepare_ligand_top() {
  local LIGAND_NAME="$1"
  local ligand_md_id="$2"
  local output_name="${LIGAND_NAME}_amber"

  cat > tleap_topology.in <<EOF
source leaprc.gaff2
loadamberparams ${output_name}.frcmod
LIGAND = loadmol2 ${output_name}.mol2
check LIGAND
saveamberparm LIGAND ${output_name}.prmtop ${output_name}.rst7
quit
EOF

  log " >  Calculating charges for ${LIGAND_NAME} ($ligand_md_id)..."
  antechamber -i "${LIGAND_NAME}.sdf" -fi sdf -o "${output_name}.mol2" -fo mol2 -c bcc -s 2 -at gaff2 -nc "${NET_CHARGE}" -m 1
  sed -i "s/\<MOL\>/${ligand_md_id}/g" "${output_name}.mol2"
  rm -rf "${LIGAND_NAME}.antechamber"; mkdir "${LIGAND_NAME}.antechamber"
  mv ANTECHAMBER_* ATOMTYPE* sqm.* "${LIGAND_NAME}.antechamber"

  log " >  Generating Amber topology for ${LIGAND_NAME}..."
  parmchk2 -i "${output_name}.mol2" -f mol2 -o "${output_name}.frcmod" -s gaff2
  tleap -f tleap_topology.in

  log " >  Converting ${LIGAND_NAME} to GROMACS topology..."
  acpype -i "${output_name}.mol2" 2> >(grep -Ev \
  "OpenBabel|Open Babel|Cannot perform atom type translation|
  This Mol2 file is non-standard|Cannot interpret atom type|
  ==============================" >&2)
  rm -rf "${LIGAND_NAME}.acpype"; mv -f "${output_name}.acpype" "${LIGAND_NAME}.acpype"
  rename -f 's/_amber//' "${LIGAND_NAME}.acpype"/*

  awk '/\[ atomtypes \]/,/\[ moleculetype \]/ {if ($0 !~ /\[ moleculetype \]/) print}' "${LIGAND_NAME}.acpype/${LIGAND_NAME}_GMX.itp" > "${LIGAND_NAME}.prm"
  awk '/\[ moleculetype \]/,/\[ system \]/ {if ($0 !~ /\[ system \]/) print}' "${LIGAND_NAME}.acpype/${LIGAND_NAME}_GMX.itp" > "${LIGAND_NAME}.itp"
  sed -i "s/\b${output_name}\b/${ligand_md_id}/g" "${LIGAND_NAME}.itp"

  mv -f tleap_topology.in "${output_name}".* leap.log "${LIGAND_NAME}.acpype/"
  cp "${LIGAND_NAME}.acpype/${LIGAND_NAME}_GMX.gro" "${LIGAND_NAME}.gro"

  # Replace GRO residue name with ligand_md_id
  awk -v new_resn="${ligand_md_id}" '
  NR <= 2 { print; next }
  NR == total { print; next }
  {
    resn = sprintf("%-5.5s", new_resn)
    print substr($0, 1, 5) resn substr($0, 11)
  }
' total="$(wc -l < "${LIGAND_NAME}.gro")" "${LIGAND_NAME}.gro" > "${LIGAND_NAME}.gro.tmp"
  mv "${LIGAND_NAME}.gro.tmp" "${LIGAND_NAME}.gro"

  obabel "${LIGAND_NAME}.sdf" -O "${LIGAND_NAME}.pdb" --writeconect
  sed -i "s/\<UNNAMED\>/${ligand_md_id}/g" "${LIGAND_NAME}.pdb"
  sed -i "s/\<UNL\>/${ligand_md_id}/g" "${LIGAND_NAME}.pdb"
}

# Iterate over ligands
while IFS=$'\t' read -r LIGAND_NAME LIGAND_MD_ID; do
  [[ -z "$LIGAND_NAME" ]] && continue
  if [[ -z "${LIGAND_MD_ID:-}" || "${LIGAND_MD_ID}" == "null" ]]; then
    echo "Error: missing md_id for ligand '${LIGAND_NAME}' in ${CONFIG_YAML}"
    exit 1
  fi
  prepare_ligand_top "$LIGAND_NAME" "$LIGAND_MD_ID"
done < <(yq -r '.ligands[] | [.name, .md_id] | @tsv' "${CONFIG_YAML}")

#!/bin/bash
set -euo pipefail
shopt -s nullglob

# Usage:
#   ./make_pocket_analysis.sh [--config CONFIG_YAML]
# Ref: https://www.cell.com/biophysj/fulltext/S0006-3495(26)00154-2

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


MD_LENGTH=$(yq -r '.md_length // empty' "${CONFIG_YAML}")
MDNAME="md_${MD_LENGTH}"
TRJNAME="trj_${MD_LENGTH}"

STRIP_XTC="${TRJNAME}_strip.xtc"
STRIP_GRO="${TRJNAME}_strip.gro"
NOLIG_XTC="${TRJNAME}_nolig.xtc"
NOLIG_GRO="${TRJNAME}_nolig.gro"


printf "q\n" | gmx make_ndx -f "${STRIP_GRO}" -o strip.ndx
printf "C-alpha\n" | gmx trjconv -f "${STRIP_GRO}" -s "${STRIP_GRO}" \
  -o "${NOLIG_GRO}" -n strip.ndx
printf "C-alpha\n" | gmx trjconv -f "${STRIP_XTC}" -s "${STRIP_GRO}" \
  -o "${NOLIG_XTC}" -n strip.ndx -ndec 2

# Calculate C-alpha statistics
printf "q\n" | gmx make_ndx -f "${NOLIG_GRO}" -o nolig.ndx
printf "C-alpha\nC-alpha\n" | gmx rms -f "${NOLIG_XTC}" -s "${NOLIG_GRO}" \
  -o trj_rmsd_prot.xvg -n nolig.ndx -tu ns
printf "C-alpha\n" | gmx rmsf -f "${NOLIG_XTC}" -s "${NOLIG_GRO}" \
  -o trj_rmsf_prot.xvg -n nolig.ndx -res
gmx dssp -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" -o trj_dssp.dat -tu ns
printf "group \"C-alpha\"\n" | gmx gyrate -f "${NOLIG_XTC}" -s "${NOLIG_GRO}" \
  -o trj_gyrate_prot.xvg -n nolig.ndx
# printf "C-alpha\n" | gmx sasa -s "${NOLIG_GRO}" -f "${NOLIG_XTC}" \
#   -n nolig.ndx -o trj_sasa_prot.xvg -or trj_sasa_residues.xvg

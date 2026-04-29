#!/bin/bash
set -euo pipefail

# Usage:
#   ./make_trj.sh [--config CONFIG_YAML]

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

# Use temporary storage for intermediate trajectory files
WORK_TMPDIR="${TMPDIR:-/tmp}"
RUN_TMPDIR="$(mktemp -d "${WORK_TMPDIR%/}/make_trj.XXXXXX")"
SAMPLED_XTC="${TRJNAME}_sampled.xtc"
WHOLE_XTC="${RUN_TMPDIR}/${TRJNAME}_whole.xtc"
NOJUMP_XTC="${RUN_TMPDIR}/${TRJNAME}_nojump.xtc"
NOROT_XTC="${RUN_TMPDIR}/${TRJNAME}_norot.xtc"
FIT_XTC="${RUN_TMPDIR}/${TRJNAME}_fit.xtc"
STRIP_XTC="${TRJNAME}_strip.xtc"
STRIP_TPR="${RUN_TMPDIR}/${TRJNAME}_strip.tpr"

if [[ ! -s "${MDNAME}.gro" ]]; then
  echo "0" | gmx trjconv -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" -o "${MDNAME}.gro" -dump -1
fi

# Prepare stripped files
echo -e "q" | gmx make_ndx -f wrapped.gro -o index.ndx
if ! grep -Eq '^\[[[:space:]]*K[[:space:]]*\]$' index.ndx || ! grep -Eq '^\[[[:space:]]*CL[[:space:]]*\]$' index.ndx; then
  echo -e "r K\nr CL\nq" | gmx make_ndx -f wrapped.gro -o index.ndx
fi
echo -e '"non-Water" & ! "K" & ! "CL"\nq' | gmx make_ndx -f wrapped.gro \
  -o index.ndx -n index.ndx
echo 'non-Water_&_!K_&_!CL' | gmx trjconv -f wrapped.gro -s "${MDNAME}.tpr" \
  -o "${TRJNAME}_strip.gro" -n index.ndx
echo "non-Water_&_!K_&_!CL" | gmx convert-tpr -s "${MDNAME}.tpr" \
  -o "${STRIP_TPR}" -n index.ndx
echo -e "q" | gmx make_ndx -f "${TRJNAME}_strip.gro" -o strip.ndx

# Downsample full system
printf "System\n" | gmx trjconv -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" \
  -o "${SAMPLED_XTC}" -dt "${MD_LENGTH}" -novel -ndec 2
cp wrapped.gro "${TRJNAME}_sampled.gro"

# Extract fitted trajectory
printf "non-Water_&_!K_&_!CL\n" | gmx trjconv -f "${SAMPLED_XTC}" -s "${MDNAME}.tpr" \
  -o "${WHOLE_XTC}" -n index.ndx -pbc whole -novel -ndec 2
printf "System\nSystem\n" | gmx trjconv -f "${WHOLE_XTC}" -s "${STRIP_TPR}" \
  -o "${NOJUMP_XTC}" -n strip.ndx -pbc nojump -center -ndec 2
printf "System\nSystem\nSystem\n" | gmx trjconv -f "${NOJUMP_XTC}" -s "${STRIP_TPR}" \
  -o "${FIT_XTC}" -n strip.ndx -fit progressive -center -ndec 2
printf "Backbone\nSystem\n" | gmx trjconv -f "${FIT_XTC}" -s "${STRIP_TPR}" \
  -o "${NOROT_XTC}" -n strip.ndx -fit rot+trans -ndec 2
printf "System\n" | gmx trjconv -f "${NOROT_XTC}" -s "${STRIP_TPR}" \
  -o "${STRIP_XTC}" -n strip.ndx -ndec 2
printf "System\n" | gmx trjconv -f "${STRIP_XTC}" -s "${STRIP_TPR}" \
  -o "${TRJNAME}_strip.gro" -dump -1

rm -rf "${RUN_TMPDIR}" strip.ndx

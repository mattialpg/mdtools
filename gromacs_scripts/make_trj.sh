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
WHOLE_XTC="${RUN_TMPDIR}/${TRJNAME}_whole.xtc"
NOJUMP_XTC="${RUN_TMPDIR}/${TRJNAME}_nojump.xtc"
FIT_XTC="${RUN_TMPDIR}/${TRJNAME}_fit.xtc"
NOROT_XTC="${RUN_TMPDIR}/${TRJNAME}_norot.xtc"
STRIP_XTC="${TRJNAME}_strip.xtc"

if [[ ! -s "${MDNAME}.gro" ]]; then
  echo "0" | gmx trjconv -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" -o "${MDNAME}.gro" -dump -1
fi

# Extract fitted trajectory
printf "System\n" | gmx trjconv -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" -o "${WHOLE_XTC}" \
  -n index.ndx -pbc whole -dt "${MD_LENGTH}" -novel -ndec 2
printf "System\nSystem\n" | gmx trjconv -f "${WHOLE_XTC}" -s "${MDNAME}.tpr" -o "${NOJUMP_XTC}" \
  -n index.ndx -pbc nojump -center -ndec 2
printf "System\nSystem\nSystem\n" | gmx trjconv -f "${NOJUMP_XTC}" -s "${MDNAME}.tpr" \
  -o "${FIT_XTC}" -n index.ndx -fit progressive -center -ndec 2
printf "Backbone\nSystem\n" | gmx trjconv -s "${MDNAME}.tpr" -f "${FIT_XTC}" -o "${NOROT_XTC}" \
  -n index.ndx -fit rot+trans -ndec 2

# Strip water and ions
echo -e '!"Water_and_ions"\nq'| gmx make_ndx -f ${MDNAME}.gro -o index.ndx
echo '!Water_and_ions' | gmx trjconv -f "${NOROT_XTC}" -s "${MDNAME}.tpr" \
  -n index.ndx -dump -1 -o "${TRJNAME}_strip.gro"
echo '!Water_and_ions' | gmx trjconv -f "${NOROT_XTC}" -s "${MDNAME}.tpr" \
  -n index.ndx -o "${STRIP_XTC}" -ndec 2
rm -rf "${RUN_TMPDIR}"

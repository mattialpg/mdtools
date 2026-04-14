#!/bin/bash
set -euo pipefail

# Usage:
#   ./make_extend.sh -mdname MDNAME -extend NS
# Example:
#   ./make_extend.sh -mdname md_100 -extend 200

while [[ $# -gt 0 ]]; do
  case "$1" in
    -mdname) MDNAME="$2"; shift 2 ;;
    -extend) EXTEND_NS="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 -mdname MDNAME -extend NS"
      echo "Example: $0 -mdname md_100 -extend 200"
      exit 0
      ;;
    *) echo "Usage: $0 -mdname MDNAME -extend NS"; exit 1 ;;
  esac
done

if [[ -z "${MDNAME:-}" || -z "${EXTEND_NS:-}" ]]; then
  echo "Usage: $0 -mdname MDNAME -extend NS"
  exit 1
fi

log() {
  printf "\n\033[38;2;255;255;255;48;2;15;88;157m%s\033[0m\n\n" "$1"
}

### NOTIFICATION SYSTEM ###
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NOTIFY_SCRIPT="${SCRIPT_DIR}/notify.sh"
NOTIFY_TEXT=""
trap 'status=$?; if [[ -n "${NOTIFY_TEXT:-}" ]]; then bash "$NOTIFY_SCRIPT" "$status" "$(basename "$PWD"): ${NOTIFY_TEXT}" || true; else bash "$NOTIFY_SCRIPT" "$status" || true; fi' EXIT

# Extract original length
ORIG_NS="${MDNAME##*_}"
TOTAL_NS="$(awk -v A="$ORIG_NS" -v B="$EXTEND_NS" 'BEGIN { printf "%.0f", (A+B) }')"
NEW_MDNAME="md_${TOTAL_NS}"
EXTEND_PS="$(awk -v X="$EXTEND_NS" 'BEGIN { printf "%.0f", (X*1000.0) }')"
NOTIFY_TEXT="Extended ${MDNAME} by ${EXTEND_NS} ns -> ${NEW_MDNAME}"

# Extend dynamics
log " > Continuing simulation (appending to existing files)..."
gmx convert-tpr -s "${MDNAME}.tpr" -extend "${EXTEND_PS}" -o "${NEW_MDNAME}.tpr"
if [[ -f "${MDNAME}_pullx.xvg" ]]; then
  gmx mdrun -deffnm "${MDNAME}" -s "${NEW_MDNAME}.tpr" -cpi "${MDNAME}.cpt" \
    -px "${MDNAME}_pullx.xvg" -pf "${MDNAME}_pullf.xvg" \
    -bonded gpu -pme gpu -append
else
  gmx mdrun -deffnm "${MDNAME}" -s "${NEW_MDNAME}.tpr" -cpi "${MDNAME}.cpt" -bonded gpu -pme gpu -append
fi

# Rename all original output files except the new TPR
rename -f "s/^${MDNAME}/${NEW_MDNAME}/" *

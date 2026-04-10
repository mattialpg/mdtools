#!/bin/bash
set -euo pipefail

# Parse command-line args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) MDNAME="$2"; shift 2 ;;
    -o) TRJNAME="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 -i MDNAME [-o TRJNAME]"
      exit 1
      ;;
    *) echo "Usage: $0 -i MDNAME [-o TRJNAME]"; exit 1 ;;
  esac
done

if [[ -z "${MDNAME:-}" ]]; then
  echo "Error: missing -i"
  echo "Usage: $0 -i MDNAME [-o TRJNAME]"
  exit 1
fi

if [[ -z "${TRJNAME:-}" ]]; then
  TRJNAME="${MDNAME//md/trj}"
fi

LENGTH_TAG="${MDNAME#md_}"
if [[ "${LENGTH_TAG}" =~ ^([0-9]+([.][0-9]+)?) ]]; then
  LENGTH_NS="${BASH_REMATCH[1]}"
else
  echo "Error: could not parse numeric length from MDNAME='${MDNAME}'"
  echo "Expected names like md_100 or md_100_rep1"
  exit 1
fi

if [[ ! -s "${MDNAME}.gro" ]]; then
  echo "0" | gmx trjconv -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" -o "${MDNAME}.gro" -dump -1
fi

log() {
  printf "\n\033[38;2;255;255;255;48;2;15;88;157m%s\033[0m\n\n" "$1"
}

# Use temporary storage for intermediate trajectory files
WORK_TMPDIR="${TMPDIR:-/tmp}"
RUN_TMPDIR="$(mktemp -d "${WORK_TMPDIR%/}/make_trj.XXXXXX")"
WHOLE_XTC="${RUN_TMPDIR}/${TRJNAME}_whole.xtc"
NOJUMP_XTC="${RUN_TMPDIR}/${TRJNAME}_nojump.xtc"
FIT_XTC="${RUN_TMPDIR}/${TRJNAME}_fit.xtc"
NOROT_XTC="${RUN_TMPDIR}/${TRJNAME}_norot.xtc"
STRIP_XTC="${TRJNAME}_strip.xtc"

# # Extract fitted trajectory
# printf "System\n" | gmx trjconv -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" -o "${WHOLE_XTC}" \
#   -n index.ndx -pbc whole -dt "${LENGTH_NS}" -novel -ndec 2
# printf "System\nSystem\n" | gmx trjconv -f "${WHOLE_XTC}" -s "${MDNAME}.tpr" -o "${NOJUMP_XTC}" \
#   -n index.ndx -pbc nojump -center -ndec 2
# printf "System\nSystem\nSystem\n" | gmx trjconv -f "${NOJUMP_XTC}" -s "${MDNAME}.tpr" \
#   -o "${FIT_XTC}" -n index.ndx -fit progressive -center -ndec 2
# printf "Backbone\nSystem\n" | gmx trjconv -s "${MDNAME}.tpr" -f "${FIT_XTC}" -o "${NOROT_XTC}" \
#   -n index.ndx -fit rot+trans -ndec 2

# # Strip water and ions
# echo -e '!"Water_and_ions"\nq'| gmx make_ndx -f ${MDNAME}.gro -o index.ndx
# echo "!Water_and_ions" | gmx trjconv -f ${MDNAME}.gro -s "${MDNAME}.tpr" -o ${TRJNAME}_strip.gro -n index.ndx
# echo '!Water_and_ions' | gmx trjconv -f "${NOROT_XTC}" -s "${MDNAME}.tpr" \
#   -n index.ndx -o "${STRIP_XTC}" -ndec 2
# rm -rf "${RUN_TMPDIR}"

# # Calculate protein statistics
# printf "3 3" | gmx rms -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" -o trj_rmsd_prot.xvg -tu ns
# printf "3" | gmx rmsf -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" -o trj_rmsf_prot.xvg -res
# gmx dssp -f ${MDNAME}.xtc -s ${MDNAME}.tpr -o trj_dssp.dat -tu ns -b $((length - 10)) -e ${LENGTH}

# Calculate ligand RMSD and pocket-ligand distance on full trajectory
LIG="LIG1"
echo -e "q" | gmx make_ndx -f "${MDNAME}.gro" -o index.ndx
printf "${LIG}\n" | gmx rms -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" \
  -n index.ndx -o trj_rmsd_lig.xvg -tu ns -fit none
gmx select -s "${MDNAME}.tpr" -n index.ndx -on pocket.ndx \
  -select "group \"Protein\" and same residue as (within 0.45 of group \"${LIG}\"); group \"${LIG}\""
sed -i '0,/^\[.*\]$/s//[ Pocket ]/' pocket.ndx
gmx pairdist -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" -n pocket.ndx \
  -o trj_dist_lig.xvg -tu ns -ref 'group "Pocket"' -sel "group \"${LIG}\"" -type min

python ~/mdtools/interaction_analysis.py --trj_name "${TRJNAME}_strip"

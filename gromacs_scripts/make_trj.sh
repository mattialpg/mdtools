#!/bin/bash
set -euo pipefail

# Usage:
#   ./postprocess.sh -i md_NAME [-o trj_name]

# Parse command-line args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) MDNAME="$2"; shift 2 ;;
    -o) TRJNAME="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 -i md_NAME [-o trj_name]"
      exit 1
      ;;
    *) echo "Usage: $0 -i md_NAME [-o trj_name]"; exit 1 ;;
  esac
done

if [[ -z "${MDNAME:-}" ]]; then
  echo "Error: missing -i"
  echo "Usage: $0 -i md_NAME [-o trj_name]"
  exit 1
fi

if [[ -z "${TRJNAME:-}" ]]; then
  TRJNAME="${MDNAME//md/trj}"
fi

length=${MDNAME#md_}

log() {
  printf "\n\033[38;2;255;255;255;48;2;15;88;157m%s\033[0m\n\n" "$1"
}

# Check if index already exists
if ! grep -q "!Water_and_ions" index.ndx; then
	echo -e '!"Water_and_ions"\nq'| gmx make_ndx -f ${MDNAME}.gro -o index.ndx 
fi

# Strip water and ions
echo "!Water_and_ions" | gmx trjconv -f system.gro -s system.gro -o ${TRJNAME}_strip.gro -n index.ndx

# Extract trajectory again
echo "!Water_and_ions" | gmx trjconv -f ${MDNAME}.xtc -s system.gro -o ${TRJNAME}_whole.xtc -n index.ndx -pbc whole -dt ${length} -novel -ndec 2
echo "!Water_and_ions" "!Water_and_ions" | gmx trjconv -f ${TRJNAME}_whole.xtc -s system.gro -o ${TRJNAME}_nojump.xtc -n index.ndx -pbc nojump -center -ndec 2
echo "!Water_and_ions" "!Water_and_ions" "!Water_and_ions" | gmx trjconv -f ${TRJNAME}_nojump.xtc -s system.gro -o ${TRJNAME}_fit.xtc -n index.ndx -fit progressive -center -ndec 2 
echo "Backbone" "!Water_and_ions" | gmx trjconv -s system.gro -f ${TRJNAME}_fit.xtc -o ${TRJNAME}_strip.xtc -n index.ndx -fit rot+trans -ndec 2
rm ${TRJNAME}_whole.xtc ${TRJNAME}_nojump.xtc ${TRJNAME}_fit.xtc

printf "3 3" | gmx rms -f ${MDNAME}.xtc -s ${MDNAME}.tpr -o trj_rmsd_prot.xvg -tu ns
printf "3" | gmx rmsf -f ${MDNAME}.xtc -s ${MDNAME}.tpr -o trj_rmsf_prot.xvg -res
gmx dssp -f ${MDNAME}.xtc -s ${MDNAME}.tpr -o trj_dssp.dat -tu ns -b $((length - 10)) -e ${length}
printf "LIG LIG" | gmx rms -f ${MDNAME}.xtc -s ${MDNAME}.tpr -o trj_rmsd_lig.xvg -tu ns

python ~/mdtools/interaction_analysis.py --trj_name ${MDNAME}
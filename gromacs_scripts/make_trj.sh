#!/bin/bash
set -e

mdname=""
trjname=""

while getopts ":i:o:" opt; do
  case ${opt} in
    i ) mdname=$OPTARG ;;
    o ) trjname=$OPTARG ;;
    \? )
      echo "Usage: $0 -i md_NAME [-o output_name]"
      exit 1
      ;;
  esac
done

if [ -z "$mdname" ]; then
  echo "Error: missing -i"
  echo "Usage: $0 -i md_NAME [-o output_name]"
  exit 1
fi

if [ -z "$trjname" ]; then
  trjname="${mdname//md/trj}"
fi

length=${mdname#md_}

# Check if index already exists
if ! grep -q "!Water_and_ions" index.ndx; then
	echo -e '!"Water_and_ions"\nq'| gmx make_ndx -f ${mdname}.gro -o index.ndx 
fi

# Strip water and ions
echo "!Water_and_ions" | gmx trjconv -f system.gro -s system.gro -o ${trjname}_strip.gro -n index.ndx

# Extract trajectory again
echo "!Water_and_ions" | gmx trjconv -f ${mdname}.xtc -s system.gro -o ${trjname}_whole.xtc -n index.ndx -pbc whole -dt ${length} -novel -ndec 2
echo "!Water_and_ions" "!Water_and_ions" | gmx trjconv -f ${trjname}_whole.xtc -s system.gro -o ${trjname}_nojump.xtc -n index.ndx -pbc nojump -center -ndec 2
echo "!Water_and_ions" "!Water_and_ions" "!Water_and_ions" | gmx trjconv -f ${trjname}_nojump.xtc -s system.gro -o ${trjname}_fit.xtc -n index.ndx -fit progressive -center -ndec 2 
echo "Backbone" "!Water_and_ions" | gmx trjconv -s system.gro -f ${trjname}_fit.xtc -o ${trjname}_strip.xtc -n index.ndx -fit rot+trans -ndec 2
rm ${trjname}_whole.xtc ${trjname}_nojump.xtc ${trjname}_fit.xtc

printf "3 3" | gmx rms -f ${mdname}.xtc -s ${mdname}.tpr -o trj_rmsd_prot.xvg -tu ns
printf "3" | gmx rmsf -f ${mdname}.xtc -s ${mdname}.tpr -o trj_rmsf_prot.xvg -res
gmx dssp -f ${mdname}.xtc -s ${mdname}.tpr -o trj_dssp.dat -tu ns -b $((length - 10)) -e ${length}
printf "LIG LIG" | gmx rms -f ${mdname}.xtc -s ${mdname}.tpr -o trj_rmsd_lig.xvg -tu ns

python ~/mdtools/interaction_analysis.py --trj_name ${mdname}

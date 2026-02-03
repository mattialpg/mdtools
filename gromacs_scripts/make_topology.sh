#!/bin/bash
set -euo pipefail
shopt -s nullglob

# Check if an input file is provided
if [ -z "${1:-}" ] || [ "${1##*.}" != "mol2" ]; then
  echo "Usage: $0 <input_file.mol2>"
  exit 1
fi

input_file="$1"
base_name="${input_file%.*}"
output_name="${input_file%.*}_amber"
ligand_name=$(sed -n '2p' "${input_file}" | tr -d '\r' | xargs)

cat > tleap_topology.in <<EOF
source leaprc.gaff2
loadamberparams ${output_name}.frcmod
LIG = loadmol2 ${output_name}.mol2
check LIG
saveamberparm LIG ${output_name}.prmtop ${output_name}.rst7
quit
EOF

printf "\e[38;2;255;255;255m\e[48;2;15;88;157m\n >  Calculating charges...\e[0m\n\n"
antechamber -i "$input_file" -fi mol2 -o "${output_name}.mol2" -fo mol2 -c bcc -s 2 -at gaff2 -nc 0 -m 1
rm -rf "${base_name}.antechamber" && mkdir "${base_name}.antechamber"
mv ANTECHAMBER_* ATOMTYPE* sqm.* "${base_name}.antechamber"

printf "\e[38;2;255;255;255m\e[48;2;15;88;157m\n >  Generating Amber topology...\e[0m\n\n"
parmchk2 -i "${output_name}.mol2" -f mol2 -o "${output_name}.frcmod" -s gaff2
tleap -f tleap_topology.in

printf "\e[38;2;255;255;255m\e[48;2;15;88;157m\n >  Converting to GROMACS topology...\e[0m\n\n"
acpype -i "${output_name}.mol2" 2> >(grep -Ev \
  "OpenBabel|Open Babel|Cannot perform atom type translation|This Mol2 file is non-standard|Cannot interpret atom type" >&2)
mv -f "${output_name}.acpype" "${base_name}.acpype"
rename -f 's/_amber//' "${base_name}.acpype"/*

awk '/\[ atomtypes \]/,/\[ moleculetype \]/ {if ($0 !~ /\[ moleculetype \]/) print}' "${base_name}.acpype/${base_name}_GMX.itp" > "${base_name}.prm"
awk '/\[ moleculetype \]/,/\[ system \]/ {if ($0 !~ /\[ system \]/) print}' "${base_name}.acpype/${base_name}_GMX.itp" > "${base_name}.itp"
sed -i "s/\b${output_name}\b/${ligand_name}/g" "${base_name}.itp"

mv -f tleap_topology.in "${output_name}".* leap.log "${base_name}.acpype/"
cp "${base_name}.acpype/${base_name}_GMX.gro" "${base_name}.gro"

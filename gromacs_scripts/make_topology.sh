#!/bin/bash
set -euo pipefail
shopt -s nullglob

# Check ligand file
if [ -z "${1:-}" ] || [[ "$1" == *.* ]] || [ ! -s "${1}.sdf" ]; then
  echo "Usage: $0 <ligand_name>"
  exit 1
fi

ligand_name="$1"
output_name="${ligand_name}_amber"
ligand_md_id=$(sed -n '2p' "${ligand_name}.sdf" | tr -d '\r' | xargs)

cat > tleap_topology.in <<EOF
source leaprc.gaff2
loadamberparams ${output_name}.frcmod
LIGAND = loadmol2 ${output_name}.mol2
check LIGAND
saveamberparm LIGAND ${output_name}.prmtop ${output_name}.rst7
quit
EOF

printf "\e[38;2;255;255;255m\e[48;2;15;88;157m\n >  Calculating charges...\e[0m\n\n"
antechamber -i "${ligand_name}.sdf" -fi sdf -o "${output_name}.mol2" -fo mol2 -c bcc -s 2 -at gaff2 -nc 0 -m 1
sed -i "s/\<MOL\>/${ligand_md_id}/g" "${output_name}.mol2"
rm -rf "${ligand_name}.antechamber"; mkdir "${ligand_name}.antechamber"
mv ANTECHAMBER_* ATOMTYPE* sqm.* "${ligand_name}.antechamber"

printf "\e[38;2;255;255;255m\e[48;2;15;88;157m\n >  Generating Amber topology...\e[0m\n\n"
parmchk2 -i "${output_name}.mol2" -f mol2 -o "${output_name}.frcmod" -s gaff2
tleap -f tleap_topology.in

printf "\e[38;2;255;255;255m\e[48;2;15;88;157m\n >  Converting to GROMACS topology...\e[0m\n\n"
acpype -i "${output_name}.mol2" 2> >(grep -Ev \
  "OpenBabel|Open Babel|Cannot perform atom type translation|
  This Mol2 file is non-standard|Cannot interpret atom type|
  ==============================" >&2)
rm -rf "${ligand_name}.acpype"; mv -f "${output_name}.acpype" "${ligand_name}.acpype"
rename -f 's/_amber//' "${ligand_name}.acpype"/*

awk '/\[ atomtypes \]/,/\[ moleculetype \]/ {if ($0 !~ /\[ moleculetype \]/) print}' "${ligand_name}.acpype/${ligand_name}_GMX.itp" > "${ligand_name}.prm"
awk '/\[ moleculetype \]/,/\[ system \]/ {if ($0 !~ /\[ system \]/) print}' "${ligand_name}.acpype/${ligand_name}_GMX.itp" > "${ligand_name}.itp"
sed -i "s/\b${output_name}\b/${ligand_md_id}/g" "${ligand_name}.itp"

mv -f tleap_topology.in "${output_name}".* leap.log "${ligand_name}.acpype/"
cp "${ligand_name}.acpype/${ligand_name}_GMX.gro" "${ligand_name}.gro"

obabel "${ligand_name}.sdf" -O "${ligand_name}.pdb" --writeconect
sed -i "s/\<UNNAMED\>/${ligand_md_id}/g" "${ligand_name}.pdb"
sed -i "s/\<UNL\>/${ligand_md_id}/g" "${ligand_name}.pdb"

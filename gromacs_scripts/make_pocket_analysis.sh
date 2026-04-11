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


MDNAME="md_200"
TRJNAME="trj_200"
STRIP_XTC="${TRJNAME}_strip.xtc"
STRIP_GRO="${TRJNAME}_strip.gro"
NOLIG_XTC="${TRJNAME}_nolig.xtc"
NOLIG_GRO="${TRJNAME}_nolig.gro"


printf "q\n" | gmx make_ndx -f "${STRIP_GRO}" -o "strip.ndx"
printf "Protein\n" | gmx trjconv -s "${STRIP_GRO}" -f "${STRIP_GRO}" \
  -n "strip.ndx" -o "${NOLIG_GRO}"
printf "Protein\n" | gmx trjconv -s "${STRIP_GRO}" -f "${STRIP_XTC}" \
  -n "strip.ndx" -o "${NOLIG_XTC}" -ndec 2

# Calculate protein statistics
printf "q\n" | gmx make_ndx -f "${NOLIG_GRO}" -o "nolig.ndx"
printf "Protein\nProtein\n" | gmx rms -f "${NOLIG_XTC}" -s "${NOLIG_GRO}" \
  -n "nolig.ndx" -o trj_rmsd_prot.xvg -tu ns
printf "Protein\n" | gmx rmsf -f "${NOLIG_XTC}" -s "${NOLIG_GRO}" \
  -n "nolig.ndx" -o trj_rmsf_prot.xvg -res
gmx dssp -f "${MDNAME}.xtc" -s "${MDNAME}.tpr" -o trj_dssp.dat -tu ns
printf "Protein\n" | gmx gyrate -s "${NOLIG_GRO}" -f "${NOLIG_XTC}" \
  -n "nolig.ndx" -o trj_gyrate_prot.xvg
# printf "Protein\n" | gmx sasa -s "${NOLIG_GRO}" -f "${NOLIG_XTC}" \
#   -n "nolig.ndx" -o trj_sasa_prot.xvg -or trj_sasa_residues.xvg

# Calculate pocket statistics
printf "Protein\nProtein\n" | gmx rms -s "${NOLIG_GRO}" -f "${NOLIG_XTC}" \
  -n "nolig.ndx" -o trj_rmsd_pocket.xvg -tu ns
printf "Protein\n" | gmx rmsf -s "${NOLIG_GRO}" -f "${NOLIG_XTC}" \
  -n "nolig.ndx" -o trj_rmsf_pocket.xvg -res

# Define pocket and expanded pocket
gmx select -s "${STRIP_GRO}" -n "strip.ndx" -on pocket.ndx \
  -select "group \"Protein\" and same residue as (within 0.45 of group \"LIG\")"
sed -i '0,/^\[.*\]$/s//[ Pocket ]/' pocket.ndx
gmx select -s "${STRIP_GRO}" -n "strip.ndx" -on pocket.ndx \
  -select "group \"Protein\" and same residue as (within 1.0 of group \"LIG\")"
sed -i '0,/^\[.*\]$/s//[ Exp_Pocket ]/' pocket.ndx

# Compute pocket centroid and max radius
echo "Pocket" | gmx trjconv -f "${NOLIG_GRO}" -s "${NOLIG_GRO}" -n pocket.ndx -o pocket.gro
read -r POCKET_X POCKET_Y POCKET_Z POCKET_R <<EOF
$(awk '
  NR > 2 {
    x = substr($0, 21, 8) + 0
    y = substr($0, 29, 8) + 0
    z = substr($0, 37, 8) + 0
    n++
    sx += x; sy += y; sz += z
    X[n] = x; Y[n] = y; Z[n] = z
  }
  END {
    if (n == 0) {
      print "0 0 0 10"
      exit 0
    }
    cx = sx / n; cy = sy / n; cz = sz / n; r = 0
    for (i = 1; i <= n; i++) {
      dx = X[i] - cx; dy = Y[i] - cy; dz = Z[i] - cz
      d = sqrt(dx*dx + dy*dy + dz*dz)
      if (d > r) r = d
    }
    # Convert nm->A and add 2 A margin to ensure pocket coverage.
    printf("%.2f %.2f %.2f %.2f\n", 10*cx, 10*cy, 10*cz, 10*r + 2.0)
  }
' pocket.gro)
EOF
echo "POVME DefineSphere: ${POCKET_X} ${POCKET_Y} ${POCKET_Z} ${POCKET_R}"

# Exctract frames on expanded pocket
mkdir -p protein.povme
printf "Exp_Pocket\n" | gmx trjconv -s "${NOLIG_GRO}" -f "${NOLIG_XTC}" -n pocket.ndx -o protein.povme/frames.pdb

# Run POVME
POVME_CONFIG="protein.povme/povme_config.ini"
cat > "$POVME_CONFIG" << EOF
PDBFileName protein.povme/frames.pdb
OutputFilenamePrefix protein.povme/povme
GridSpacing 1.0
DistanceCutoff 1.09
PointsInclusionSphere ${POCKET_X} ${POCKET_Y} ${POCKET_Z} ${POCKET_R}
SavePocketVolumesTrajectory true
SaveTabbedVolumeFile true
CompressOutput false
EOF
python ~/POVME/POVME2.py "$POVME_CONFIG"

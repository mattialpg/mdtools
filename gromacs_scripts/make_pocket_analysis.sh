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
mapfile -t BSITE < <(yq -r '(.BSITE // .binding_site // [])[]' "${CONFIG_YAML}")
if [[ ${#BSITE[@]} -eq 0 ]]; then
  echo "Error: set BSITE or binding_site list in ${CONFIG_YAML}" >&2
  exit 1
fi
ANCHOR_RESID="${BSITE[0]//[^0-9]/}"

MDNAME="md_${MD_LENGTH}"
TRJNAME="trj_${MD_LENGTH}"
STRIP_XTC="${TRJNAME}_strip.xtc"
STRIP_GRO="${TRJNAME}_strip.gro"
MDPOCKET_DIR="protein.mdpocket"
PROT_XTC="${MDPOCKET_DIR}/${TRJNAME}_prot.xtc"
PROT_PDB="${MDPOCKET_DIR}/${TRJNAME}_prot.pdb"
POVME_DIR="protein.povme"

# # Define binding site
# printf "q\n" | gmx make_ndx -f "${STRIP_GRO}" -o strip.ndx
# gmx select -s "${STRIP_GRO}" -n strip.ndx -on bsite.ndx \
#   -select "group \"Protein\" and same residue as (within 0.45 of resid ${ANCHOR_RESID})"
# sed -i '0,/^\[.*\]$/s//[ Binding-Site ]/' bsite.ndx

# # Calculate binding-site statistics
# printf "Binding-Site\nBinding-Site\n" | gmx rms -s "${STRIP_GRO}" -f "${STRIP_XTC}" \
#   -n "bsite.ndx" -o trj_rmsd_bsite.xvg -tu ns
# printf "Binding-Site\n" | gmx rmsf -s "${STRIP_GRO}" -f "${STRIP_XTC}" \
#   -n "bsite.ndx" -o trj_rmsf_bsite.xvg -res

# # -----------------------------------------------------

# # Extract protein-only trajectory
# mkdir -p "${MDPOCKET_DIR}"
# printf "Protein\n" | gmx trjconv -s "${STRIP_GRO}" -f "${STRIP_GRO}" \
#   -n strip.ndx -o "${PROT_PDB}"
# printf "Protein\n" | gmx trjconv -s "${STRIP_GRO}" -f "${STRIP_XTC}" \
#   -n strip.ndx -o "${PROT_XTC}" -ndec 2

# mdpocket --trajectory_file "${PROT_XTC}" --trajectory_format xtc \
#   -f "${PROT_PDB}" --output_prefix "${MDPOCKET_DIR}/mdpocket"

# # Define pocket
# mkdir -p "${POVME_DIR}"
# gmx select -s "${STRIP_GRO}" -n strip.ndx -on "${POVME_DIR}/pocket.ndx" \
#   -select "group \"Protein\" and same residue as (within 1.0 of resid ${ANCHOR_RESID})"
# sed -i '0,/^\[.*\]$/s//[ Pocket ]/' "${POVME_DIR}/pocket.ndx"

# # Extract pocket frames
# printf "Pocket\n" | gmx trjconv -s "${STRIP_GRO}" -f "${STRIP_XTC}" \
#   -o "${POVME_DIR}/frames.pdb" -n "${POVME_DIR}/pocket.ndx"

# # Compute pocket centroid and max radius
# echo "Pocket" | gmx trjconv -f "${STRIP_GRO}" -s "${STRIP_GRO}" \
#  -o "${POVME_DIR}/pocket.gro" -n "${POVME_DIR}/pocket.ndx" 
# read -r POCKET_X POCKET_Y POCKET_Z POCKET_R <<EOF
# $(awk '
#   NR == 2 {
#     n_atoms = int($1)
#     last_atom_line = 2 + n_atoms
#   }
#   NR > 2 && NR <= last_atom_line {
#     x = substr($0, 21, 8) + 0
#     y = substr($0, 29, 8) + 0
#     z = substr($0, 37, 8) + 0
#     n++
#     sx += x; sy += y; sz += z
#     X[n] = x; Y[n] = y; Z[n] = z
#   }
#   END {
#     if (n == 0) {
#       print "0 0 0 10"
#       exit 0
#     }
#     cx = sx / n; cy = sy / n; cz = sz / n; r = 0
#     for (i = 1; i <= n; i++) {
#       dx = X[i] - cx; dy = Y[i] - cy; dz = Z[i] - cz
#       d = sqrt(dx*dx + dy*dy + dz*dz)
#       if (d > r) r = d
#     }
#     # Convert nm->A and add 2 A margin to ensure pocket coverage.
#     printf("%.2f %.2f %.2f %.2f\n", 10*cx, 10*cy, 10*cz, 10*r + 2.0)
#   }
# ' "${POVME_DIR}/pocket.gro")
# EOF
# echo "POVME DefineSphere: ${POCKET_X} ${POCKET_Y} ${POCKET_Z} ${POCKET_R}"

# # Calculate pocket volume
# POVME_CONFIG="${POVME_DIR}/povme_config.yml"
# cat > "$POVME_CONFIG" << EOF
# grid_spacing: 1.0
# distance_cutoff: 1.09
# points_inclusion_sphere:
#   - center: [${POCKET_X}, ${POCKET_Y}, ${POCKET_Z}]
#     radius: ${POCKET_R}
# save_pocket_volumes_trajectory: true
# save_individual_pocket_volumes: true
# save_volumetric_density_map: true
# compress_output: false
# EOF
# /home/mattia/.miniforge/envs/my-chem/bin/povme volume \
#   -c "$POVME_CONFIG" -i "${POVME_DIR}/frames.pdb" -o "${POVME_DIR}/"
# rm -f "${POVME_DIR}"/frame_*.pdb "${POVME_DIR}"/frame_*.pdb.gz

sed -i '/^REMARK Volume/d' "${POVME_DIR}/volume_trajectory.pdb"
sed -i 's/^REMARK Frame /MODEL        /' "${POVME_DIR}/volume_trajectory.pdb"
sed -i ':a;N;$!ba;s/^END$/TER\nENDMDL/gm' "${POVME_DIR}/volume_trajectory.pdb"

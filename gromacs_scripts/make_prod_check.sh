#!/bin/bash
set -euo pipefail

# Parse command-line args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -length) LENGTH="$2"; shift 2 ;;
    *) echo "Usage: $0 [-length NS]"; exit 1 ;;
  esac
done

DT=0.002     # ps (2 fs)
OUT_NAME="md_${LENGTH}"
NSTEPS="$(awk -v L="$LENGTH" -v DT="$DT" 'BEGIN { printf "%.0f", (L*1000.0/DT) }')"

log() {
  printf "\n\033[38;2;255;255;255;48;2;15;88;157m%s\033[0m\n\n" "$1"
}

# Generate MDP
cat > prod.mdp <<EOF
; Force field = ff99SB-ILDN

; Run parameters
integrator              = md
nsteps                  = ${NSTEPS}    ; ${DT} * ${NSTEPS} = ${LENGTH} ns
dt                      = ${DT}

; Standard output 
nstenergy               = 5000
nstlog                  = 5000
nstxout-compressed      = 5000

; Full precision trajectory (TRR)
nstxout             	= 20000
nstvout             	= 0
nstfout             	= 0

; Bond parameters
continuation            = yes
constraint-algorithm    = LINCS
constraints             = h-bonds
lincs-iter              = 1
lincs-order             = 4

; Neighbor searching
cutoff-scheme           = Verlet
nstlist                 = 20
pbc                     = xyz
rlist                   = 1.0

; Electrostatics
coulombtype             = PME
rcoulomb                = 1.0
pme-order               = 4
fourierspacing          = 0.125

; Van der Waals
vdwtype                 = cutoff
vdw-modifier            = Potential-shift
rvdw-switch             = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres

; Temperature coupling
tcoupl                  = V-rescale
tc-grps                 = System
tau-t                   = 1.0
ref-t                   = 300

; Pressure coupling
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau-p                   = 5.0
ref-p                   = 1.0
compressibility         = 4.5e-5

; Pull code for monitoring Protein–LIG COM distance (no bias)
pull                    = yes
pull-ngroups            = 2
pull-group1-name        = Protein
pull-group2-name        = LIG
pull-ncoords            = 1
pull-coord1-geometry    = distance
pull-coord1-groups      = 1 2
pull-coord1-dim         = Y Y Y
pull-coord1-start       = yes
pull-coord1-k 		= 0
pull-nstxout            = 200
EOF

### PRODUCTION ###
log " > Running production dynamics (${LENGTH} ns)..."
gmx grompp -f prod.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o "${OUT_NAME}.tpr" -maxwarn 100
gmx mdrun -deffnm "${OUT_NAME}" -bonded gpu -pme gpu & MDPID=$!

while kill -0 "${MDPID}" 2>/dev/null; do
  sleep 30
  
  dist=$(awk '$1 !~ /^[@#]/ && $2+0==$2 {!first && (first=$2); last=$2} END{d=last-first; print d<0?-d:d}' "${OUT_NAME}_pullx.xvg")
  log " > Protein-LIG distance: ${dist} nm"
done

wait "${MDPID}"


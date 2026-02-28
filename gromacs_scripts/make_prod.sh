#!/bin/bash
set -euo pipefail

# Usage:
#   ./run.sh -length NS

# Parse command-line args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -length) LENGTH="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 -length NS"
      exit 0
      ;;
    *) echo "Usage: $0 -length NS"; exit 1 ;;
  esac
done

LENGTH="${LENGTH:-100}"   # ns
DT=0.002                  # ps (2 fs)
OUT_NAME="md_${LENGTH}"
NSTEPS="$(awk -v L="$LENGTH" -v DT="$DT" 'BEGIN { printf "%.0f", (L*1000.0/DT) }')"

log() {
    printf "\033[38;2;255;255;255;48;2;15;88;157m%s\033[0m\n" "$1"
}

# Generate MDP
cat > prod.mdp <<EOF
; Force field = ff99SB-ILDN

; Run parameters
integrator              = md
nsteps                  = ${NSTEPS}    ; ${DT} * ${NSTEPS} = ${LENGTH} ns
dt                      = ${DT}

; Output control
nstenergy               = 5000
nstlog                  = 5000
nstxout-compressed      = 5000

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
EOF

### PRODUCTION ###
log " > Running production dynamics..."
echo "Simulation length: ${LENGTH} ns (${DT} ps * ${NSTEPS} steps)"
gmx grompp -f prod.mdp -c npt.gro -t npt.cpt -p topol.top -o "${OUT_NAME}.tpr" -maxwarn 100
gmx mdrun -deffnm "${OUT_NAME}" -bonded gpu -pme gpu
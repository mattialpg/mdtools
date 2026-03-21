#!/bin/bash
set -euo pipefail

# Parse command-line args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -len) LENGTH="$2"; shift 2 ;;
    *) echo "Usage: $0 [-len NS]"; exit 1 ;;
  esac
done

DT=0.002     # ps (2 fs)
OUT_NAME="md_${LENGTH}"
NSTEPS="$(awk -v L="$LENGTH" -v DT="$DT" 'BEGIN { printf "%.0f", (L*1000.0/DT) }')"
MONITOR_INTERVAL=300   # 5 minutes
RMSD_CUTOFF=1.5        # nm

log() {
  printf "\n\033[38;2;255;255;255;48;2;15;88;157m%s\033[0m\n\n" "$1"
}

# Generate MDP
cat > prod.mdp <<EOF2
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
nstxout              	  = 20000
nstvout              	  = 0
nstfout              	  = 0

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
EOF2

### PRODUCTION ###
log " > Running production dynamics (${LENGTH} ns)..."
gmx grompp -f prod.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o "${OUT_NAME}.tpr" -maxwarn 100
gmx mdrun -deffnm "${OUT_NAME}" -bonded gpu -pme gpu & MDPID=$!
killed_by_rmsd=0

while kill -0 "${MDPID}" 2>/dev/null; do
  sleep "${MONITOR_INTERVAL}"

  [[ -s "${OUT_NAME}.xtc" ]] || {
    log " > RMSD check waiting for ${OUT_NAME}.xtc to be written..."
    continue
  }

  # Compute ligand RMSD from current trajectory snapshot.
  if ! printf 'LIG\nLIG\n' | gmx rms \
    -s "${OUT_NAME}.tpr" \
    -f "${OUT_NAME}.xtc" \
    -n index.ndx \
    -o "${OUT_NAME}_lig_rmsd_tmp.xvg" \
    -quiet >/dev/null 2>&1; then
    log " > RMSD check skipped (gmx rms failed this cycle)."
    continue
  fi

  rmsd=$(awk '$1 !~ /^[@#]/ && $2+0==$2 {last=$2} END{if(last=="") print "NA"; else print last}' "${OUT_NAME}_lig_rmsd_tmp.xvg")
  if awk -v r="${rmsd}" -v c="${RMSD_CUTOFF}" 'BEGIN{exit !((r+0) > c)}'; then
    log " > LIG RMSD (${rmsd} nm) > ${RMSD_CUTOFF} nm: stopping dynamics."
    kill -TERM "${MDPID}" 2>/dev/null || true
    killed_by_rmsd=1
    break
  fi
done

if [[ "${killed_by_rmsd}" -eq 1 ]]; then
  wait "${MDPID}" || true
else
  wait "${MDPID}"
fi

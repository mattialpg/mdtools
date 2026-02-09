#!/bin/bash
set -euo pipefail

# Parse command-line args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -length) LENGTH="$2"; shift 2 ;;
    *) echo "Usage: $0 [-length NS]"; exit 1 ;;
  esac
done

LENGTH=100   # ns
DT=0.002     # ps (2 fs)
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

THRESH_NM=1.2        # stop if distance > this (nm)
POLL_SEC=300         # check every 10 minutes wall time
CPT_MIN=5            # write checkpoint every 5 minutes wall time

# safety: avoid monitor hanging forever inside a gmx tool
TIMEOUT_SEC=180      # max seconds for each gmx analysis command

gmx grompp -f prod.mdp -c npt.gro -t npt.cpt -p topol.top -o "${OUT_NAME}.tpr" -maxwarn 100

# run MD in background with frequent checkpoints
gmx mdrun -deffnm "${OUT_NAME}" -bonded gpu -pme gpu -cpt ${CPT_MIN} &
MDPID=$!

# monitor loop
while kill -0 "${MDPID}" 2>/dev/null; do
  sleep "${POLL_SEC}"

  [ -f "${OUT_NAME}.xtc" ] || continue
  [ -f "${OUT_NAME}.tpr" ] || continue

  # copy xtc to avoid racing with a file being written
  cp -f "${OUT_NAME}.xtc" monitor.xtc 2>/dev/null || continue

  # get time (ps) of last frame currently present in xtc
  t_last=$(gmx check -f monitor.xtc 2>/dev/null | awk '/Last frame/ {print $(NF-1)}')
  [ -n "${t_last}" ] || { log " > Could not read last frame time from monitor.xtc"; continue; }

  # extract last frame to a gro
  timeout "${TIMEOUT_SEC}" bash -lc 'printf "0\n" | gmx trjconv -s "'"${OUT_NAME}.tpr"'" -f monitor.xtc -dump "'"${t_last}"'" -o last.gro' \
    >/dev/null 2>trjconv.err || { log " > trjconv failed/hung: $(tail -5 trjconv.err | tr "\n" " ")"; continue; }

  # measure COM distance between Protein and LIG using an index file
  timeout "${TIMEOUT_SEC}" gmx distance -s "${OUT_NAME}.tpr" -f last.gro -n index.ndx \
    -select 'com of group "Protein" plus com of group "LIG"' \
    -oall dist.xvg >/dev/null 2>distance.err || { log " > distance failed/hung: $(tail -5 distance.err | tr "\n" " ")"; continue; }

  d=$(awk 'NF==2 && $1 !~ /^[@#]/ {v=$2} END{print v+0}' dist.xvg)

  log " > Protein-LIG distance: ${d} nm (last frame at ${t_last} ps)"

  awk -v d="${d}" -v t="${THRESH_NM}" 'BEGIN{exit !(d>t)}'
  if [ $? -eq 0 ]; then
    log " > Threshold exceeded (${d} > ${THRESH_NM}), stopping mdrun"
    kill -TERM "${MDPID}"
    break
  fi
done

wait "${MDPID}"
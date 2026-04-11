#!/bin/bash
set -euo pipefail

# Parse command-line args
EARLYSTOP=0
LENGTH=500
NREPS=1
BASE_SEED=1000
while [[ $# -gt 0 ]]; do
  case "$1" in
    -len) LENGTH="$2"; shift 2 ;;
    -nreps) NREPS="$2"; shift 2 ;;
    --earlystop) EARLYSTOP=1; shift ;;
    -h|--help)
      echo "Usage: $0 [-len NS] [-nreps N] [--earlystop]"
      echo "Defaults: -len 100 -nreps 5"
      exit 0
      ;;
    *) echo "Usage: $0 [-len NS] [-nreps N] [--earlystop]"; exit 1 ;;
  esac
done

DT=0.002     # ps (2 fs)
NSTEPS="$(awk -v L="$LENGTH" -v DT="$DT" 'BEGIN { printf "%.0f", (L*1000.0/DT) }')"
MONITOR_INTERVAL=300   # seconds
DIST_CUTOFF=1.5        # nm
EARLYSTOP_GRACE_SEC=600  # continue MD for 10 minutes after trigger

log() {
  printf "\n\033[38;2;255;255;255;48;2;15;88;157m%s\033[0m\n\n" "$1"
}

### NOTIFICATION SYSTEM ###
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NOTIFY_SCRIPT="${SCRIPT_DIR}/notify.sh"
NOTIFY_TEXT=""
trap 'status=$?; if [[ -n "${NOTIFY_TEXT:-}" ]]; then bash "$NOTIFY_SCRIPT" "$status" "$(basename "$PWD"): ${NOTIFY_TEXT}" || true; else bash "$NOTIFY_SCRIPT" "$status" || true; fi' EXIT

if [[ "${EARLYSTOP}" -eq 1 ]]; then
  # Build pocket.ndx
  echo -e "0\nq" | gmx make_ndx -f npt.gro -o pocket.ndx >/dev/null
  gmx select -s npt.gro -n pocket.ndx -on pocket_tmp.ndx -select 'group "Protein" and same residue as (within 0.45 of group "LIG")'
  sed -i '0,/^\[.*\]$/s//[ Pocket ]/' pocket_tmp.ndx
  cat pocket_tmp.ndx >> pocket.ndx
  rm -f pocket_tmp.ndx

  # Select first atom in Pocket group
  PULL_GROUP1_PBCATOM="$(awk 'NR>1 {print $1; exit}' pocket.ndx)"
  if [[ -z "${PULL_GROUP1_PBCATOM:-}" ]]; then
    echo "Error: Pocket group is empty; cannot set pull-group1-pbcatom." >&2
    exit 1
  fi
fi

for rep in $(seq 1 "${NREPS}"); do
  if [[ "${NREPS}" -eq 1 ]]; then
    OUT_NAME="md_${LENGTH}"
    MDP_NAME="prod.mdp"
  else
    OUT_NAME="md_${LENGTH}_rep${rep}"
    MDP_NAME="prod_rep${rep}.mdp"
  fi
  SEED=$((BASE_SEED + rep))

  # Generate MDP
  cat > "${MDP_NAME}" <<EOF
; Force field = ff99SB-ILDN

; Run parameters
integrator              = md
nsteps                  = ${NSTEPS}    ; ${LENGTH} ns
dt                      = ${DT}

; Standard output
nstenergy               = 5000
nstlog                  = 5000
nstxout-compressed      = 5000

; Bond parameters
continuation            = no
gen-vel                 = yes
gen-temp                = 300
gen-seed                = ${SEED}
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

  if [[ "${EARLYSTOP}" -eq 1 ]]; then
    cat >> "${MDP_NAME}" <<EOF
; Pull code to monitor Pocket-LIG COM distance (no bias)
pull                    = yes
pull-ngroups            = 2
pull-group1-name        = Pocket
pull-group2-name        = LIG
pull-ncoords            = 1
pull-coord1-geometry    = distance
pull-coord1-groups      = 1 2
pull-coord1-dim         = Y Y Y
pull-coord1-start       = yes
pull-coord1-k           = 0
pull-group1-pbcatom     = ${PULL_GROUP1_PBCATOM}
pull-pbc-ref-prev-step-com = yes
pull-nstxout            = 5000    ; 10 ps
EOF
  fi

  ### PRODUCTION ###
  log " > Running ${OUT_NAME} (${LENGTH} ns)..."
  if [[ "${EARLYSTOP}" -eq 1 ]]; then
    gmx grompp -f "${MDP_NAME}" -c npt.gro -p topol.top -n pocket.ndx -o "${OUT_NAME}.tpr" -maxwarn 100
  else
    gmx grompp -f "${MDP_NAME}" -c npt.gro -p topol.top -o "${OUT_NAME}.tpr" -maxwarn 100
  fi
  gmx mdrun -deffnm "${OUT_NAME}" -bonded auto -pme auto -update auto & MDPID=$!
  killed_by_distance=0
  grace_started=0
  trigger_time=0
  trigger_dist="NA"

  if [[ "${EARLYSTOP}" -eq 1 ]]; then
    while kill -0 "${MDPID}" 2>/dev/null; do
      sleep "${MONITOR_INTERVAL}"
      kill -0 "${MDPID}" 2>/dev/null || break

      [[ -s "${OUT_NAME}_pullx.xvg" ]] || continue

      dist=$(awk '$1 !~ /^[@#]/ && $1+0==$1 && $2+0==$2 && $1 > 100 { !first && (first=$2); last=$2 }
             END { if(first=="" || last=="") print "NA";
                   else { d=last-first; print d<0?-d:d } }' "${OUT_NAME}_pullx.xvg")
      echo " > Pocket-LIG distance drift: ${dist} nm"
      if awk -v d="${dist}" -v c="${DIST_CUTOFF}" 'BEGIN{exit !((d+0) > c)}'; then
        if [[ "${grace_started}" -eq 0 ]]; then
          grace_started=1
          trigger_time=$(date +%s)
          trigger_dist="${dist}"
        else
          now=$(date +%s)
          elapsed=$(( now - trigger_time ))
          if (( elapsed >= EARLYSTOP_GRACE_SEC )); then
            NOTIFY_TEXT="MD stopped by distance (initial drift = ${trigger_dist} nm > ${DIST_CUTOFF} nm)"
            echo " > ${NOTIFY_TEXT}"
            kill -TERM "${MDPID}" 2>/dev/null || true
            killed_by_distance=1
            break
          fi
        fi
      fi
    done
  fi

  if [[ "${killed_by_distance}" -eq 1 ]]; then
    wait "${MDPID}" || true
    echo "0" | gmx trjconv -f "${OUT_NAME}.xtc" -s "${OUT_NAME}.tpr" -o "${OUT_NAME}.gro" -dump -1
  else
    wait "${MDPID}"
  fi
done

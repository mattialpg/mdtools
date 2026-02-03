#!/bin/bash
set -euo pipefail

# Derived values
OUT_NAME="md_pull"

# Generate MDP
cat > pull.mdp <<EOF
; Working with ff99SB-ILDN

; Run parameters
integrator              = md
nsteps                  = 250000    ; 500 ps
dt                      = 0.002

; Output control
nstenergy               = 500
nstlog                  = 500
nstxout-compressed      = 500

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
pcoupl                  = Parrinello-Rahman
pcoupltype              = anisotropic
tau-p                   = 5.0
ref-p                   = 1 1 1 0 0 0
compressibility         = 4.5e-5 4.5e-5 0 0 0 0

; Pull code
pull                    = yes
pull_ncoords            = 1           ; only one reaction coordinate
pull_ngroups            = 2           ; two groups defining one reaction coordinate
pull_group1_name        = GROUP1
pull_group2_name        = GROUP2 
pull_coord1_type        = umbrella
pull_coord1_geometry    = direction-periodic 
pull_coord1_dim         = N N Y
pull_coord1_vec         = 0 0 1
pull_coord1_groups      = 1 2         ; groups 1 and 2 define the reaction coordinate
pull_coord1_start       = yes         ; define initial COM distance > 0
pull_coord1_rate        = 0.005       ; 0.01 nm per ps = 10 nm per ns
pull_coord1_k           = 1000        ; kJ mol^-1 nm^-2
pull-group1-pbcatom     = 1170 
pull-pbc-ref-prev-step-com    = yes
EOF

### PRODUCTION ###
printf "\e[48;2;224;64;251m\e[30m\n >  Running steered dynamics...\e[0m\n\n"
gmx grompp -f pull.mdp -c npt.gro -t npt.cpt -p topol.top -n pull.ndx -o "${OUT_NAME}.tpr" -maxwarn 100
gmx mdrun -deffnm "${OUT_NAME}" -bonded gpu -pme gpu

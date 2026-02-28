#!/bin/bash
set -euo pipefail

# Usage:
#   bash make_equil.sh [ligand1 ligand2 ...]

# Parse command-line args
LIGANDS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      echo "Usage: bash make_equil.sh [ligand1 ligand2 ...]"
      exit 0
      ;;
    *) LIGANDS+=("$1"); shift ;;
  esac
done

ligands="${LIGANDS[*]}"

log() {
  printf "\n\033[38;2;255;255;255;48;2;15;88;157m%s\033[0m\n\n" "$1"
}

log " >  Preparing protein..."
gmx pdb2gmx -f protein.pdb -o protein.gro -merge all -ignh -ff amber99sb-ildn -water tip3p
sleep 2

### CREATE COMPLEX ###
if [ -n "$ligands" ]; then
    head -n -1 protein.gro > complex.gro

    ### ADD LIGANDS TO .GRO FILE ###
    for ligand in $ligands; do
        tail -n +3 "$ligand.gro" | head -n -1 > ligand.tmp
        sed -i "2s/.*/$(( $(sed -n '2p' complex.gro) + $(wc -l < ligand.tmp) ))/" complex.gro
        cat ligand.tmp >> complex.gro
        rm ligand.tmp
    done
    tail -1 protein.gro >> complex.gro

    ### ADD LIGANDS TO TOPOLOGY ###
    line=$(grep -n "Include forcefield parameters" topol.top | cut -d: -f1)
    {
        echo ""
        echo "; Include ligand parameters"
        for ligand in $ligands; do
            echo "#include \"${ligand}.prm\""
        done
    } | sed -i "$((line + 1))r /dev/stdin" topol.top

    line=$(grep -n "Include topology for ions" topol.top | cut -d: -f1)
    {
        echo ""
        echo "; Include ligand topology"
            for ligand in $ligands; do
            echo "#include \"${ligand}.itp\""
        done
    } | sed -i "$((line + 1))r /dev/stdin" topol.top

    for ligand in $ligands; do
          ligand_md_id=$(sed -n '2p' "${ligand}.sdf" | tr -d '\r' | xargs)  # Extract the second line from ${ligand}.sdf
        sed -i -e "\$a${ligand_md_id}                 1" topol.top  # Append the ligand name at the end of topol.top
    done    
else
    cp protein.gro complex.gro
fi

### CREATE BOX AND SOLVATE ###
# Should be -d >= rvdw but d != 1.0 fails to minimise
log " >  Solvating complex..."
cat << EOF > ions.mdp
integrator       = steep           ; Algorithm (steep = steepest descent minimization)
emtol            = 1000.0          ; Stop when the max force < 100 kJ/mol/nm
emstep           = 0.01            ; Energy step size
nsteps           = 50000           ; Maximum number of (minimization) steps to perform

cutoff-scheme    = Verlet          ; Cut-off scheme
nstlist          = 20              ; Frequency to update the neighbor list and long range forces
rlist            = 1.0             ; Cut-off for making neighbor list (short range forces)
coulombtype      = PME             ; Treatment of long range electrostatic interactions
rcoulomb         = 1.0             ; Long range electrostatic cut-off
vdwtype          = cutoff
rvdw             = 1.0             ; Long range Van der Waals cut-off
DispCorr         = EnerPres
pbc              = xyz             ; Periodic Boundary Conditions
EOF

gmx editconf -f complex.gro -o complex_box.gro -bt dodecahedron -d 1
gmx solvate -cp complex_box.gro -cs spc216.gro -p topol.top -o complex_solv.gro
gmx grompp -f ions.mdp -c complex_solv.gro -p topol.top -o ions.tpr -maxwarn 100
echo 'SOL' | gmx genion -s ions.tpr -o system.gro -p topol.top -pname NA -nname CL -neutral
sleep 2
rm complex_box.gro complex_solv.gro

if [ -n "$ligands" ]; then
    ### CREATE LIGAND RESTRAINTS ###
    log " >  Creating ligand restraints..."
    line=$(sed -n '/Include ligand topology/=' topol.top)
    
    for ligand in $ligands; do
        echo -e "2 & a C*\ndel 0-2\nq" | gmx make_ndx -f "${ligand}.gro" -o "index_${ligand}.ndx"
        gmx genrestr -f "${ligand}.gro" -n "index_${ligand}.ndx" -o "posre_${ligand}.itp" -fc 1000 1000 1000
    done

    ### ADD RESTRAINTS TO TOPOLOGY ###
    line=$(sed -n '/Include ligand topology/=' topol.top)
    {
        echo ""
        echo "; Ligand position restraints"
        echo "#ifdef POSRES_LIG"
        for ligand in $ligands; do
            echo "#include \"posre_${ligand}.itp\""
        done
        echo "#endif"
    } | sed -i "$((line + 1))r /dev/stdin" topol.top
fi

### CREATE C-ALPHA RESTRAINTS ###
log " >  Creating protein restraints..."
echo -e "3\nq" | gmx make_ndx -f system.gro -o index.ndx
echo "3" | gmx genrestr -f system.gro -n index.ndx -o posre.itp -fc 1000 1000 1000
sleep 2

### MINIMISE COMPLEX ###
log " >  Minimising complex..."
cat << EOF > em.mdp
; Working with ff99SB-ILDN

integrator       = steep           ; Algorithm (steep = steepest descent minimization)
emtol            = 100.0           ; Stop minimization when the maximum force < 10.0 kJ/mol
emstep           = 0.001           ; Energy step size
nsteps           = 50000           ; Maximum number of (minimization) steps to perform

cutoff-scheme    = Verlet
nstlist          = 20              ; Frequency to update the neighbor list and long range forces
rlist            = 1.0             ; Cut-off for making neighbor list (short range forces)
coulombtype      = PME             ; Treatment of long range electrostatic interactions
rcoulomb         = 1.0             ; Long range electrostatic cut-off
vdwtype          = cutoff
rvdw             = 1.0             ; Long range Van der Waals cut-off
DispCorr         = EnerPres
pbc              = xyz             ; Periodic Boundary Conditions
EOF

gmx grompp -f em.mdp -c system.gro -r system.gro -p topol.top -o em.tpr -maxwarn 100
gmx mdrun -deffnm em -ntomp 24 -ntmpi 1 -v
sleep 2

### NVT/NPT EQUILIBRATION ###
log " >  Equilibrating in NVT/NPT ensemble..."
cat << EOF > nvt.mdp
; Working with ff99SB-ILDN

; Run control
define                  = -DPOSRES -DPOSRES_LIG
integrator              = md             ; leap-frog integrator
nsteps                  = 250000         ; 2 * 250000 = 500 ps
dt                      = 0.002          ; 2 fs

; Output control
nstenergy               = 500            ; save energies every 10 ps
nstlog                  = 500            ; update log file every 10 ps
nstxout-compressed      = 500            ; save coordinates every 10 ps

; Bond parameters
continuation            = no             ; first dynamics run
constraint-algorithm    = LINCS          ; holonomic constraints
constraints             = h-bonds        ; bonds to H are constrained 
lincs-iter              = 1              ; accuracy of LINCS
lincs-order             = 4              ; also related to accuracy

; Neighbor searching
cutoff-scheme           = Verlet
nstlist                 = 20
pbc                     = xyz
rlist                   = 1.0

; Electrostatics
coulombtype             = PME            ; particle-Mesh Ewald electrostatics
rcoulomb                = 1.0            ; short-range electrostatic cutoff (nm)
pme-order               = 4              ; cubic interpolation
fourierspacing          = 0.125          ; grid spacing for FFT

; Van der Waals
vdwtype                 = cutoff
vdw-modifier            = Potential-shift
rvdw                    = 1.0            ; short-range van der Waals cutoff (nm)
DispCorr                = EnerPres

; Temperature coupling
tcoupl                  = V-rescale      ; thermostat
tc-grps                 = System
tau-t                   = 1.0            ; time constant (ps)
ref-t                   = 300            ; reference temperature (K)

; Pressure coupling
pcoupl                  = no             ; no pressure coupling in NVT

; Velocity generation
gen-vel                 = yes            ; assign velocities from Maxwell distribution
gen-temp                = 300            ; temperature for Maxwell distribution
gen-seed                = -1             ; generate a random seed
EOF

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 100
gmx mdrun -deffnm nvt -bonded gpu -pme gpu 

cat << EOF > npt.mdp
; Working with ff99SB-ILDN

; Run control
define                  = -DPOSRES -DPOSRES_LIG
integrator              = md             ; leap-frog integrator
nsteps                  = 250000         ; 2 * 250000 = 500 ps
dt                      = 0.002          ; 2 fs

; Output control
nstenergy               = 500            ; save energies every 10 ps
nstlog                  = 500            ; update log file every 10 ps
nstxout-compressed      = 500            ; save coordinates every 10 ps

; Bond parameters
continuation            = yes            ; continuing from NVT
constraint-algorithm    = LINCS          ; holonomic constraints
constraints             = h-bonds        ; bonds to H are constrained 
lincs-iter              = 1              ; accuracy of LINCS
lincs-order             = 4              ; also related to accuracy

; Neighbor searching
cutoff-scheme           = Verlet
nstlist                 = 20
pbc                     = xyz
rlist                   = 1.0

; Electrostatics
coulombtype             = PME            ; particle-Mesh Ewald electrostatics
rcoulomb                = 1.0            ; short-range electrostatic cutoff (nm)
pme-order               = 4              ; cubic interpolation
fourierspacing          = 0.125          ; grid spacing for FFT

; Van der Waals
vdwtype                 = cutoff
vdw-modifier            = Potential-shift
rvdw                    = 1.0            ; short-range van der Waals cutoff (nm)
DispCorr                = EnerPres

; Temperature coupling
tcoupl                  = V-rescale      ; thermostat
tc-grps                 = System
tau-t                   = 1.0            ; time constant (ps)
ref-t                   = 300            ; reference temperature (K)

; Pressure coupling
pcoupl                  = C-rescale      ; pressure coupling
pcoupltype              = isotropic      ; scaling of box vectors
tau-p                   = 2.0            ; time constant (ps)
ref-p                   = 1.0            ; reference pressure (bar)
compressibility         = 4.5e-5         ; isothermal compressibility of water (bar^-1)

; Velocity generation
gen-vel                 = no             ; velocity generation off after NVT
EOF

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt_1.tpr -maxwarn 100
gmx mdrun -deffnm npt_1 -bonded gpu -pme gpu

echo "3" | gmx genrestr -f system.gro -n index.ndx -o posre.itp -fc 100 100 100
gmx grompp -f npt.mdp -c npt_1.gro -t npt_1.cpt -r npt_1.gro -p topol.top -o npt_2.tpr -maxwarn 100
gmx mdrun -deffnm npt_2 -bonded gpu -pme gpu

echo "3" | gmx genrestr -f system.gro -n index.ndx -o posre.itp -fc 10 10 10
gmx grompp -f npt.mdp -c npt_2.gro -t npt_2.cpt -r npt_2.gro -p topol.top -o npt_3.tpr -maxwarn 100
gmx mdrun -deffnm npt_3 -bonded gpu -pme gpu

echo "3" | gmx genrestr -f system.gro -n index.ndx -o posre.itp -fc 0 0 0
for ligand in $ligands; do
  gmx genrestr -f "${ligand}.gro" -n "index_${ligand}.ndx" -o "posre_${ligand}.itp" -fc 100 100 100
done
gmx grompp -f npt.mdp -c npt_3.gro -t npt_3.cpt -r npt_3.gro -p topol.top -o npt_4.tpr -maxwarn 100
gmx mdrun -deffnm npt_4 -bonded gpu -pme gpu

for ligand in $ligands; do
  gmx genrestr -f "${ligand}.gro" -n "index_${ligand}.ndx" -o "posre_${ligand}.itp" -fc 10 10 10
done
gmx grompp -f npt.mdp -c npt_4.gro -t npt_4.cpt -r npt_4.gro -p topol.top -o npt_5.tpr -maxwarn 100
gmx mdrun -deffnm npt_5 -bonded gpu -pme gpu

for ligand in $ligands; do
  gmx genrestr -f "${ligand}.gro" -n "index_${ligand}.ndx" -o "posre_${ligand}.itp" -fc 0 0 0
done
gmx grompp -f npt.mdp -c npt_5.gro -t npt_5.cpt -r npt_5.gro -p topol.top -o npt_6.tpr -maxwarn 100
gmx mdrun -deffnm npt_6 -bonded gpu -pme gpu

# Wrap unit cell for visual inspection
echo 0 | gmx trjconv -s em.tpr -f npt_6.gro -o wrap_tmp.gro -pbc whole
echo 1 0 | gmx trjconv -s em.tpr -f wrap_tmp.gro -o wrap_tmp2.gro -pbc nojump -center
echo 0 | gmx trjconv -s em.tpr -f wrap_tmp2.gro -o wrapped.gro -pbc mol -ur compact

rename -f 's/npt_6/npt/' *
rm wrap_tmp* ions* em* nvt* npt_*
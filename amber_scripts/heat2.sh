#!/bin/sh
#SBATCH --error=slurm-%j.err
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cat > heat2 <<EOF
   # Equilibration at 298 (5 ns) NPT
   &cntrl
   imin=0, irest=1, ntx=5,
   ntxo=1, ioutfm=1,
   ntf=2, ntc=2,
   ig=-1, cut=10, tol=0.0005,
   ntb=2, ntp=1, pres0 = 1.0, taup=2.0,
   ntt=3, gamma_ln=1.0,
   nscm=100, iwrap=1,
   tempi=298, temp0=298, vlimit=10.0,
   nstlim=2500000, dt=0.002,
   ntpr=250, ntwr=2500, ntwx=2500,
   &end
EOF

srun pmemd.cuda -O -i heat2 -o heat2.out -r heat2.rst7 -x heat2.mdcrd -c heat1.rst7 -ref heat1.rst7 -p orig.parm7

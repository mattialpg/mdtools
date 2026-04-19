#!/bin/sh
#SBATCH --error=slurm-%j.err
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cat > heat1 <<EOF
   # Heating 100 -> 298 (20 ns) NVT
   &cntrl
   imin=0, irest=0, ntx=1,
   ntxo=1, ioutfm=1,
   ntf=2, ntc=2,
   ig=-1, cut=10, tol=0.0005,
   ntb=1, ntp=0,
   ntt=3, gamma_ln=1.0,
   nscm=100, iwrap=1,
   tempi=100, temp0=298, vlimit=10.0,
   nstlim=10000000, dt=0.002,
   ntpr=1000, ntwr=10000, ntwx=10000,
   nmropt=1, ntr=1,
   restraint_wt=2.0, restraintmask='(@CA,C,O,N)&!:WAT',
   &end
   &wt TYPE='TEMP0', ISTEP1=1, ISTEP2=5000000, VALUE1=100.0, VALUE2=298.0, /
   &wt TYPE='END' /
EOF

srun pmemd.cuda -O -i heat1 -o heat1.out -r heat1.rst7 -x heat1.mdcrd -c heat0.rst7 -ref heat0.rst7 -p orig.parm7

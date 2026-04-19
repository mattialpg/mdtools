#!/bin/sh
#SBATCH --error=slurm-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cat > prod_in <<EOF
    &cntrl
    imin=0, irest=1, ntx=5,
    nstlim=5000000, dt=0.002,
    ntc=2, ntf=2,
    cut=8, ntb=2, ntp=1, taup=2.0,
    ntpr=500, ntwx=5000, ntwr=5000,
    ntt=3, gamma_ln=2.0, iwrap=1,
    temp0=298, vlimit=10.0,
    &end
EOF

srun pmemd.cuda -O -i prod_in -o prod_n.out -r prod_n.rst7 -x prod_n.mdcrd -c prod_m.rst7 -ref prod_m.rst7 -p orig.parm7

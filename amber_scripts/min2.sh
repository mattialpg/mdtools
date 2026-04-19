#!/bin/bash
#SBATCH --error=slurm-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=160
#SBATCH --time=01:00:00

cat > min2 <<EOF  
   # Water minimisation
   &cntrl
   imin = 1,      ntb    = 1,       cut   = 10,    
   ncyc = 5000,   maxcyc = 10000,   iwrap = 1,
   ntpr = 500,    ntwr   = 1000,
   ntr  = 1,
   restraint_wt = 50, restraintmask = '!:WAT&!@Na+',
   &end
EOF

srun sander.MPI -O -i min2 -o min2.out -r min2.rst7 -c min1.rst7 -ref min1.rst7 -p orig.parm7

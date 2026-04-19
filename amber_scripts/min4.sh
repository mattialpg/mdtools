#!/bin/bash
#SBATCH --error=slurm-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=160
#SBATCH --time=01:00:00

cat > min4 <<EOF  
   # Water minimisation
   &cntrl
   imin = 1,      ntb    = 1,       cut   = 10,    
   ncyc = 5000,   maxcyc = 10000,   iwrap = 1,
   ntpr = 500,    ntwr   = 1000,
   ntr  = 1,
   restraint_wt = 5, restraintmask = '!:WAT&!@Na+',
   &end
EOF

srun sander.MPI -O -i min4 -o min4.out -r min4.rst7 -c min3.rst7 -ref min3.rst7 -p orig.parm7

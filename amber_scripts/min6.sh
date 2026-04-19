#!/bin/bash
#SBATCH --error=slurm-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=160
#SBATCH --time=01:00:00

cat > min6 <<EOF  
   # Full minimisation
   &cntrl
   imin = 1,      ntb    = 1,       cut   = 10,    
   ncyc = 5000,   maxcyc = 10000,   iwrap = 1,
   ntpr = 500,    ntwr   = 1000,
   &end
EOF

srun sander.MPI -O -i min6 -o min6.out -r min6.rst7 -c min5.rst7 -ref min5.rst7 -p orig.parm7

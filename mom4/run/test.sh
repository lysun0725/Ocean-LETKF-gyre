#!/bin/sh --login
#SBATCH -n 80         
#SBATCH -t 02:00:00  
#SBATCH -A aosc-hi 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

mpirun -n 20 letkf.DRIFTERS.040

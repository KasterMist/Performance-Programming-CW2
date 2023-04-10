#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1:00:0
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH --account=m22oc-s2312252

module load intel-compilers-19

# Switch to current working directory
cd C

# Run the parallel program
./MD




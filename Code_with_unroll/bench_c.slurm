#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1:00:0
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH --output=slurm-delete_Ndim_5.out
#SBATCH --error=slurm-delete_Ndim.err


#SBATCH --account=m22oc-s2312252

module load intel-compilers-19
module load perf

# Switch to current working directory
cd C

# Run the parallel program
perf stat -e cycles,instructions,cache-misses,branches,branch-misses,page-faults ./MD
# ./MD




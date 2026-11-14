#!/bin/sh

#SBATCH --job-name=CBE735
#SBATCH -q standard         #standard, debug, long, high_priority
#SBATCH -t 24:00:00         #job length
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -n 16
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err

date > a

source /home/software/intel_hpc_toolkit/oneapi/setvars.sh
mpirun /home/software/vasp/euclid-2/bin/vasp-6.4.2-vtst/vasp_std > out

date > b

pwd >> ~/finished_jobs

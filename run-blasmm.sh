#!/bin/bash  
#SBATCH --job-name="blasmm"
#SBATCH --output="job.blasmm.%j.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --export=ALL
#SBATCH -t 00:03:00
#SBATCH --account=csb175

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
./blasmm

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
./blasmm


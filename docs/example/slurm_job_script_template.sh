#!/bin/bash
#SBATCH --account=ACCOUNT
#SBATCH --ntasks=NUMPROCS               # number of MPI processes
#SBATCH --mem-per-cpu=MEMPERCPU      # memory; default unit is megabytes
#SBATCH --time=RUNTIME           # time (DD-HH:MM)
module purge
module load nixpkgs/16.09 gcc/5.4.0 openmpi/1.8.8 fftw-mpi/3.3.6 petsc/3.1-p8 openblas/0.2.20 blacs/1.1 scalapack/1.8.0 python/2.7.14 scipy-stack/2017b
export PATH=$HOME/apps/spherls/1.0/bin:$PATH
mpirun SPHERLS
# Elmag_204
Monte Carlo simulation of three-dimensional electromagnetic in the extragalactic void. Sky image of very high energy point sources on the gamma ray sky.

Written mainly in Fortran 90.
Module for MPI compile: 
module load mpi/openmpi-x86_64

MPI compile command: 
make mpi

Execute for X desired processors: 
mpirun -np X a.out

Single processor compile command: 
make single

Execute: 
./a.out

Running time varies with the desired parameters used in the user204.90 file.

For an extensive description of the program, please see my Master thesis (link will come).

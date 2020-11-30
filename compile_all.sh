#!/bin/bash

echo '------------------------------'
gcc -O2 Project-Serial-Original.c -lm -o code-SO.o
echo 'Compiled Serial-Original'

gcc -O2 Project-Serial-Updated.c -lm -o code-SU.o
echo 'Compiled Serial-Updated'

gcc -O2 -fopenmp Project-Openmp-Coupled.c -lm -o code-OMP-C.o
echo 'Compiled OpenMP-Coupled'

gcc -O2 -fopenmp Project-Openmp-Domain-Decomposition.c -lm -o code-OMP-D.o
echo 'Compiled OpenMP-Domain-Decomposition'

gcc -O2 -fopenmp Project-Openmp-Both.c -lm -o code-OMP-B.o
echo 'Compiled OpenMP-Coupled'

mpicc -O2 Project-MPI.c -lm -o code-MPI.o
echo 'Compiled MPI'



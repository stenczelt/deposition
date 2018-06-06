#/bin/bash -l

mkdir -p bin
cd src/

gfortran -c neighbors.f90
gfortran -fopenmp -o ../bin/prep_for_impact prep_for_impact.f90 *.o

cd ..

# Parallel
Parallel Algorithm Implementations
To compile the program, run the following command:
mpic++ -std=c++11 -o P Assignment4-MPIMatrixMultiplication.cpp

Please note that you must include the -std=c++11 flag as the random number engine requires the program be compiled with at least C++ version 11.


Then run the script on the hpc shell as 
mpiexec -np <NumberOfProcessors> ./P <MatrixSize>

Example:
mpiexec -np 16 P 2000

to run with 16 threads and 2000 as the row and column size of each matrix (2000x2000).

The script will output the row numbers split between the threads as well as the run time.

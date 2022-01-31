# Overview
AGM is a program for computing the arithmetic-geometric mean of complex numbers. It also computes a new function called khe (խ(z)) such that the AGM of խ(z) and խ(z+πi) is 1 for all z in the left half-plane, and one step of the AGM of these numbers produces խ(2z) as the arithmetic mean and խ(2z+πi) as the geometric mean.

# Compiling
To compile, if you're not developing the program:

1. Create a subdirectory build/ inside the directory where you untarred the source code.
2. `cd build`
3. `cmake ..`
4. `make`

If you are developing the program:

1. Create a directory build/agm outside the directory where you cloned the source code.
2. `cd build/agm`
3. `cmake <directory where the source code is>`
4. `make`

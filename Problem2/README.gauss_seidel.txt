Instructions for executing the Gauss-Seidel solver simulation for Problem 2.

Components:
    gauss_seidel.java (source file)
    gauss_seidel.class (compiled class file)

Other files referred from within the folder:
    Jama library (within Problem 2 folder, DO NOT MOVE)
        Jama.Matrix

Java external libraries imported (from outside folder):
    java.io.File
    java.util.Scanner
    java.util.ArrayList
    java.lang.Math

PURPOSE:
    solves a linear system in the matrix form Ax=b using Gauss Seidel iterative method.

TO RUN:
    Change current directory in the command module to that of the Problem 2 folder.
    type/run "java gauss_seidel" in the command line

    Type in the name of .dat or other simple text (such as .txt) file (complete address) that holds the desired augmented (nx(n+1)) [A|x] matrix. This should be written with each line containing a row, with separate entries separated by spaces. An example has been provided in the Problem2 folder, named "test_aug_matrix.dat".

    Type the desired error bound for solving the system, in decimal form (such as ".00000001" would be 1e-8).

    The final estimated x vector will be printed, along with the number of iterations used to solve the method within the error bound.
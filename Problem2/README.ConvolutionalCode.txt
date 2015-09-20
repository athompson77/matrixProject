Instructions for executing the Convolutional Code simulation for Problem 2.

Components:
    CompoutationalCode.java (source file)
    CompoutationalCode.class (compiled class file)

Other files imported from within the folder:
    Jama library (within Problem 2 folder, DO NOT MOVE)
        Jama.Matrix

Java external libraries imported (from outside folder):
    java.io.File
    java.util.Scanner
    java.util.ArrayList
    java.lang.Random

PURPOSE:
    encodes and decodes binary streams created by convolutional codes

TO RUN:
    Change current directory in the command module to that of the Problem 2 folder.
    type/run "java CompoutationalCode" in the command module

    type/run "e" to have the code encode a random binary stream
        type a number (ex: 5), which is to be the length of x (the input stream)
        the randomly created input will be printed, underneath will be the encoded y stream (output, composed of the two streams y0 and y1)

        type "n" if you wish to end the program

        -OR-

        type "y" if you wish to re-decode the created y stream using examples of Jacobi and Gauss-Seidel iterative methods
            enter a decimal representing the desired error bound (ex: ".00000001")

            results from decoding using A0 and y0, and A1 and y1, with both Jacobi and Gauss-Seidel will be printed, along with the number of iterations taken (or 100, if a norm to compute error could not be computed)

    -OR-

    type/run "d" to skip to decode a binary matrix of own design
    (note: this will only work if the stream is possible to have been the result of encoding in the first place)

    type in the desired output stream, following the format listed on screen (no extraneous spaces except those inbetween entry pairs, and no commas)

    type in the desired guess stream, following the format listed on screen (no extraneous spaces except those inbetween entries, and no commas, must be same row length as the previously inputted y vector)

    Type the desired error bound for solving the system, in decimal form (such as ".00000001" would be 1e-8).

    results from decoding (as vectors representing x) using A0 and y0, and A1 and y1, with both Jacobi and Gauss-Seidel will be printed, along with the number of iterations taken (or 100, if a norm to compute error could not be computed)

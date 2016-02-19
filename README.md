Documentation:
There are two functions, strassen does strassen matrix multiplication and “Multiply” does iterative multiplication.
Dimensions are random in the power of 2 and the matrix is populated randomly.
To run simply run the program and the main function will be called and both the strassen and iterative multiplication will run.

Both matrix are printed.
Code was tested at a matrix size of 64 x 64.
Since this is a small matrix, strassen takes 6 seconds and iterative takes 0.07 seconds.

 

Public repo link : https://github.com/sohaib3k/Matrix-multiplication_Strassen/tree/master


Strassen code taken from :
https://github.com/MartinThoma/matrix-multiplication/blob/master/C%2B%2B/strassen-algorithm.cpp

Testing has been done using a class named testMatrix which has a function checkMatrixMutiplication.
This calls both Strassen and Iterative functions and subtracts their 2d vectors. If the resulting 2d vector is all zeroes, the test case passes. Otherwise it fails.


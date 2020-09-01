# gaussian-elimination
Program for solving system of equations

Gaussian elimination is a method for solving matrix equations of the form ax=b. 
First step is to convert the system of equations into an augmented matrix (the matrix contains all of the information in the system of equations without the x, y, z etc. labels to carry around).
Next step is to perform elementary row operations to put the augmented matrix into the upper triangular form, eliminating variables until only one variable is left.
Once this final variable is determined, its value is substituted back into the other equations in order to evaluate the remaining unknowns.

This program implements Gaussian elimination by reading coefficients from text files, performing the algorithm described above and in the end giving either the solution or a reason why solution can not be found.

# Numerical_Algebra
Build a matrix class here to implement routine operation
In this program, I create a Matrix class to store vital information and operations.
As for solving linear equations, I pack all dependent functions into below function,

    sle(Matrix b, int method)

Here the method is a selective parameter can be defined as

    nonpivot: corresponds 0(selective)
    colpivot: corresponds 1(default)
    allpivot: corresponds 2(not ready yet)

So the syntax is like

    Matrix A(row1, col1, array1[]);
    Matrix B(row2, 1, array2[]);
    Matrix X = A.sle(B, nonpivot);


Analysis:
To solve linear equations efficiently, like AX = B, where A is not singular and B is
a column matrix, we must decompose A into a lower triangle matrix L and an upper
matrix U. There are three ways to implement this operation, so users can choose an
appropriate one when calling LU function. The LU funtion returns three Matrices in
one vector, the first one stores the result, the second and third one are permutation
matrices which represent row and column interchange respectively. This problem turns
out to solve LUX = PBQ. Then we denote y as the solution to Ly = PBQ and UX = y.
They are much simplier form, named triangle matrix, to solve.

Dependence tree as below:

Class Matrix:

    ----properties:
        ----unsigned int row
        ----unsigned int col
        ----vector< vector<double> > mat

    ----function:
        ----Matrix sle(Matrix X)
            ----Matrix lowtri_sle(Matrix y)
            ----Matrix uptri_sle(Matrix x)
            ----vector<Matrix> LU(int method)
            ----int max_index(const vector<double>& vec)
        ----print()
        ----trans()
        ----Matrix operator*(const Matrix& A)

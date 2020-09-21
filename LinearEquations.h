#ifndef _JOVI_WONG_
#define _JOVI_WONG_
#define forward_method 1
#define backward_method 0
#include <iostream>
#include <vector>

using namespace std;
/*Index of this matrix is equavalent with math expression*/
class Matrix
{
public:
    unsigned int row;
    unsigned int col;
    vector< vector<double> > mat;
    int trans();
    int gauss();
    int print();

    Matrix(double init[], unsigned int m, unsigned int n);
    ~Matrix(){};
};

/*principle of nominating function 
 * first argument is the type of number
    -d:double
    -s:single
    -l:long
    -i:int
 *second argument is the type of matrix
    -lowtri: lower triangular matrix
    -uptri: upper triangular matrix 
 * third argument is the abbreviation of the operation
    -sle: solve linear equations
*/

Matrix d_lowtri_sle(unsigned int leading_dimension, double L[], double b[]); 
/* It must receive an argument describing the scale of katrix 
 * because we can't know how long is the input number array*/

Matrix d_lowtri_sle(Matrix L, Matrix b);

Matrix d_uptri_sle(Matrix L, Matrix b);

#endif

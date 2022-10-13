#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <vector>
#include <exception>
#include <cmath>

#define zeros 0
#define identity 1

#define LL 1
#define LDL 2

enum class LUPivotType
{
	NONPIVOT = (int)0,
	COLPIVOT = (int)1,
	ALLPIVOT = (int)2
};

enum class MatrixInitType
{
	ZEROS = (int) 0,
	IDENTITY = (int) 1
};

using namespace std;


/*Index of this matrix is equavalent with math expression*/
class Matrix
{
public:
    unsigned int row;
    unsigned int col;
    vector< vector<double> > mat;
    
    Matrix(unsigned int m, unsigned int n, double init[]);
    Matrix(unsigned int m, unsigned int n, MatrixInitType type = MatrixInitType::ZEROS);
    ~Matrix(){};

    Matrix trans();
    int print();
    virtual Matrix sle(Matrix b, LUPivotType LU_method = LUPivotType::COLPIVOT);
    Matrix operator*(const Matrix& A);
    Matrix operator-(const Matrix& B);
    friend Matrix operator*(double cof, Matrix& B);

    Matrix lowtri_sle(Matrix b);//function in sle(Matrix b)
    Matrix uptri_sle(Matrix y);//function in sle(Matrix b)
    Matrix diag_sle(Matrix b); // function for sle
    virtual vector<Matrix> LU(LUPivotType method = LUPivotType::COLPIVOT);//function in sle(Matrix b)
    int max_index(const vector<double>& vec);//return the index of maximum

};



class SymmetricalPositiveMatrix: public Matrix
{
public: 
   SymmetricalPositiveMatrix(unsigned int m, unsigned int n, double init[]);
   virtual Matrix sle(Matrix b, int method=LDL);
   virtual vector<Matrix> LU(int method=LDL);
};







/*principle of nominating function 
 * first argument is the type of number
    -d:double
    -s:singlo
    -l:long
    -i:int
 *second argument is the type of matrix
    -lowtri: lower triangular matrix
    -uptri: upper triangular matrix 
 * third argument is the abbreviation of the operation
    -sle: solve linear equations
*/

//Matrix lowtri_sle(unsigned int leading_dimension, double L[], double b[]); 
/* It must receive an argument describing the scale of katrix 
 * because we can't know how long is the input number array*/

#endif

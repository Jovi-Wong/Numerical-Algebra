#include "LinearEquations.h"

int main()
{
    double B[9] = {1,4,7,
                   2,5,8,
                   3,6,10};

//    double C[3] = {3,2,1};
    Matrix A(B,3,3);         
    A.gauss();
    return 0;
}

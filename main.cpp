#include "Matrix.hpp"
#include "time.h"

int main()
{
/*HOMEWORK_1_
    double a[84*84]={0};
    double b[84];

    for(int i=0; i<84; ++i)
    {
        switch(i)
        {
            case 0:
            {
                b[i] = 7;
                a[0] = 6;
                a[1] = 1;
                break;
            }

            case 83:
            {
                b[i] = 14;
                a[84*i+82] = 8;
                a[84*i+83] = 6;
                break;
            }

            default:
            {
                b[i] = 15;
                int pos = 85*i;
                a[pos-1] = 8;
                a[pos] = 6;
                a[pos+1] = 1;
                break;
            }
        }   
    }

//  double c[3] = {3,2,1};
    Matrix A(84,84,a);         
    Matrix B(84,1,b);
//    A.print();
    cout << "Here are answers by selecting no pivot: " << endl;
    Matrix X = A.sle(B, nonpivot);
    X.print();
    cout << "Here are answers by selecting column pivot: " << endl;
    X = A.sle(B, colpivot);
    X.print();
//  A.gauss();
*/

/*example
    double a[16] = 
    {
     16, 4, 8, 4,
      4,10, 8, 4,
      8, 8,12,10,  
      4, 4,10,12
    };
    
    double b[4] = 
    {
        32,
        26,
        38,
        30
    };

    SymmetricalPositiveMatrix A(4,4,a);
    Matrix B(4,1,b);
    Matrix x = A.sle(B,2);
    x.print();
*/

/*
    double a[40*40] = {0};
    for(int i=0; i<40; ++i)
    {
        switch(i)
        {
            case 0:
            {
                a[0] = 10;
                a[1] = 1;
                break;
            }

            case 39:
            {
                a[39*38] = 1;
                a[39*40+39] = 10;
            }

            default:
            {
                a[i*40+i-1] = 1;
                a[i*40+i] = 10;
                a[i*40+i+1] = 1;
            }
        }
    }

    double b[40] = {0};
    for(int i=0; i<40; ++i)
    {
        for(int j=1; j<=40; ++j)
        {
            b[i] = b[i] + (1.0/(i+j));
        }
    }


    clock_t start,end;
    cout << "solve linear equations by column pivot method:" << endl;
    Matrix A(40,40,a);
    Matrix B(40,1,b);
    start = clock();
    Matrix x = A.sle(B,1);
    end = clock();
    x.print();
    double dur = end - start;
    cout << "using time " << dur << endl << endl;
    cout << "solve linear equations by " << endl;
    SymmetricalPositiveMatrix SA(40,40,a);
    start = clock();
    Matrix Sx = SA.sle(B,2);
    end = clock();
    dur = end - start;
    Sx.print(); 
    cout << "using time " << dur << endl << endl;
    Matrix diff = x - Sx;
    cout << "First answer substracts the second answer: " << endl;
    diff.print();
*/

    double a[4] = {375,374,752,750};
    double b[2] = {1,1};
    Matrix A(2,2,a);
    Matrix B(2,1,b);
    Matrix x = A.sle(B);
    x.print();
    return 0;
}

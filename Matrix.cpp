#include "Matrix.h"

#define nonpivot 0
#define colpivot 1
#define allpivot 2

#define zeros 0
#define identity 1

#define LL 1
#define LDL 2


//construct a m x n Matrix instance with initial numbers stored in init[]
Matrix::Matrix(unsigned int m, unsigned int n, double init[]) 
{
    row = m;
    col = n;
    vector<double> temp(1,0);
    vector<double> empty(col+1,0);
    mat.push_back(empty);
    for(int i=1; i<=row; i++)
    {
        mat.push_back(temp);
        for(int j=1; j<=col; j++)
        {
            mat[i].push_back(init[(i-1)*col+j-1]);
        }
        mat[i].shrink_to_fit();
    }
//  print();
}



//construct a m x m identity matrix
Matrix::Matrix(unsigned int m, unsigned int n, int type)
{
    switch(type)
    {
        case 0:
        {
            row = m;
            col = n;
            vector<double> empty(n+1,0);
            empty.shrink_to_fit();
            for(int i=0; i<=m; ++i)
            {
                mat.push_back(empty);
            }
            break;
        }
        
        case 1:
        {
            if(m != n)
            {
                cerr << "Identity matrix' row must equal to its column!";
            }
            row = m;
            col = m;
            vector<double> empty(m+1,0);
            empty.shrink_to_fit();
            mat.push_back(empty);
            for(int i=1; i<=m; ++i)
            {
                mat.push_back(empty);
                mat[i][i] = 1;
            }
            break;
        }

        default:
        {
            cerr << "Please input legal matrix type!"<< endl;
            break;
        }
    }
}


//output the transposed matrix without interchanging values itsef
Matrix Matrix::trans()
{
    Matrix newmat(row,col,zeros);
    for(int i=1; i<=col; i++)
    {
        for(int j=1; j<=row; j++)
        {
            newmat.mat[i][j] = mat[j][i];
        }
    }
//  new.print();
    return newmat;
}


//print values in the matrix instance on the screen
int Matrix::print()
{
    for(vector< vector<double> >::iterator pos1 = mat.begin()+1; pos1 != mat.end(); pos1++)
    {
        for(vector<double>::iterator pos2 = (*pos1).begin()+1; pos2 != (*pos1).end(); pos2++)
        {
            cout << (*pos2) << " ";
        }
        cout << endl;
    }
    cout << endl;
    return 0;
}



//solve linear equations represented by low triangle matrix 
Matrix Matrix::lowtri_sle(Matrix b)
{
    if(row != col)
    {
        cerr << "input matrix is illegal"<< endl;
    }

    unsigned int n = row;
    for(int j=1; j<=n-1; ++j) 
    {
        b.mat[j][1] = b.mat[j][1]/mat[j][j];
        for(int k=j+1; k<=n; ++k)
        {
            b.mat[k][1]=b.mat[k][1]-b.mat[j][1]*(this->mat[k][j]);
        }
    }
    b.mat[n][1] = b.mat[n][1]/(this->mat[n][n]);
//  b.print();
    return b;
}

//solve linear equations represented by up triangle matrix

Matrix Matrix::uptri_sle(Matrix y)
{
    unsigned int n = row;
    for(int j=n; j>=2; j=j-1)
    {
        y.mat[j][1] = y.mat[j][1]/mat[j][j];
        for(int k=1; k<j; ++k)
        {
            y.mat[k][1] = y.mat[k][1]-y.mat[j][1]*mat[k][j];
        }
    }
    y.mat[1][1] = y.mat[1][1]/mat[1][1];
    return y;
}



//output the Low triangle and up triangle matrix stored in one matrix
vector<Matrix> Matrix::LU(int method)
{
    if(row != col)
    {
        cerr << "Gauss Method only applies to square matrix" << endl;
    }
    
    vector<double> u(row+1,0);
    vector<double> t(col+1,0);
    vector<Matrix> LUPQ;

    for(int i=1; i<row; ++i)
    {
        u[i] = i;
        t[i] = i;
    }

    switch(method)
    {
        case 0:
        {
//            cout << "enter nonpivot LU" << endl;
            int n = row;
            Matrix rst = *this;
            for(int k=1; k<n; ++k)
            {
                for(int i=k+1; i<=n; ++i)
                {
                    try
                    {
                        rst.mat[i][k] = rst.mat[i][k]/rst.mat[k][k];
                    }
                    catch(exception& e)
                    {
                        cerr << "This matrix doesn't support non-pivot method!" <<endl;
                    }

                    for(int j=k+1; j<=n; ++j)
                    {
                        rst.mat[i][j] = rst.mat[i][j]-rst.mat[i][k]*rst.mat[k][j];
                    }
                }
            }
//            cout << "print LU matrix by nonpivot method" << endl;
//            rst.print();
            LUPQ.push_back(rst);
            break;
        }

        case 1:
        {
//            cout <<"enter colpivot method" <<endl;
            int n = row;
            Matrix rst = *this;
            //firstly, the max value in the column
            for(int k=1; k<n; ++k)
            {
                int p = 0;
                double temp_max = 0;
                for(int i=k; i<=n; ++i)
                {
                    if(abs(rst.mat[i][k]) > temp_max)
                    {
                        temp_max = rst.mat[i][k];
                        p = i;
                    }
                }
                swap(rst.mat[p], rst.mat[k]);
                u[k] = p; //store permutation matrix 
                
                if(rst.mat[k][k] != 0)
                {
                    for(int i=k+1; i<=n; ++i)
                    {
                        rst.mat[i][k] = rst.mat[i][k]/rst.mat[k][k];
                    }
                    for(int i=k+1; i<=n; ++i)
                    {
                        for(int j=k+1; j<=n; ++j)
                        {
                            rst.mat[i][j] = rst.mat[i][j]-rst.mat[i][k]*rst.mat[k][j];
                        }
                    }
                }
                else
                {
                    cerr <<"Error! This is a singular matrix";
                }
//                cout << k << " step in colpivot method" << endl;
//                rst.print();
            }
//            cout << "print LU matrix by colpivot method" << endl;
//            rst.print();
            LUPQ.push_back(rst);
            break;
        }

        case 2:
        {
            cerr << "All-pivot method hasn't been ready yet!" << endl;
            break;
        }
        
        default:
        {
            cerr << "Please input valid method name!" << endl;
            break;
        }
    }
    
    Matrix P(row,row,identity);
    Matrix Q(row,row,identity);

    for(int i=1; i<row; ++i)
    {
        swap(P.mat[i],P.mat[u[i]]);
        for(int j=1; j<col; j=j-1)
        {
             swap(Q.mat[j][i],Q.mat[j][t[i]]);
        } 
    }
//    cout << "print matrix P:" << endl;
//    P.print();
//    cout << "print matrix Q:" << endl;
//    Q.print();
    LUPQ.push_back(P);
    LUPQ.push_back(Q);
    return LUPQ;
}


int Matrix::max_index(const vector<double>& vec)
{
    int idx = 1;
    double m = 0;

    for(int i=1; i<vec.size();++i)
    {
        if (abs(vec[i]) > m)
        {
            idx = i;
            m = vec[i];
        }
    }
    return idx;
}


//solve linear equations 
//syntax such as x = A.sle(b) to solve Ax=b
Matrix Matrix::sle(Matrix b, int LU_method)
{
    vector<double> u(row+1,0);
    vector<Matrix> plu = LU(LU_method);
    Matrix L(row, col, identity);
    Matrix U(row, col, zeros);
//    cout << "print intial U" << endl;
//    U.print();
//    Matrix P = plu[1];
    //extract low triangle matrix
    for(int i=1; i<=row;++i)
    {
        for(int j=1; j<=col; ++j)
        {
            if(j<i)
            {
                L.mat[i][j]=plu[0].mat[i][j];
            }
            else
            {
                U.mat[i][j]=plu[0].mat[i][j];
            }
        }
    }
//    cout << "final L" << endl;
//    L.print();
//    cout << "final U" << endl;
//    U.print();
//    cout << "Show matrix Pb" << endl;
//    Pb.print();
    Matrix y = L.lowtri_sle(plu[1]*b);
//    cout << "print y" << endl;
//    y.print();
    Matrix x = U.uptri_sle(y);
    return x;
}

Matrix Matrix::operator*(const Matrix& A)
{
//    cout << "the input A is " << endl;
//    A.print();
//    cout << "this -> row = " << this->row << endl;
//    cout << "A.col = " << A.col << endl;
    Matrix rst(this->row, A.col, zeros);
//    cout << "intial result matrix" << endl;
//    rst.print();
//    (*this).print();
    for(int i=1; i<=rst.row; ++i)
    {
        for(int j=1; j<=rst.col; ++j)
        {
            for(int c=1; c<=(this->col); ++c)
            {
                rst.mat[i][j] = rst.mat[i][j] + (this->mat[i][c])*(A.mat[c][j]);
            }
        }
    }
    return rst;
}

Matrix Matrix::diag_sle(Matrix b)
{
    for(int i=1; i<=row; ++i)
    {     
        b.mat[i][1] = b.mat[i][1] / mat[i][i];
    }
    return b;
}


//Matrix lowtri_sle(unsigned int leading_dimension, double L[], double b[])
//{
//    for(int j=0; j<leading_dimension-1; ++j) 
//    {
//        b[j] = b[j]/L[j*leading_dimension+j];
//        for(int k=j+1; k<leading_dimension; ++k)
//        {
//            b[k]=b[k]-b[j]*L[k*leading_dimension+j];
//        }
//    }
//    b[leading_dimension-1] = b[leading_dimension-1]/L[leading_dimension*leading_dimension-1];

//    Matrix rst(leading_dimension,1,b);
//    return rst;
//}




Matrix operator*(double cof, Matrix& B)
{
    Matrix rst(B.row, B.col, zeros);
    for(int i=1; i<=B.row; ++i)
    {
        for(int j=1; j<=B.col; ++j)
        {
            rst.mat[i][j] = cof*B.mat[i][j];
        }
    }

    return rst;
}

SymmetricalPositiveMatrix::SymmetricalPositiveMatrix(unsigned int m, unsigned int n, double init[]): Matrix(m, n, init){}

Matrix SymmetricalPositiveMatrix::sle(Matrix b, int method)
{
   Matrix x = b;
   switch(method)
   {
       case 1:
       {
//           cout << "implement method 1" << endl;
           vector<Matrix> LLT = LU(method);
           Matrix y = LLT[0].lowtri_sle(b);
           x = LLT[1].uptri_sle(y);
           break;
       }

       case 2:
       {
//           cout << "implement method 2" << endl;
           vector<Matrix> LDLT = LU(method);
//           cout << "LU disposition complete!" << endl;
           Matrix z = LDLT[0].lowtri_sle(b);
           Matrix y = LDLT[1].diag_sle(z);
           x = LDLT[2].uptri_sle(y);
           break;
       }
    }
   return x;
}


vector<Matrix> SymmetricalPositiveMatrix::LU(int method)
{
    int n = row;
    Matrix A = (*this);
    vector<Matrix> rst;

    switch(method)
    {
        case 1:
        {
            for(int k=1; k<=n; ++k)
            {
                A.mat[k][k] = sqrt(A.mat[k][k]);
                for(int i=k+1; i<=n; ++i)
                {
                    A.mat[i][k] = (A.mat[i][k]) / (A.mat[k][k]);
                }

                for(int j=k+1; j<=n; ++j)
                {
                    for(int t=j; t<=n; ++t)
                    {
                        A.mat[t][j] = A.mat[t][j] - A.mat[t][k]*A.mat[j][k];
                    }
                }
            }
            
            Matrix L(n,n,identity);

            for(int i=1; i<=n; ++i)
            {
                for(int j=1; j<=i; ++j)
                {
                    L.mat[i][j] = A.mat[i][j];
                }
            }
//            cout << "Matrix L as below" << endl;
//            L.print();
            Matrix LT(n,n,identity);
//            cout << "Matrix LT as below" << endl;
            LT = L.trans();
//            LT.print();
            rst.push_back(L);
            rst.push_back(LT);
            break;
        }
        

        case 2:
        {
//            cout << "LU case 2" << endl;
            for(int j=1; j<=n; ++j)
            {
//                cout << "j = " << j << endl;
                Matrix v(1,j-1,zeros);
                if(j != 1)
                {    
                    for(int i=1; i<=j-1; ++i)
                    {
                        v.mat[1][i] = A.mat[j][i]*A.mat[i][i];
                    }

                    for(int k=1; k<=j-1; ++k)
                    {
                        A.mat[j][j] = A.mat[j][j]-A.mat[j][k]*v.mat[1][k];
                    }
//                    cout << "j = " << j << " block 1 finished!" << endl;
                
                    for(int k=j+1; k<=n; ++k)
                    {
                        double temp = 0;
                        for(int t=1; t<=j-1; ++t)
                        {
                            temp = temp + A.mat[k][t]*v.mat[1][t];
                        }
                        A.mat[k][j] = (A.mat[k][j]-temp)/A.mat[j][j];
                    }
//                    cout << "j = " << j << " block 2 complete!" << endl;
                }
                else
                {
                    for(int k=2; k<=n; ++k)
                    {
                        A.mat[k][1] = A.mat[k][1]/A.mat[1][1];
                    }
                }
            }

            Matrix L(n,n,identity);
            Matrix D(n,n, zeros);

            for(int i=1; i<=n; ++i)
            {
                for(int j=1; j<i; ++j)
                {
                    L.mat[i][j] = A.mat[i][j];
                }
                D.mat[i][i] = A.mat[i][i];
            }
//            cout << "Matrix L as below by LDL" << endl;
            rst.push_back(L);
//            L.print();
//            cout << "Matrix D as below by LDL " << endl;
            rst.push_back(D);
//            D.print();
            rst.push_back(L.trans());
            break;
        }
    }
    return rst;
}


Matrix Matrix::operator-(const Matrix& B)
{
    Matrix rst(row,col,zeros);
    for(int i=1; i<=row; ++i)
    {
        for(int j=1; j<=col; ++j)
        {
            rst.mat[i][j] = mat[i][j]-B.mat[i][j];
        }
    }
    return rst;
}

#include "LinearEquations.h"

Matrix::Matrix(double init[], unsigned int m, unsigned int n)
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
    }
    print();
}

int Matrix::trans()
{
    vector< vector<double> > newmat;
    vector<double> temp(1,0);
    vector<double> empty(col+1,0);
    newmat.push_back(empty);
    
    for(int i=1; i<=col; i++)
    {
        newmat.push_back(temp);
        for(int j=1; j<=row; j++)
        {
            newmat[i].push_back(mat[j][i]);
        }
    }
    swap(row,col);
    mat = newmat;
    print();
    return 0;
}

int Matrix::gauss()
{
    if(row != col)
    {
        cerr << "Gauss Method only applies to square matrix" << endl;
    }

    int n = row;

    for(int k=1; k<n; ++k)
    {
        for(int i=k+1; i<=n; ++i)
        {
            mat[i][k] = mat[i][k]/mat[k][k];
            for(int j=k+1; j<=n; ++j)
            {
                mat[i][j] = mat[i][j]-mat[i][k]*mat[k][j];
            }
        }
    }
    print();
    return 0;
}

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

Matrix d_lowtri_sle(Matrix L, Matrix b)
{
    if(L.row != L.col)
    {
        cerr << "input matrix is illegal"<< endl;
    }

    unsigned int n = L.row;
    for(int j=1; j<=n-1; ++j) 
    {
        b.mat[j][1] = b.mat[j][1]/L.mat[j][j];
        for(int k=j+1; k<=n; ++k)
        {
            b.mat[k][1]=b.mat[k][1]-b.mat[j][1]*L.mat[k][j];
        }
    }
    b.mat[n][1] = b.mat[n][1]/L.mat[n][n];
    b.print();
    return b;
}

Matrix d_lowtri_sle(unsigned int leading_dimension, double L[], double b[])
{
    for(int j=0; j<leading_dimension-1; ++j) 
    {
        b[j] = b[j]/L[j*leading_dimension+j];
        for(int k=j+1; k<leading_dimension; ++k)
        {
            b[k]=b[k]-b[j]*L[k*leading_dimension+j];
        }
    }
    b[leading_dimension-1] = b[leading_dimension-1]/L[leading_dimension*leading_dimension-1];

    Matrix rst(b,leading_dimension,1);
    return rst;
}

Matrix d_uptri_sle(Matrix U, Matrix y)
{
    unsigned int n = U.row;
    for(int j=n; j>=2; --j)
    {
        y.mat[j][1] = y.mat[j][1]/U.mat[j][j];
        for(int k=1; k<j; ++k)
        {
            y.mat[k][1] = y.mat[k][1]-y.mat[j][1]*U.mat[k][j];
        }
    }
    y.mat[1][1] = y.mat[1][1]/U.mat[1][1];
    y.print(); 
    return y;
}

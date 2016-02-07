#include <iostream>
#include <cmath>
#include "lib/matrixman.h"
#define EPS 0.000001
using namespace MatrixMan;

int main()
{
    int m, n, i, j;
    std::cout << "n=";
    std::cin >> n;
    std::cout << "Matrix:\n";
    Matrix<double> A(n,n,0), b(n,1,0);
    for (i=0;i<n;i++)
        for (j=0;j<n;j++)
            std::cin >> A.get(i,j);
    std::cout << "b:\n";
    for (i=0;i<n;i++)
        std::cin >> b.get(i,0);

    Matrix<double> X;


    for (int k=0;k<n-1;k++)
        if (fabs(A.get(k,k)) > EPS)
        {
            A.slice(k+1,n-1,k,k) /= A.get(k,k);
            A.slice(k+1,n-1,k+1,n-1) -= A.slice(k+1,n-1,k,k) * A.slice(k,k,k+1,n-1);
            b.slice(k+1,n-1,0,0) -= A.slice(k+1,n-1,k,k)*b.get(k,0);
            X = Matrix<double>(n-1-k,1,0);
            A.slice(k+1,n-1,k,k) = X;
        }
        else
        {
            std::cout << "Pivot too small.";
            return 0;
        }

    for (int i=0;i<n;i++)
    {
        if (fabs(A.get(i,i)) < EPS)
        {
            std::cout << "Badly-conditioned matrix.";
            return 0;
        }
    }

    b.get(n-1,0) = b.get(n-1,0) / A.get(n-1,n-1);
    double sum;
    for (int i=n-2;i>=0;i--)
    {
        sum = (A.slice(i,i,i+1,n-1) * b.slice(i+1,n-1,0,0)).get(0,0);
        b.get(i,0) = (b.get(i,0) - sum) / A.get(i,i);
    }

    std::cout << "The solutions are:\n";
    b.display();
}

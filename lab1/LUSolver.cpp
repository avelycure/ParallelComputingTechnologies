#include "header.hpp"

LUSolver::LUSolver(int sizeParam)
{
    size = sizeParam;
    l = new double[size * size];
    u = new double[size * size];
    m = new double[size * size];
}

void LUSolver::resetMatrix(double *&x)
{
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            x[i * size + j] = 0.0;
}

void LUSolver::setDiagonal(double *&x)
{
    for (int i = 0; i < size; i++)
        x[i * size + i] = 1.0;
}

/**
 * Multiplication of two matrices with result in variable res
 */
void LUSolver::multiply(double *&x, double *&y, double *&res)
{
    double sum;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
        {
            sum = 0.0;
            for (int k = 0; k < size; k++)
                sum = sum + x[i * size + k] * y[k * size + j];

            res[i * size + j] = sum;
        }
}

void LUSolver::findDecomposition(double *&a)
{
    //set matrix non diagonal elements of matrix L to zeros, and diagonal to 1
    resetMatrix(l);
    setDiagonal(l);

    //init matrix U as matrix A
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            u[i * size + j] = a[i * size + j];

    for (int i = 0; i < size; i++)
    {
        resetMatrix(m);
        setDiagonal(m);

        for (int k = i + 1; k < size; k = k + 1)
            m[k * size + i] = -u[k * size + i] / u[i * size + i];

        multiply(m, u, u);

        for (int k = i + 1; k < size; k = k + 1)
            m[k * size + i] *= -1;

        multiply(l, m, l);
    }
}

void LUSolver::printMatrix(double *&x, std::string matrixName)
{
    cout << matrixName << endl;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
            cout << x[i * size + j] << " ";
        cout << endl;
    }
    cout << endl;
}

void LUSolver::printDecomposition()
{
    printMatrix(l, "Matrix L");
    printMatrix(u, "Matrix U");
    multiply(l, u, m);
    printMatrix(m, "Matrix A");
}

LUSolver::~LUSolver()
{
    delete[] l;
    delete[] u;
    delete[] m;
}

double* LUSolver::solveWith(double *&b)
{
    double *x = new double[size];
    double *y = new double[size];
    forwardSubstitution(b, y);
    backwordSubstitution(b, x, y);
    return x;
}

void LUSolver::forwardSubstitution(double *&b, double *&y)
{
    double sum;
    for (int i = 0; i < size; i++)
    {
        sum = 0.0;
        for (int j = 0; j < i; j++)
            sum = sum + l[i * size + j] * y[j];
        y[i] = (b[i] - sum) ; // divided on l[i * size + i], but it is 1.0
    }
}

void LUSolver::backwordSubstitution(double *&b, double *&x, double *&y)
{
    double sum;
    for (int i = size - 1; i >= 0; i--)
    {
        sum = 0.0;
        for (int j = size - 1; j > i; j--)
            sum = sum + u[i * size + j] * x[j];
        x[i] = (y[i] - sum) / u[i * size + i];
    }
}

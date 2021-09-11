#include "header.hpp"

LUSolver::LUSolver(int sizeParam)
{
    this->size = sizeParam;
    l = new double[size * size];
    u = new double[size * size];
    x = new double[size * size];
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

void LUSolver::printMatrix(double *&x, string matrixName)
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
}
#include "../header.hpp"

CommonLUDecomposer::CommonLUDecomposer(int sizeParam)
{
    size = sizeParam;
    a = new double[size * size];
}

CommonLUDecomposer::~CommonLUDecomposer()
{
    delete[] a;
}

void CommonLUDecomposer::setMatrix(double *x)
{
    for (int i = 0; i < size * size; i++)
        a[i] = x[i];
}

/**
 * LU decomposition in case of square matrix
 */
void CommonLUDecomposer::findDecomposition()
{
    for (int i = 0; i < size; i++)
        for (int j = i + 1; j < size; j++)
        {
            a[j * size + i] /= a[i * size + i];
            for (int k = i + 1; k < size; k++)
                a[j * size + k] -= a[j * size + i] * a[i * size + k];
        }
}

/**
 * LU decomposition in case of square matrix
 */
void CommonLUDecomposer::findDecompositionParallel()
{
    int i, j, k;
    for (i = 0; i < size; i++)
    {
#pragma omp parallel for private(j)
        for (j = i + 1; j < size; j++)
        {
            a[j * size + i] /= a[i * size + i];
#pragma omp parallel for private(k)
            for (k = i + 1; k < size; k++)
            {
                a[j * size + k] -= a[j * size + i] * a[i * size + k];
            }
        }
    }
}

void CommonLUDecomposer::printDecomposition()
{
    cout << "Demmel Lu decomposition" << endl;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            cout << a[i * size + j] << " ";
        }
        cout << endl;
    }
}
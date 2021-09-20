#include "header.hpp"

DemmelLuSolver::DemmelLuSolver(int sizeParamN, int sizeParamM)
{
    sizeN = sizeParamN;
    sizeM = sizeParamM;
    a = new double[sizeM * sizeN];
}

DemmelLuSolver::~DemmelLuSolver()
{
    delete[] a;
}

void DemmelLuSolver::setA(double *x)
{
    for (int i = 0; i < sizeM * sizeN; i++)
        a[i] = x[i];
}

/**
 * LU decomposition in case of square matrix
 */ 
void DemmelLuSolver::findDecompositionSquare()
{
    for (int i = 0; i < sizeM; i++)
        for (int j = i + 1; j < sizeM; j++)
        {
            a[j * sizeM + i] /= a[i * sizeM + i];

            for (int k = i + 1; k < sizeM; k++)
                a[j * sizeM + k] -= a[j * sizeM + i] * a[i * sizeM + k];
        }
}

/**
 * LU decomposition. Straight algorithm
 */
void DemmelLuSolver::findDecomposition()
{
    int sizeX = std::min(sizeM - 1, sizeN);

    for (int i = 0; i < sizeX; i++)
    {
        for (int j = i + 1; j < sizeM; j++)
        {
            a[j * sizeM + i] /= a[i * sizeM + i];
        }

        if (i < sizeN)
        {
            for (int k = i + 1; k < sizeX; k++)
                for (int j = i + 1; j < sizeN; j++)
                    a[j * sizeM + k] -= a[j * sizeM + i] * a[i * sizeM + k];
        }
    }
}

/**
 * Straight LU decomposition for block, whick starts at (iStart, jStart) and ends with (iEnd,jEnd)
 */
void DemmelLuSolver::findDecomposition(int iStart, int jStart, int iEnd, int jEnd)
{
    int currentSize = jEnd - jStart;
    for (int i = iStart; i < iEnd; i++)
        for (int j = i + 1; j < jEnd; j++)
        {
            a[j * currentSize + i] /= a[i * currentSize + i];

            for (int k = i + 1; k < iEnd; k++)
                a[j * currentSize + k] -= a[j * currentSize + i] * a[i * currentSize + k];
        }
}

void DemmelLuSolver::printDecomposition()
{
    cout << "Demmel Lu decomposition" << endl;
    for (int i = 0; i < sizeM; i++)
    {
        for (int j = 0; j < sizeN; j++)
        {
            cout << a[i * sizeM + j] << " ";
        }
        cout << endl;
    }
}
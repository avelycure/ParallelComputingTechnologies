#include "header.hpp"

DemmelLuSolver::DemmelLuSolver(int sizeParam)
{
    size = sizeParam;
    a = new double[size * size];
}

DemmelLuSolver::~DemmelLuSolver()
{
    delete[] a;
}

void DemmelLuSolver::setA(double *x)
{
    for (int i = 0; i < size * size; i++)
        a[i] = x[i];
}

void DemmelLuSolver::findDecomposition()
{
    for (int i = 0; i < size; i++)
        for (int j = i + 1; j < size; j++)
        {
            a[j * size + i] /= a[i * size + i];

            for (int k = i + 1; k < size; k++)
                a[j * size + k] -= a[j * size + i] * a[i * size + k];
        }
}

void DemmelLuSolver::printDecomposition()
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
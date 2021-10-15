#include "../header.hpp"

void fillMatrixRandom(double *a, int size)
{
    srand(time(0));
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            a[i * size + j] = rand() % 1000;
}

void fillVectorRandom(double *x, int size)
{
    srand(time(0));
    for (int i = 0; i < size; i++)
        x[i] = rand() % 1000;
}

void printVector(double *x, int size, std::string name)
{
    cout << "Vector " + name << endl;
    for (int i = 0; i < size; i++)
        cout << x[i] << " ";
    cout << endl;
}

double findNorm(double *x, double *y, int size)
{
    double max = 0.0;
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            if (fabs(x[i * size + j] - y[i * size + j]) > max)
                max = fabs(x[i * size + j] - y[i * size + j]);
    return max;
}

void copyMatrix(double *x, double *y, int size)
{
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            y[i * size + j] = x[i * size + j];
}

double *makeLFromA(double *a, int size)
{
    double *l = new double[size * size];
    for (int i = 0; i < size; i++)
    {
        l[i * size + i] = 1.0;
        for (int j = 0; j < i; j++)
            l[i * size + j] = a[i * size + j];
    }
    return l;
}

double *makeUFromA(double *a, int size)
{
    double *u = new double[size * size];
    for (int i = 0; i < size; i++)
    {
        for (int j = i; j < size; j++)
            u[i * size + j] = a[i * size + j];
    }
    return u;
}

double *matr_product(double *A, double *B, int size)
{
    double *C = new double[size * size];
    for (int i = 0; i < size; ++i)
    {
        for (int k = 0; k < size; ++k)
        {
            for (int j = 0; j < size; ++j)
            {
                C[i * size + j] += A[i * size + k] * B[k * size + j];
            }
        }
    }
    return C;
}
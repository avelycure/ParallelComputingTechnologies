#include "header.hpp"

void readFromFile(double *a, string fileName)
{
    ifstream indata;
    //int num;
    int i = 0;
    indata.open(fileName);
    if (!indata)
    {
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    //indata >> num;
    while (!indata.eof())
    {
        indata >> a[i];
        i++;
    }
    indata.close();
}

void fillMatrixRandom(double *a, int sizeM, int sizeN)
{
    srand(time(0));
    for (int i = 0; i < sizeM; i++)
        for (int j = 0; j < sizeN; j++)
            a[i * sizeM + j] = ((double)rand()) / rand();
}

void fillVectorRandom(double *x, int size)
{
    srand(time(0));
    for (int i = 0; i < size; i++)
        x[i] = ((double)rand()) / rand();
}

void printVector(double *x, int size, std::string name)
{
    cout << "Vector " + name << endl;
    for (int i = 0; i < size; i++)
        cout << x[i] << " ";
    cout << endl;
}

void read_matr(string namefile, double *&matr, int n, int m)
{
    ifstream fin;
    fin.open(namefile.c_str(), std::ios::in);
    if (!fin.is_open())
    {
        std::cout << "Could not open file!\n";
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                fin >> matr[i * n + j];
            }
        }
    }
    fin.close();
}

void write_matr(string namefile, double *&matr, int n, int m)
{
    ofstream fout;
    fout.open(namefile.c_str());
    if (!fout.is_open())
    {
        std::cout << "Could not open file!\n";
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                fout << matr[i * n + j] << "\t";
            }
            fout << "\n";
        }
    }
    fout.close();
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

double * matr_product(double *A, double *B, int size)
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
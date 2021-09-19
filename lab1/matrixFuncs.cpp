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
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
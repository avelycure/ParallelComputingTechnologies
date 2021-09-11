#include "header.hpp"

int main()
{
    int size = 3;
    double* a = new double[size*size];
    a[0] = 2.0;
    a[1] = -1.0;
    a[2] = -2.0;

    a[3] = -4.0;
    a[4] = 6.0;
    a[5] = 3.0;

    a[6] = -4.0;
    a[7] = -2.0;
    a[8] = 8.0;
    LUSolver luSolver = LUSolver(size);

    luSolver.findDecomposition(a);
    luSolver.printDecomposition();
}

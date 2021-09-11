#include "header.hpp"

int main()
{
    int size = 3;
    double* a = new double[size*size];
    readFromFile(a, "A.txt");

    LUSolver luSolver = LUSolver(size);

    luSolver.findDecomposition(a);
    luSolver.printDecomposition();
}

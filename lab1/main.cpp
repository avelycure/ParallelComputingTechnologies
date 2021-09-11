#include "header.hpp"

int main()
{
    int size = 3;
    double* a = new double[size*size];
    //readFromFile(a, "A.txt");
    fillRandom(a, size);

    LUSolver luSolver = LUSolver(size);

    luSolver.findDecomposition(a);
    luSolver.printDecomposition();
}

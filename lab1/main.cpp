#include "header.hpp"

int main()
{
    int size = 3;
    double *a = new double[size * size];
    double *b = new double[size];
    readFromFile(a, "input_data/a.txt");
    readFromFile(b, "input_data/b.txt");
    //fillMatrixRandom(a, size);
    //fillVectorRandom(b, size);

    LUSolver luSolver = LUSolver(size);

    luSolver.findDecomposition(a);
    luSolver.printDecomposition();

    double *res = luSolver.solveWith(b);
    printVector(res, size, "Result");
}

#include "header.hpp"

int main()
{
    int sizeN = 5;
    int sizeM = 5;
    double *a = new double[sizeM * sizeN];
    double *b = new double[sizeM];
    //readFromFile(a, "input_data/a.txt");
    //readFromFile(b, "input_data/b.txt");
    fillMatrixRandom(a, sizeM, sizeN);
    //fillVectorRandom(b, size);

    DemmelLuSolver demmelLuSolver = DemmelLuSolver(sizeN, sizeM);
    demmelLuSolver.setA(a);
    //demmelLuSolver.findDecomposition(0, 0, 4, 4);
    demmelLuSolver.findDecomposition();
    demmelLuSolver.printDecomposition();

    LUSolver luSolver = LUSolver(sizeN);

    luSolver.findDecomposition(a);
    luSolver.printDecomposition();

    //double *res = luSolver.solveWith(b);
    //printVector(res, size, "Result");
}

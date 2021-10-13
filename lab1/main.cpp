#include "header.hpp"

/**
 * Program which search lu decomposition of the matrix
 */
int main()
{
    //init
    int sizeN = 16;
    int sizeM = 16;

    double *a = new double[sizeM * sizeN];
    double *aCopy = new double[sizeM * sizeN];
    double *b = new double[sizeM];

    //readFromFile(a, "input_data/a.txt");
    //readFromFile(b, "input_data/b.txt");

    fillMatrixRandom(a, sizeM, sizeN);
    copyMatrix(a, aCopy, sizeM);

    //fillVectorRandom(b, size);

    //Decompose

    //Part 1. Common LU decomposition. Demmel variant.
    DemmelLuSolver demmelLuSolver = DemmelLuSolver(sizeN, sizeM);
    demmelLuSolver.setA(a);
    //demmelLuSolver.findDecomposition();
    demmelLuSolver.findDecompositionSquare();//square variant
    // demmelLuSolver.findDecomposition(0, 0, 4, 4);//block with size of matrix
    // demmelLuSolver.printDecomposition();

    LUBlocked luBlocked = LUBlocked(16, 4);
    luBlocked.block_LU_decomposition(aCopy, 16, 4);

    double *l1 = makeLFromA(aCopy, sizeM);
    double *u1 = makeUFromA(aCopy, sizeM);

    double *l2 = makeLFromA(demmelLuSolver.getMatrix(), sizeM);
    double *u2 = makeUFromA(demmelLuSolver.getMatrix(), sizeM);

    cout << "Norm0: " << findNorm(matr_product(l1, u1, sizeM), a, sizeM) << endl;
    cout << "Norm1: " << findNorm(matr_product(l2, u2, sizeM), a, sizeM) << endl;

    cout << "Norm3: " << findNorm(demmelLuSolver.getMatrix(), aCopy, sizeM) << endl;

    //Part 2. LU decomposition. Kalitkin variant.
    // LUSolver luSolver = LUSolver(sizeM);
    // luSolver.findDecomposition(a);
    // luSolver.printDecomposition();

    // double *res = luSolver.solveWith(b);
    // printVector(res, sizeM, "Result");

    //Part 3. Blocked LU decomposition
    /*auto LU_step_by_step = LUBlocked(128 * 16, 512);

    auto begin = std::chrono::steady_clock::now();
    LU_step_by_step.LUDecomposition();
    auto end = std::chrono::steady_clock::now();
    // Mesuring time
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin); 
    std::cout << "elapsed ms: " << elapsed_ms.count() << "\n";
*/
}

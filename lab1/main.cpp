#include "header.hpp"

/**
 * Program which search lu decomposition of the matrix
 */
int main()
{
    omp_set_num_threads(4);

    //init
    int sizeN = 5;
    int sizeM = 5;

    double *a = new double[sizeM * sizeN];
    double *b = new double[sizeM];

    //readFromFile(a, "input_data/a.txt");
    //readFromFile(b, "input_data/b.txt");

    fillMatrixRandom(a, sizeM, sizeN);
    //fillVectorRandom(b, size);

    //Decompose

    //Part 1. Common LU decomposition. Demmel variant.
    // DemmelLuSolver demmelLuSolver = DemmelLuSolver(sizeN, sizeM);
    // demmelLuSolver.setA(a);
    // demmelLuSolver.findDecomposition();
    //demmelLuSolver.findDecompositionSquare();//square variant
    //demmelLuSolver.findDecomposition(0, 0, 4, 4);//block with size of matrix
    // demmelLuSolver.printDecomposition();

    //Part 2. LU decomposition. Kalitkin variant.
    // LUSolver luSolver = LUSolver(sizeM);
    // luSolver.findDecomposition(a);
    // luSolver.printDecomposition();

    // double *res = luSolver.solveWith(b);
    // printVector(res, sizeM, "Result");

    //Part 3. Blocked LU decomposition
    auto LU_step_by_step = LUBlocked(128 * 16, 512);

    auto begin = std::chrono::steady_clock::now();
    LU_step_by_step.LUDecomposition();
    auto end = std::chrono::steady_clock::now();
    // Mesuring time
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin); 
    std::cout << "elapsed ms: " << elapsed_ms.count() << "\n";

}

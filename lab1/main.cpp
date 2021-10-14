#include "header.hpp"

/**
 * Program which search lu decomposition of the matrix
 */
void countError(double* a, double *a1, double *a2, int sizeN, int sizeM);

int main()
{
    double timeStraight;
    double timeBlocked;
    double timeBlockedParallel;
    
    //init
    int sizeN = 512;
    int sizeM = 512;
    int blockSize = 32; //128 best on 1024

    double *a = new double[sizeM * sizeN];
    double *aCopy = new double[sizeM * sizeN];
    double *b = new double[sizeM];

    fillMatrixRandom(a, sizeM, sizeN);
    copyMatrix(a, aCopy, sizeM);

    // Decompose
    // Part 1. Common LU decomposition
    DemmelLuSolver demmelLuSolver = DemmelLuSolver(sizeN, sizeM);
    demmelLuSolver.setA(a);
    auto begin = std::chrono::steady_clock::now();
    demmelLuSolver.findDecompositionSquare();
    auto end = std::chrono::steady_clock::now();
    timeStraight = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

    // Part 2. Blocked decomposition
    LUBlocked luBlocked = LUBlocked(sizeN, blockSize);
    luBlocked.setMatrix(a);
    begin = std::chrono::steady_clock::now();
    luBlocked.decompose();
    end = std::chrono::steady_clock::now();
    timeBlocked = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

    // Counting norm of error to compare results of two methods
    countError(a, luBlocked.getMatrix(), demmelLuSolver.getMatrix(), sizeN, sizeM);

    // Part 3. Parallel LU decomposition straight
    luBlocked.setMatrix(a);
    begin = std::chrono::steady_clock::now();
    luBlocked.decomposeParallel(4);
    end = std::chrono::steady_clock::now();
    timeBlockedParallel = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    countError(a, luBlocked.getMatrix(), demmelLuSolver.getMatrix(), sizeN, sizeM);

    cout << "Straight: " << timeStraight << endl;
    cout << "Blocked: " << timeBlocked << endl;
    cout << "BlockedParallel: " << timeBlockedParallel << endl;
}

void countError(double* a, double *a1, double *a2, int sizeN, int sizeM)
{
    double *l1 = makeLFromA(a1, sizeM);
    double *u1 = makeUFromA(a1, sizeM);

    double *l2 = makeLFromA(a2, sizeM);
    double *u2 = makeUFromA(a2, sizeM);

    cout << "Norm0: " << findNorm(matr_product(l1, u1, sizeM), a, sizeM) << endl;
    cout << "Norm1: " << findNorm(matr_product(l2, u2, sizeM), a, sizeM) << endl;
    cout << "Norm2: " << findNorm(a1, a2, sizeM) << endl;
}
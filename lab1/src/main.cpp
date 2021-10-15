#include "../header.hpp"

/**
 * Program which search lu decomposition of the matrix
 */
void countError(double *a, double *a1, double *a2, int size);

int main()
{
    double timeCommonSequatial;
    double timeBlockedSequantial;
    double timeCommonParallel;
    double timeBlockedParallel;

    //init
    int size = 512;
    int blockSize = 32; //128 best on 1024
    int numThreads = 8;

    double *a = new double[size * size];
    double *aCopy = new double[size * size];
    double *b = new double[size];

    fillMatrixRandom(a, size);
    copyMatrix(a, aCopy, size);

    // Decompose
    // Part 1. Common LU decomposition sequentially
    CommonLUDecomposer commonLuDecomposer = CommonLUDecomposer(size);
    commonLuDecomposer.setMatrix(a);
    auto begin = std::chrono::steady_clock::now();
    commonLuDecomposer.findDecomposition();
    auto end = std::chrono::steady_clock::now();
    timeCommonSequatial = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

    // Part 2. Blocked decomposition sequentially
    BlockedLUDecomposer blockedLuDecomposer = BlockedLUDecomposer(size, blockSize);
    blockedLuDecomposer.setMatrix(a);
    begin = std::chrono::steady_clock::now();
    blockedLuDecomposer.findDecomposition();
    end = std::chrono::steady_clock::now();
    timeBlockedSequantial = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

    // Counting norm of error to compare results of two methods
    //countError(a, blockedLuDecomposer.getMatrix(), commonLuDecomposer.getMatrix(), size);

    // Part 3. Common LU decomposition parallel
    commonLuDecomposer.setMatrix(a);
    begin = std::chrono::steady_clock::now();
    commonLuDecomposer.findDecompositionParallel();
    end = std::chrono::steady_clock::now();
    timeCommonParallel = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

    // Part 4. Blocked LU decomposition parallel
    blockedLuDecomposer.setMatrix(a);
    begin = std::chrono::steady_clock::now();
    blockedLuDecomposer.findDecompositionParallel(numThreads);
    end = std::chrono::steady_clock::now();
    timeBlockedParallel = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    countError(a, blockedLuDecomposer.getMatrix(), commonLuDecomposer.getMatrix(), size);

    cout << "CommonSequantial: " << timeCommonSequatial << endl;
    cout << "BlockedSequantial: " << timeBlockedSequantial << endl;
    cout << "CommonParallel: " << timeCommonParallel << endl;
    cout << "BlockedParallel: " << timeBlockedParallel << endl;
}

void countError(double *a, double *a1, double *a2, int size)
{
    double *l1 = makeLFromA(a1, size);
    double *u1 = makeUFromA(a1, size);

    double *l2 = makeLFromA(a2, size);
    double *u2 = makeUFromA(a2, size);

    cout << "Norm0: " << findNorm(matr_product(l1, u1, size), a, size) << endl;
    cout << "Norm1: " << findNorm(matr_product(l2, u2, size), a, size) << endl;
    cout << "Norm2: " << findNorm(a1, a2, size) << endl;
}
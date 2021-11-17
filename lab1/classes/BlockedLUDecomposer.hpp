class BlockedLUDecomposer
{
private:
    double *matrix;
    int matrixSize;
    int blockSize;
    void clearMemory(double *obj);

public:
    BlockedLUDecomposer(int matrixSize, int blockSize);
    ~BlockedLUDecomposer();

    void LUDecomposition();

    void findDecomposition();
    void findDecompositionParallel(int numTh);

    void setMatrix(double *x);
    double *getMatrix() { return matrix; };
    int getSize() { return matrixSize; };
    int getBlockSize() { return blockSize; };
};
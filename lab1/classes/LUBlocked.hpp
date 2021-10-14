class LUBlocked
{
public:
    LUBlocked(int matrixSize, int blockSize);
    ~LUBlocked();
    void LUDecomposition();
    void decompose();
    void decomposeParallel(int numTh);

    void setMatrix(double *x);

    double *getMatrix(){return matrix;};
    int getSize(){return matrixSize;};
    int getBlockSize(){return blockSize;};

private:
    double *matrix;
    int matrixSize;
    int blockSize;
    void clearMemory(double *obj);
};
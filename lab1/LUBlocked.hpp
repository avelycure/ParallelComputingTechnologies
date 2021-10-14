class LUBlocked
{
public:
    double *matrix;
    int matrixSize;
    int blockSize;
    int numBlocks;

    LUBlocked(int matrixSize_, int blockSize_);
    ~LUBlocked();
    void LUDecomposition();
    void decompose();
    void decomposeParallel(int numTh);
    
    void setMatrix(double *x);
    double* getMatrix()
    {
        return matrix;
    };

private:
    void A22(int);
    void A23(int, int);
    void A32(int, int);
    void A33(int, int, int);

    void clearMemory(double *obj);
};
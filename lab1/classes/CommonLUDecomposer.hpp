class CommonLUDecomposer
{
private:
    double *a;
    int size;
public:
    CommonLUDecomposer(int sizeParam);
    ~CommonLUDecomposer();

    void findDecomposition();
    void findDecompositionParallel(int numTh);
    int getSize() { return size; };
    void setMatrix(double *a);
    double getMatrixElement(int i, int j) { return a[i * size + j]; };
    double *getMatrix() { return a; };

    void printDecomposition();
};
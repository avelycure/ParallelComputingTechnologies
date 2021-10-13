class DemmelLuSolver
{
private:
    double *a;
    int sizeM;
    int sizeN;

    void decompose(double *a);
    void printMatrix(double *&x, std::string matrixName);

public:
    int getSizeM() { return sizeM; };
    int getSizeN() { return sizeN; };
    DemmelLuSolver(int sizeParamN, int sizeParamM);
    ~DemmelLuSolver();
    void findDecompositionSquare();
    void findDecomposition();
    void findDecomposition(int iStart, int jStart, int iEnd, int jEnd);
    void printDecomposition();
    void setA(double *a);
    //todo refactor class to work with square matrixes
    double getMatrixElement(int i, int j)
    {
        return a[i * sizeM + j];
    };

    double* getMatrix()
    {
        return a;
    };
};
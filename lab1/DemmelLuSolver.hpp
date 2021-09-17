class DemmelLuSolver
{
private:
    double *a;
    int size;

    void decompose(double *A);
    void printMatrix(double *&x, std::string matrixName);
public:
    int getSize() { return size; };
    DemmelLuSolver(int size);
    ~DemmelLuSolver();
    void findDecomposition();
    void printDecomposition();
    void setA(double* a);
};
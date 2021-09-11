class LUSolver
{
private:
    double *l;
    double *u;
    double *m;
    int size;

    int getSize() { return size; };

    void decompose(double *A);

    void multiply(double *&x, double *&y, double *&res);
    void resetMatrix(double *&x);
    void setDiagonal(double *&x);
    void printMatrix(double *&x, std::string matrixName);

public:
    LUSolver(int size);
    void findDecomposition(double *&a);
    void printDecomposition();
};
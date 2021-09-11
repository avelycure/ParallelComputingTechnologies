class LUSolver
{
private:
    double *l;
    double *u;
    double *x;
    double *m;
    int size;

    int getSize() { return size; };

    LUSolver(int size);
    void decompose(double *A);

    void multiply(double *&x, double *&y, double *&res);
    void resetMatrix(double *&x);
    void setDiagonal(double *&x);
    void printMatrix(double *&x, string matrixName);

public:
    void findDecomposition(double *&a);
    void printDecomposition();
};
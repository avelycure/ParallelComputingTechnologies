class LUSolver
{
private:
    double *l;
    double *u;
    double *m;
    int size;

    void decompose(double *A);

    void multiply(double *&x, double *&y, double *&res);
    void resetMatrix(double *&x);
    void setDiagonal(double *&x);
    void printMatrix(double *&x, std::string matrixName);

    void forwardSubstitution(double *&b, double *&y);
    void backwordSubstitution(double *&b, double *&x, double *&y);

public:
    int getSize() { return size; };
    LUSolver(int size);
    ~LUSolver();
    void findDecomposition(double *&a);
    double *solveWith(double *&b);
    void printDecomposition();
};
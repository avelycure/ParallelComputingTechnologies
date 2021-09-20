class LUBlocked
{
public:
    double *matrix;
    size_t matrix_size;
    size_t block_size;
    size_t num_blocks;
    
    LUBlocked(size_t matrix_size_, size_t block_size_);
    ~LUBlocked();
    void LUDecomposition();

private:
    void A22(int);
    void A23(int, int);
    void A32(int, int);
    void A33(int, int, int);
};
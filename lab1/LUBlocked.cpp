#include "header.hpp"

LUBlocked::LUBlocked(size_t matrix_size_, size_t block_size_)
{
    srand(time(0));
    matrix_size = matrix_size_;
    block_size = block_size_;
    num_blocks = matrix_size / block_size;

    matrix = new double[matrix_size * matrix_size];

    for (size_t i = 0; i < matrix_size; ++i)
    {
        for (size_t j = 0; j < matrix_size; ++j)
        {
            matrix[(i * matrix_size) + j] = rand() % 200 + 2;
        }
    }
}

LUBlocked::~LUBlocked()
{
    delete[] matrix;
}

void LUBlocked::A22(int ii)
{
    for (size_t i = ii * block_size; i < (ii * block_size) + block_size - 1; i++)
    {
        for (size_t j = i + 1; j < (ii * block_size) + block_size; j++)
        {
            matrix[(j * matrix_size) + i] /= matrix[(i * matrix_size) + i];

            for (size_t k = i + 1; k < (ii * block_size) + block_size; k++)
            {
                matrix[(j * matrix_size) + k] = matrix[(j * matrix_size) + k] - (matrix[(j * matrix_size) + i] * matrix[(i * matrix_size) + k]);
            }
        }
    }
}

void LUBlocked::A23(int ii, int jj)
{
    for (size_t i = ii * block_size; i < (ii * block_size) + (block_size - 1); i++)
    {
        for (size_t j = i + 1; j < block_size; j++)
        {
            for (size_t k = jj * block_size; k < (jj * block_size) + block_size; k++)
            {
                matrix[(j * matrix_size) + k] = matrix[(j * matrix_size) + k] - (matrix[(j * matrix_size) + i] * matrix[(i * matrix_size) + k]);
            }
        }
    }
}

void LUBlocked::A32(int II, int JJ)
{
    for (size_t i = II * block_size; i < (II * block_size) + block_size; i++)
    {
        for (size_t j = JJ * block_size; j < (JJ * block_size) + block_size; j++)
        {
            matrix[(j * matrix_size) + i] /= matrix[(i * matrix_size) + i];

            for (size_t k = i + 1; k < (II * block_size) + block_size; k++)
            {
                matrix[(j * matrix_size) + k] = matrix[(j * matrix_size) + k] - (matrix[(j * matrix_size) + i] * matrix[(i * matrix_size) + k]);
            }
        }
    }
}

void LUBlocked::A33(int ii, int jj, int kk)
{
    for (size_t i = ii * block_size; i < (ii * block_size) + block_size; i++)
    {
        for (size_t j = jj * block_size; j < (jj * block_size) + block_size; j++)
        {
            for (size_t k = kk * block_size; k < (kk * block_size) + block_size; k++)
            {
                matrix[(j * matrix_size) + k] = matrix[(j * matrix_size) + k] - (matrix[(j * matrix_size) + i] * matrix[(i * matrix_size) + k]);
            }
        }
    }
}

void LUBlocked::LUDecomposition()
{
    size_t i, j, k;
    for (i = 0; i < num_blocks; i++)
    {

        A22(i);

        #pragma omp parallel for private(j)
        for (j = i + 1; j < num_blocks; ++j)
        {
            A23(i, j);
        }

        #pragma omp parallel for private(j)
        for (j = i + 1; j < num_blocks; ++j)
        {
            A32(i, j);
        }

        #pragma omp parallel for private(j, k)
        for (j = i + 1; j < num_blocks; ++j)
        {
            for (k = i + 1; k < num_blocks; k++)
            {
                A33(i, j, k);
            }
        }
    }
}
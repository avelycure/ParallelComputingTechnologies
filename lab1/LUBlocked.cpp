#include "header.hpp"

void LUBlocked::setMatrix(double *x)
{
    for (int i = 0; i < matrixSize * matrixSize; i++)
        matrix[i] = x[i];
}

LUBlocked::LUBlocked(int matrixSize_, int blockSize_)
{
    srand(time(0));
    matrixSize = matrixSize_;
    blockSize = blockSize_;
    numBlocks = matrixSize / blockSize;

    matrix = new double[matrixSize * matrixSize];
}

LUBlocked::~LUBlocked()
{
    delete[] matrix;
}

void LUBlocked::A22(int ii)
{
    for (size_t i = ii * blockSize; i < (ii * blockSize) + blockSize - 1; i++)
    {
        for (size_t j = i + 1; j < (ii * blockSize) + blockSize; j++)
        {
            matrix[(j * matrixSize) + i] /= matrix[(i * matrixSize) + i];

            for (size_t k = i + 1; k < (ii * blockSize) + blockSize; k++)
            {
                matrix[(j * matrixSize) + k] = matrix[(j * matrixSize) + k] - (matrix[(j * matrixSize) + i] * matrix[(i * matrixSize) + k]);
            }
        }
    }
}

void LUBlocked::A23(int ii, int jj)
{
    for (size_t i = ii * blockSize; i < (ii * blockSize) + (blockSize - 1); i++)
    {
        for (size_t j = i + 1; j < blockSize; j++)
        {
            for (size_t k = jj * blockSize; k < (jj * blockSize) + blockSize; k++)
            {
                matrix[(j * matrixSize) + k] = matrix[(j * matrixSize) + k] - (matrix[(j * matrixSize) + i] * matrix[(i * matrixSize) + k]);
            }
        }
    }
}

void LUBlocked::A32(int II, int JJ)
{
    for (size_t i = II * blockSize; i < (II * blockSize) + blockSize; i++)
    {
        for (size_t j = JJ * blockSize; j < (JJ * blockSize) + blockSize; j++)
        {
            matrix[(j * matrixSize) + i] /= matrix[(i * matrixSize) + i];

            for (size_t k = i + 1; k < (II * blockSize) + blockSize; k++)
            {
                matrix[(j * matrixSize) + k] = matrix[(j * matrixSize) + k] - (matrix[(j * matrixSize) + i] * matrix[(i * matrixSize) + k]);
            }
        }
    }
}

void LUBlocked::A33(int ii, int jj, int kk)
{
    for (size_t i = ii * blockSize; i < (ii * blockSize) + blockSize; i++)
    {
        for (size_t j = jj * blockSize; j < (jj * blockSize) + blockSize; j++)
        {
            for (size_t k = kk * blockSize; k < (kk * blockSize) + blockSize; k++)
            {
                matrix[(j * matrixSize) + k] = matrix[(j * matrixSize) + k] - (matrix[(j * matrixSize) + i] * matrix[(i * matrixSize) + k]);
            }
        }
    }
}

void LUBlocked::LUDecomposition()
{
    size_t i, j, k;
    for (i = 0; i < numBlocks; i++)
    {
        A22(i);

#pragma omp parallel for private(j)
        for (j = i + 1; j < numBlocks; j++)
        {
            A23(i, j);
        }

#pragma omp parallel for private(j)
        for (j = i + 1; j < numBlocks; j++)
        {
            A32(i, j);
        }

#pragma omp parallel for private(j, k)
        for (j = i + 1; j < numBlocks; j++)
        {
            for (k = i + 1; k < numBlocks; k++)
            {
                A33(i, j, k);
            }
        }
    }
}

void LUBlocked::decompose()
{
    double *a11 = new double[blockSize * blockSize]; //diagonal block
    double *u12 = new double[blockSize * (matrixSize - blockSize)];
    double *l21 = new double[(matrixSize - blockSize) * blockSize];

    for (int bi = 0; bi < matrixSize - 1; bi += blockSize)
    {
        // variable to increase efficity of work with cache
        double temp;

        // Copy of diagonal block
        for (int i = 0; i < blockSize; i++)
            for (int j = 0; j < blockSize; j++)
                a11[i * blockSize + j] = matrix[(i + bi) * matrixSize + (j + bi)];

        // Copy u12 block
        for (int i = 0; i < matrixSize - bi - blockSize; i++)
            for (int j = 0; j < blockSize; j++)
                u12[i * blockSize + j] = matrix[(j + bi) * matrixSize + (i + bi + blockSize)];

        // Copy l21 block
        for (int i = 0; i < matrixSize - bi - blockSize; i++)
            for (int j = 0; j < blockSize; j++)
                l21[i * blockSize + j] = matrix[(i + bi + blockSize) * matrixSize + (j + bi)];

        // LU decomposition of diagonal block
        for (int i = 0; i < blockSize - 1; i++)
            for (int j = i + 1; j < blockSize; j++)
            {
                temp = a11[j * blockSize + i] / a11[i * blockSize + i];
                for (int k = i + 1; k < blockSize; k++)
                    a11[j * blockSize + k] -= temp * a11[i * blockSize + k];
                a11[j * blockSize + i] = temp;
            }

        // Filling u12 block
        for (int j = 0; j < matrixSize - bi - blockSize; j++)
            for (int i = 1; i < blockSize; i++)
            {
                temp = 0.0;
                for (int k = 0; k <= i - 1; k++)
                    temp += a11[i * blockSize + k] * u12[j * blockSize + k];
                u12[j * blockSize + i] -= temp;
            }

        // Filling l21 block
        for (int i = 0; i < matrixSize - bi - blockSize; i++)
            for (int j = 0; j < blockSize; j++)
            {
                temp = 0.0;
                for (int k = 0; k <= j - 1; k++)
                    temp += l21[i * blockSize + k] * a11[k * blockSize + j];
                l21[i * blockSize + j] = (l21[i * blockSize + j] - temp) / a11[j * blockSize + j];
            }

        // Вычитание перед рекурсией
        for (int i = 0; i < matrixSize - bi - blockSize; i++)
            for (int j = 0; j < matrixSize - bi - blockSize; j++)
            {
                temp = 0.0;
                for (int k = 0; k < blockSize; k++)
                    temp += l21[i * blockSize + k] * u12[j * blockSize + k];
                matrix[(i + bi + blockSize) * matrixSize + (j + bi + blockSize)] -= temp;
            }

        // Copy blocks back to matrix
        // Diagonal
        for (int i = 0; i < blockSize; i++)
            for (int j = 0; j < blockSize; j++)
                matrix[(i + bi) * matrixSize + (j + bi)] = a11[i * blockSize + j];

        // Block u12
        for (int i = 0; i < matrixSize - bi - blockSize; i++)
            for (int j = 0; j < blockSize; j++)
                matrix[(j + bi) * matrixSize + (i + bi + blockSize)] = u12[i * blockSize + j];

        // Block l21
        for (int i = 0; i < matrixSize - bi - blockSize; i++)
            for (int j = 0; j < blockSize; j++)
                matrix[(i + bi + blockSize) * matrixSize + (j + bi)] = l21[i * blockSize + j];
    }

    clearMemory(a11);
    clearMemory(u12);
    clearMemory(l21);
}

void LUBlocked::clearMemory(double *obj)
{
    delete[] obj;
}
#include "../header.hpp"

void BlockedLUDecomposer::setMatrix(double *x)
{
    for (int i = 0; i < matrixSize * matrixSize; i++)
        matrix[i] = x[i];
}

BlockedLUDecomposer::BlockedLUDecomposer(int matrixSize, int blockSize)
{
    this->matrixSize = matrixSize;
    this->blockSize = blockSize;

    matrix = new double[matrixSize * matrixSize];
}

BlockedLUDecomposer::~BlockedLUDecomposer()
{
    delete[] matrix;
}

void BlockedLUDecomposer::findDecomposition()
{
    double *a11 = new double[blockSize * blockSize]; //diagonal block
    double *u12 = new double[blockSize * (matrixSize - blockSize)];
    double *l21 = new double[(matrixSize - blockSize) * blockSize];

    // variable to increase efficity of work with cache
    double temp;

    for (int bi = 0; bi < matrixSize - 1; bi += blockSize)
    {
        // Copy of a11
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
        // Analog of previous cycle
        for (int i = 0; i < matrixSize - bi - blockSize; i++)
            for (int j = 0; j < blockSize; j++)
            {
                temp = 0.0;
                for (int k = 0; k <= j - 1; k++)
                    temp += l21[i * blockSize + k] * a11[k * blockSize + j];
                l21[i * blockSize + j] = (l21[i * blockSize + j] - temp) / a11[j * blockSize + j];
            }

        // Substration
        // Calculate matrix a22~
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

void BlockedLUDecomposer::clearMemory(double *obj)
{
    delete[] obj;
}

void BlockedLUDecomposer::findDecompositionParallel(int numTh)
{
    double *a11 = new double[blockSize * blockSize];
    double *u12 = new double[blockSize * (matrixSize - blockSize)];
    double *l21 = new double[(matrixSize - blockSize) * blockSize];

    double temp;
    for (int bi = 0; bi < matrixSize - 1; bi += blockSize)
    {
        // Copy a11
        for (int i = 0; i < blockSize; i++)
        {
            for (int j = 0; j < blockSize; j++)
            {
                a11[i * blockSize + j] = matrix[(i + bi) * matrixSize + (j + bi)];
            }
        }

        // Copy u12
        for (int i = 0; i < matrixSize - bi - blockSize; i++)
        {
            for (int j = 0; j < blockSize; j++)
            {
                u12[i * blockSize + j] = matrix[(j + bi) * matrixSize + (i + bi + blockSize)];
            }
        }

        // Copy l21
        for (int i = 0; i < matrixSize - bi - blockSize; i++)
        {
            for (int j = 0; j < blockSize; j++)
            {
                l21[i * blockSize + j] = matrix[(i + bi + blockSize) * matrixSize + (j + bi)];
            }
        }

        // LU decomposition of a11
        for (int i = 0; i < blockSize - 1; i++)
        {
#pragma omp parallel for private(temp) num_threads(numTh) if (numTh > 1)
            for (int j = i + 1; j < blockSize; j++)
            {
                temp = a11[j * blockSize + i] / a11[i * blockSize + i];
                for (int k = i + 1; k < blockSize; k++)
                {
                    a11[j * blockSize + k] = a11[j * blockSize + k] - temp * a11[i * blockSize + k];
                }
                a11[j * blockSize + i] = temp;
            }
        }

        // Filling u12
#pragma omp parallel for private(temp) num_threads(numTh) if (numTh > 1)
        for (int j = 0; j < matrixSize - bi - blockSize; j++)
        {
            for (int i = 1; i < blockSize; i++)
            {
                temp = 0.0;
                for (int k = 0; k <= i - 1; k++)
                {
                    temp += a11[i * blockSize + k] * u12[j * blockSize + k];
                }
                u12[j * blockSize + i] -= temp;
            }
        }

        // Filling l21
#pragma omp parallel for private(temp) num_threads(numTh) if (numTh > 1)
        for (int i = 0; i < matrixSize - bi - blockSize; i++)
        {
            for (int j = 0; j < blockSize; j++)
            {
                temp = 0.0;
                for (int k = 0; k <= j - 1; k++)
                {
                    temp += l21[i * blockSize + k] * a11[k * blockSize + j];
                }
                l21[i * blockSize + j] = (l21[i * blockSize + j] - temp) / a11[j * blockSize + j];
            }
        }

        // Computing a22~
#pragma omp parallel for private(temp) num_threads(numTh) if (numTh > 1)
        for (int i = 0; i < matrixSize - bi - blockSize; i++)
        {
            for (int j = 0; j < matrixSize - bi - blockSize; j++)
            {
                temp = 0.0;
                for (int k = 0; k < blockSize; k++)
                {
                    temp += l21[i * blockSize + k] * u12[j * blockSize + k];
                }
                matrix[(i + bi + blockSize) * matrixSize + (j + bi + blockSize)] -= temp;
            }
        }

        // Copy diagonal block
#pragma omp parallel for private(temp) num_threads(numTh) if (numTh > 1)
        for (int i = 0; i < blockSize; i++)
        {
            for (int j = 0; j < blockSize; j++)
            {
                matrix[(i + bi) * matrixSize + (j + bi)] = a11[i * blockSize + j];
            }
        }

        // Copy u12
#pragma omp parallel for private(temp) num_threads(numTh) if (numTh > 1)
        for (int i = 0; i < matrixSize - bi - blockSize; i++)
        {
            for (int j = 0; j < blockSize; j++)
            {
                matrix[(j + bi) * matrixSize + (i + bi + blockSize)] = u12[i * blockSize + j];
            }
        }

        // Copy l21
#pragma omp parallel for private(temp) num_threads(numTh) if (numTh > 1)
        for (int i = 0; i < matrixSize - bi - blockSize; i++)
        {
            for (int j = 0; j < blockSize; j++)
            {
                matrix[(i + bi + blockSize) * matrixSize + (j + bi)] = l21[i * blockSize + j];
            }
        }
    }

    clearMemory(a11);
    clearMemory(u12);
    clearMemory(l21);
}
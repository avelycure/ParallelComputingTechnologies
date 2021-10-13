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

void LUBlocked::block_LU_decomposition(double *&matr, int size, const int bs) {

    // bs  Размер блока
    // Диаг., U_12, L_21 блоки соответственно
    double *A11 = new double[bs * bs];
    double *U12 = new double[bs * (size - bs)];
    double *L21 = new double[(size - bs) * bs];
    
    for (int bi = 0; bi < size - 1; bi += bs) {
        double temp; // Временная переменная
        // Копирование диагонального блока
        for (int i = 0; i < bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                A11[i * bs + j] = matr[(i + bi) * size + (j + bi)];
            }
        }
        // Копирование блока U_12
        for (int i = 0; i < size - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                U12[i * bs + j] = matr[(j + bi) * size + (i + bi + bs)];
            }
        }
        // Копирование блока L_21
        for (int i = 0; i < size - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                L21[i * bs + j] = matr[(i + bi + bs) * size + (j + bi)];
            }
        }
        // LU разложение диагонального блока
        for (int i = 0; i < bs - 1; ++i) {
            for (int j = i + 1; j < bs; ++j) {
                temp = A11[j * bs + i] / A11[i * bs + i];
                for (int k = i + 1; k < bs; ++k) {
                    A11[j * bs + k] = A11[j * bs + k] - temp * A11[i * bs + k];
                }
                A11[j * bs + i] = temp;
            }
        }
        // Заполнение блока U_12
        for (int j = 0; j < size - bi - bs; ++j) {
            for (int i = 1; i < bs; ++i) {
                temp = 0.0;
                for (int k = 0; k <= i - 1; ++k) {
                    temp += A11[i * bs + k] * U12[j * bs + k];
                }
                U12[j * bs + i] -= temp;
            }
        }
        // Заполнение блока L_21
        for (int i = 0; i < size - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                temp = 0.0;
                for (int k = 0; k <= j - 1; ++k) {
                    temp += L21[i * bs + k] * A11[k * bs + j];
                }
                L21[i * bs + j] = (L21[i * bs + j] - temp) / A11[j * bs + j];
            }
        }
        // Вычитание перед рекурсией
        for (int i = 0; i < size - bi - bs; ++i) {
            for (int j = 0; j < size - bi - bs; ++j) {
                temp = 0.0;
                for (int k = 0; k < bs; ++k) {
                     temp += L21[i * bs + k] * U12[j * bs + k];
                }
                matr[(i + bi + bs) * size + (j + bi + bs)] -= temp;
            }
        }
        // Перенос лок. массивов в матрицу
        // Диаг. блок
        for (int i = 0; i < bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                matr[(i + bi) * size + (j + bi)] = A11[i * bs + j];
            }
        }
        // Блок U_12
        for (int i = 0; i < size - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                matr[(j + bi) * size + (i + bi + bs)] = U12[i * bs + j];
            }
        }
        // Блок L_21
        for (int i = 0; i < size - bi - bs; ++i) {
            for (int j = 0; j < bs; ++j) {
                matr[(i + bi + bs) * size + (j + bi)] = L21[i * bs + j];
            }
        }
    }

    clearMemory(A11);
    clearMemory(U12);
    clearMemory(L21);
}

void LUBlocked::clearMemory(double *obj) {
    delete[] obj;
}
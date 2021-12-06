#include "Jacobi.hpp"

/**
 * Solving the linear system of equations with 5 diagonal matrix which we get trying to solve
 * Helmholtz equation in square
 * */
void solveSystem(
    std::vector<double> &yLocal,
    std::vector<double> &yLocalPrevious,
    std::vector<double> &yLocalPreviousUpHighBorder,
    std::vector<double> &yLocalPreviousDownLowBorder,
    int numberOfProcesses,
    int processId,
    int localRows,
    int localSize,
    int localOffsetInRows,
    InitialConditions initialConditions)
{
    //Local renaming to increase readability
    double h = initialConditions.h;
    int n = initialConditions.n;
    //Just coefficient of the equation
    double c = 1.0 / (4.0 + initialConditions.h * initialConditions.h * initialConditions.k * initialConditions.k);

    if (processId != 0)
        for (int j = 1; j < n - 1; j++)
            yLocal[j] = c * (h * h * initialConditions.f(localOffsetInRows * h, j * h) +
                             yLocalPreviousUpHighBorder[j] +
                             yLocalPrevious[n + j] +
                             yLocalPrevious[j - 1] +
                             yLocalPrevious[j + 1]);
    //Else if it is first process do nothing because initial conditions in the first row are zeros

    //Calculate only rows that are not borders of process part
    for (int i = 1; i < localRows - 1; i++)
        for (int j = 1; j < n - 1; j++)
            yLocal[i * n + j] = c * (h * h * initialConditions.f((localOffsetInRows + i) * h, j * h) +
                                     yLocalPrevious[(i - 1) * n + j] +
                                     yLocalPrevious[(i + 1) * n + j] +
                                     yLocalPrevious[i * n + (j - 1)] +
                                     yLocalPrevious[i * n + (j + 1)]);

    if (processId != numberOfProcesses - 1)
        for (int j = 1; j < n - 1; j++)
            yLocal[localSize - n + j] = c * (h * h * initialConditions.f((localOffsetInRows + localRows - 1) * h, j * h) +
                                             yLocalPrevious[localSize - 2 * n + j] +
                                             yLocalPreviousDownLowBorder[j] +
                                             yLocalPrevious[localSize - n + j - 1] +
                                             yLocalPrevious[localSize - n + j + 1]);
}
#include "Seidel.hpp"

//NODE (0,0) IS BLACK

//If offset in rows is even number then that row begins with black
//If offset in rows is odd number then that row begins with red

//(localOffsetInRows % 2 == 0) ? 2 : 1 means that if we are in domain which starts from black node then begin from third element
// (0, 1, 2 in vector order), if we are in domain that begins with red node then we should begin from element on 1 position

/**
 * Method for solving system only in black nodes. This method is using PREVIOUS(on previos iteration) values of our solution
 * */
void solveBlack(
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
        for (int j = (localOffsetInRows % 2 == 0) ? 2 : 1; j < n - 1; j += 2)
            yLocal[j] = c * (h * h * initialConditions.f(localOffsetInRows * h, j * h) +
                             yLocalPreviousUpHighBorder[j] +
                             yLocalPrevious[n + j] +
                             yLocalPrevious[j - 1] +
                             yLocalPrevious[j + 1]);
    //Else if it is first process do nothing because initial conditions in the first row are zeros

    //Calculate only rows that are not borders of process part
    for (int i = 1; i < localRows - 1; i++)
        for (int j = ((localOffsetInRows + i) % 2 == 0) ? 2 : 1; j < n - 1; j += 2)
            yLocal[i * n + j] = c * (h * h * initialConditions.f((localOffsetInRows + i) * h, j * h) +
                                     yLocalPrevious[(i - 1) * n + j] +
                                     yLocalPrevious[(i + 1) * n + j] +
                                     yLocalPrevious[i * n + (j - 1)] +
                                     yLocalPrevious[i * n + (j + 1)]);

    if (processId != numberOfProcesses - 1)
        for (int j = ((localOffsetInRows + localRows - 1) % 2 == 0) ? 2 : 1; j < n - 1; j += 2)
            yLocal[localSize - n + j] = c * (h * h * initialConditions.f((localOffsetInRows + localRows - 1) * h, j * h) +
                                             yLocalPrevious[localSize - 2 * n + j] +
                                             yLocalPreviousDownLowBorder[j] +
                                             yLocalPrevious[localSize - n + j - 1] +
                                             yLocalPrevious[localSize - n + j + 1]);
}

/**
 * Method for solving system only in red nodes. This method is using CURRENT(on current iteration) values of our solution
 * */
void solveRed(
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
        for (int j = (localOffsetInRows % 2 == 0) ? 1 : 2; j < n - 1; j += 2)
            yLocal[j] = c * (h * h * initialConditions.f(localOffsetInRows * h, j * h) +
                             yLocalPreviousUpHighBorder[j] +
                             yLocal[n + j] +
                             yLocal[j - 1] +
                             yLocal[j + 1]);
    //Else if it is first process do nothing because initial conditions in the first row are zeros

    //Calculate only rows that are not borders of process part
    for (int i = 1; i < localRows - 1; i++)
        for (int j = ((localOffsetInRows + i) % 2 == 0) ? 1 : 2; j < n - 1; j += 2)
            yLocal[i * n + j] = c * (h * h * initialConditions.f((localOffsetInRows + i) * h, j * h) +
                                     yLocal[(i - 1) * n + j] +
                                     yLocal[(i + 1) * n + j] +
                                     yLocal[i * n + (j - 1)] +
                                     yLocal[i * n + (j + 1)]);

    if (processId != numberOfProcesses - 1)
        for (int j = ((localOffsetInRows + localRows - 1) % 2 == 0) ? 1 : 2; j < n - 1; j += 2)
            yLocal[localSize - n + j] = c * (h * h * initialConditions.f((localOffsetInRows + localRows - 1) * h, j * h) +
                                             yLocal[localSize - 2 * n + j] +
                                             yLocalPreviousDownLowBorder[j] +
                                             yLocal[localSize - n + j - 1] +
                                             yLocal[localSize - n + j + 1]);
}
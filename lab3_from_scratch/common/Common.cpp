#include "Common.hpp"

/**
 * Divide responsibilities between processes
 * Which should do which work
 * @yLocal is part of solution vector which we search
 * @yLocalPrevious is the same as @yLocal, but on previous iteration
 * */
void divideResponsibilities(
    std::vector<double> &y,
    std::vector<double> &yLocal,
    std::vector<double> &yLocalPrevious,
    int numberOfProcesses,
    int processId,
    int &localSize,
    int &localDisplacement,
    int &localRows,
    int &localOffsetInRows,
    std::vector<int> &processesLocationSizes,
    std::vector<int> &processesDisplacement,
    InitialConditions initialConditions)
{
    //@processesLocalRows stores data about number of rows in part of every process
    std::vector<int> processesLocalRows;
    //@processesOffsetInRows stores data about offset in rows from the beginning of our solution
    std::vector<int> processesOffsetInRows;

    if (processId == 0)
    {
        setUpLocations(
            processesLocationSizes,
            processesDisplacement,
            processesLocalRows,
            processesOffsetInRows,
            numberOfProcesses,
            initialConditions);
    }

    MPI_Scatter(processesLocationSizes.data(), 1, MPI_INT, &localSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(processesDisplacement.data(), 1, MPI_INT, &localDisplacement, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(processesLocalRows.data(), 1, MPI_INT, &localRows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(processesOffsetInRows.data(), 1, MPI_INT, &localOffsetInRows, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

/**
 * Init variables before calculation
 * */
void init(
    int processId,
    std::vector<double> &y,
    std::vector<double> &yLocal,
    std::vector<double> &yLocalPrevious,
    std::vector<double> &yLocalHighBorder,
    std::vector<double> &yLocalLowBorder,
    std::vector<double> &buf1,
    std::vector<double> &buf2,
    InitialConditions initialConditions,
    int localSize)
{
    if (processId == 0)
        fillVectorWithZeros(y);

    yLocal.resize(localSize);
    fillVectorWithZeros(yLocal);

    yLocalPrevious.resize(localSize);
    fillVectorWithZeros(yLocalPrevious);

    yLocalHighBorder.resize(initialConditions.n);
    fillVectorWithZeros(yLocalHighBorder);

    yLocalLowBorder.resize(initialConditions.n);
    fillVectorWithZeros(yLocalLowBorder);

    buf1.resize(initialConditions.n);
    buf2.resize(initialConditions.n);
}

void fillVectorWithZeros(std::vector<double> &x)
{
    for (int i = 0; i < x.size(); i++)
        x[i] = 0.0;
}

void setUpLocations(std::vector<int> &processesLocationSizes,
                    std::vector<int> &processesDisplacement,
                    std::vector<int> &processesLocalRows,
                    std::vector<int> &processesOffsetInRows,
                    int numberOfProcesses,
                    InitialConditions initialConditions)
{
    processesLocationSizes.resize(numberOfProcesses);
    processesDisplacement.resize(numberOfProcesses);
    processesLocalRows.resize(numberOfProcesses);
    processesOffsetInRows.resize(numberOfProcesses);

    //Divide y between processes
    for (int i = 0; i < numberOfProcesses; i++)
        processesLocationSizes[i] = (initialConditions.n / numberOfProcesses) * initialConditions.n;

    //Case of not integer division
    for (int i = 0; i < initialConditions.n - (initialConditions.n / numberOfProcesses) * (numberOfProcesses); i++)
        processesLocationSizes[i] += initialConditions.n;

    //Set up number of rows in every location
    for (int i = 0; i < numberOfProcesses; i++)
        processesLocalRows[i] = processesLocationSizes[i] / initialConditions.n;

    //Calculating offset for every process
    //For example we have 8x8 matrix, 4 processes
    //Then parts will be 0..15, 16..31, 32..47, 48..63
    int processOffset;
    for (int i = 0; i < numberOfProcesses; i++)
    {
        processOffset = 0;
        for (int j = 0; j < i; j++)
            processOffset += processesLocationSizes[j];
        processesDisplacement[i] = processOffset;
    }

    int rowsOffset;
    for (int i = 0; i < numberOfProcesses; i++)
    {
        rowsOffset = 0;
        for (int j = 0; j < i; j++)
            rowsOffset += processesLocalRows[j];
        processesOffsetInRows[i] = rowsOffset;
    }

    printProcessLocations(numberOfProcesses,
                          processesDisplacement,
                          processesLocalRows,
                          processesOffsetInRows,
                          initialConditions.isDebugMode);
}

double infiniteNorm(std::vector<double> &x, std::vector<double> &y)
{
    double norm = 0.0;
    double elem = 0.0;
    for (int i = 0; i < x.size(); i++)
    {
        elem = fabs(x[i] - y[i]);
        if (norm < elem)
            norm = elem;
    }
    return norm;
};
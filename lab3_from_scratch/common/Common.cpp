#include "Common.hpp"

void printLog(
    int numberOfProcesses,
    int sumSize,
    double epsilon,
    bool isDebugMode)
{
    if (isDebugMode)
    {
        std::cout << "\033[1;35m *** Program statistic *** \033[0m" << std::endl;
        std::cout << "Number of processes: " << numberOfProcesses << std::endl;
        std::cout << "Size of sum: " << sumSize << std::endl;
        std::cout << "Epsilon: " << epsilon << std::endl;
    }
}

void printProcessData(
    std::vector<double> yLocal,
    std::vector<double> yLocalPrevious,
    int processId,
    int localSize,
    int localDisplacement,
    int localRows,
    int localOffsetInRows,
    bool isDebugMode)
{
    if (isDebugMode)
    {
        std::cout << "\033[1;33m *** Process statistic *** \033[0m" << std::endl;
        std::cout << "Process id: " << processId << std::endl;
        std::cout << "Local size: " << localSize << std::endl;
        std::cout << "Local rows: " << localRows << std::endl;
        std::cout << "Rows offset: " << localOffsetInRows << std::endl;
        std::cout << "Local displacement: " << localDisplacement << std::endl;
        std::cout << "yLocal size: " << yLocal.size() << std::endl;
        std::cout << "yLocalPrevious size: " << yLocalPrevious.size() << std::endl;
    }
}

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

void printProcessLocations(int numberOfProcesses,
                           std::vector<int> processesDisplacement,
                           std::vector<int> processesLocalRows,
                           std::vector<int> &processesOffsetInRows,
                           bool isDebugMode)
{
    if (isDebugMode)
    {
        for (int i = 0; i < numberOfProcesses; i++)
            std::cout << "Process " << i << ", displacement: " << processesDisplacement[i] << std::endl;

        for (int i = 0; i < numberOfProcesses; i++)
            std::cout << "Process " << i << ", rows: " << processesLocalRows[i] << std::endl;

        for (int i = 0; i < numberOfProcesses; i++)
            std::cout << "Process " << i << ", rows offset: " << processesOffsetInRows[i] << std::endl;
    }
}

void printMethodStatistic(
    std::string methodName,
    int finalNorm,
    int iterationsNumber,
    double timeStart,
    double timeEnd,
    double differenceWithAnalyticSolution,
    bool isDebugMode)
{
    if (isDebugMode)
    {
        std::cout << "\033[1;32m" << methodName << "\033[0m" << std::endl;
        std::cout << "Time: " << timeEnd - timeStart << std::endl;
        //std::cout << "Difference: " << differenceWithAnalyticSolution << std::endl;
        std::cout << "Number of iterations: " << iterationsNumber << std::endl;
    }
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

/**
 * Copy first row of local part of solution
 * This function is needed for transmition
 * */
void copyFirstRow(std::vector<double> &yLocal,
                  std::vector<double> &localHighBorder,
                  InitialConditions initialConditions)
{
    for (int i = 0; i < initialConditions.n; i++)
        localHighBorder[i] = yLocal[i];
}

/**
 * Copy last row of local part of solution
 * This function is needed for transmition
 * */
void copyLastRow(std::vector<double> &yLocal,
                 std::vector<double> &localLowBorder,
                 InitialConditions initialConditions)
{
    for (int i = 0; i < initialConditions.n; i++)
        localLowBorder[i] = yLocal[yLocal.size() - initialConditions.n + i];
}

/**
 * We have to exchange two rows: the first and the last. We have two buffers to simplify readability: @buf1, @buf2.
 * Firstly we send from top of the matrix to the bottom. The schema is like this: 0->1->2->3->... where numbers are processes ranks
 * Then we send data from bottom to the top: 4->3->2->1->0. Each process(except 0 and numberOfProcesses - 1,  which are
 * special cases) has one row before(or up) it and one row after. So we have two vectors to work with: @yLocalPreviousUpHighBorder
 * and @yLocalPreviousDownLowBorder. @yLocalPrevious is vector of solution on current iteration
 * */
void exchangeDataV1(std::vector<double> &yLocalPrevious,
                    std::vector<double> &yLocalPreviousUpHighBorder,
                    std::vector<double> &yLocalPreviousDownLowBorder,
                    std::vector<double> &buf1,
                    std::vector<double> &buf2,
                    int numberOfProcesses,
                    int processId,
                    InitialConditions initialConditions)
{
    //Variables to get transaction status
    MPI_Status statU, statL;

    //Send data from lower rank processes to higher
    if (processId != numberOfProcesses - 1)
    {
        copyLastRow(yLocalPrevious, buf1, initialConditions);
        MPI_Send(buf1.data(), initialConditions.n, MPI_DOUBLE, processId + 1, 800, MPI_COMM_WORLD);
    }
    if (processId != 0)
        MPI_Recv(yLocalPreviousUpHighBorder.data(), initialConditions.n, MPI_DOUBLE, processId - 1, 800, MPI_COMM_WORLD, &statL);

    //Send data from higher rank processes to lower
    if (processId != 0)
    {
        copyFirstRow(yLocalPrevious, buf2, initialConditions);
        MPI_Send(buf2.data(), initialConditions.n, MPI_DOUBLE, processId - 1, 25, MPI_COMM_WORLD);
    }
    if (processId != numberOfProcesses - 1)
        MPI_Recv(yLocalPreviousDownLowBorder.data(), initialConditions.n, MPI_DOUBLE, processId + 1, 25, MPI_COMM_WORLD, &statU);
}
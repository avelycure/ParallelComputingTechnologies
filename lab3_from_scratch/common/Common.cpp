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
    bool isDebugMode)
{
    if (isDebugMode)
    {
        std::cout << "\033[1;33m *** Process statistic *** \033[0m" << std::endl;
        std::cout << "Process id: " << processId << std::endl;
        std::cout << "Local size: " << localSize << std::endl;
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
    InitialConditions initialConditions)
{
    //@processesLocationSizes stores data about number of nodes which belongs to each process
    std::vector<int> processesLocationSizes;
    //@processesLocationSizes stores data about displacement from the beginning of y
    std::vector<int> processesDisplacement;

    if (processId == 0)
    {
        setUpLocations(
            processesLocationSizes,
            processesDisplacement,
            numberOfProcesses,
            initialConditions);
    }

    MPI_Scatter(processesLocationSizes.data(), 1, MPI_INT, &localSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(processesDisplacement.data(), 1, MPI_INT, &localDisplacement, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

/**
 * Init variables before calculation
 * */
void init(
    int processId,
    std::vector<double> &y,
    std::vector<double> &yLocal,
    std::vector<double> &yLocalPrevious,
    int localSize)
{
    if (processId == 0)
        fillVectorWithZeros(y);

    yLocal.resize(localSize);
    yLocalPrevious.resize(localSize);
}

void fillVectorWithZeros(std::vector<double> &x)
{
    for (int i = 0; i < x.size(); i++)
        x[i] = 0.0;
}

void setUpLocations(std::vector<int> &processesLocationSizes,
                    std::vector<int> &processesDisplacement,
                    int numberOfProcesses,
                    InitialConditions initialConditions)
{
    processesLocationSizes.resize(numberOfProcesses);
    processesDisplacement.resize(numberOfProcesses);

    //Divide y between processes
    for (int i = 0; i < numberOfProcesses; i++)
        processesLocationSizes[i] = (initialConditions.n / numberOfProcesses) * initialConditions.n;

    //Case of not integer division
    for (int i = 0; i < initialConditions.n - (initialConditions.n / numberOfProcesses) * (numberOfProcesses); i++)
        processesLocationSizes[i] += initialConditions.n;

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

    //for (int i = 0; i < numberOfProcesses; i++)
    //    std::cout << i << " " << processesDisplacement[i] << std::endl;
}
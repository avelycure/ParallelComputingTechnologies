#include "Common.hpp"

void printLog(
    int numberOfProcesses,
    int sumSize,
    double epsilon)
{
    std::cout << "\033[1;33m *** Program statistic *** \033[0m" << std::endl;
    std::cout << "Number of processes: " << numberOfProcesses << std::endl;
    std::cout << "Size of sum: " << sumSize << std::endl;
    std::cout << "Epsilon = " << epsilon << std::endl;
}

/**
 * Divide responsibilities between processes
 * Which should do which work
 * @yPart is part of solution vector which we search
 * @yPartPrevious is the same as @yPart, but on previous iteration
 * */
void divideResponsibilities(
    std::vector<double> &y,
    std::vector<double> &yPart,
    int numberOfProcesses,
    int processId,
    int &locationSize,
    int &displacement,
    InitialConditions initialConditions)
{
    //@processesLocationSizes stores data about number of nodes which belongs to each process
    std::vector<int> processesLocationSizes;
    //@processesLocationSizes stores data about displacement from the beginning of y
    std::vector<int> processesDisplacement;

    if (processId == 0)
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

        for (int i = 0; i < numberOfProcesses; i++)
            std::cout << i << " " << processesDisplacement[i] << std::endl;
    }

    MPI_Scatter(processesLocationSizes.data(), 1, MPI_INT, &locationSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(processesDisplacement.data(), 1, MPI_INT, &displacement, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
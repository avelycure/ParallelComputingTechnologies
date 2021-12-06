#include "Common.hpp"

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
    bool isDebugMode)
{
    if (isDebugMode)
    {
        std::cout << "\033[1;32m" << methodName << "\033[0m" << std::endl;
        std::cout << "Time: " << timeEnd - timeStart << std::endl;
        std::cout << "Number of iterations: " << iterationsNumber << std::endl;
    }
}

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
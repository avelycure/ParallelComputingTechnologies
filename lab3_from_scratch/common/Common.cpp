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
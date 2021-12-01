#include "Common.hpp"

void printLog(
    int numberOfProcesses,
    int sumSize,
    double khRelation,
    double eps)
{
    std::cout << "LOG" << std::endl;
    std::cout << "numberOfProcesses = " << numberOfProcesses << std::endl;
    std::cout << "sumSize = " << sumSize << std::endl;
    std::cout << "k^2 / h^2 = " << khRelation << std::endl;
    std::cout << "eps = " << eps << std::endl;
}
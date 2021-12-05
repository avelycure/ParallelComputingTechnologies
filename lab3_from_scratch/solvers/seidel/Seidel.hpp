#include <vector>
#include "../../common/Common.hpp"
#include "mpi.h"
#include <iostream>
#include "../../init_conds/InitialConditions.hpp"

void seidelV1(std::vector<double> &y,
              InitialConditions initialConditions,
              int numberOfProcesses,
              int processId);

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
    InitialConditions initialConditions);

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
    InitialConditions initialConditions);
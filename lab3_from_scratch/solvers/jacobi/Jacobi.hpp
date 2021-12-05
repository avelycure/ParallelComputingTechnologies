#include <vector>
#include <iostream>
#include "../../mpi.h"
#include "../../common/Common.hpp"
#include "../../init_conds/InitialConditions.hpp"

void jacobiV1(
    std::vector<double> &y,
    InitialConditions initialConditions,
    int numberOfProcesses,
    int processId);

void copyFirstRow(std::vector<double> &yLocal,
                  std::vector<double> &localHighBorder,
                  InitialConditions initialConditions);

void copyLastRow(std::vector<double> &yLocal,
                 std::vector<double> &localLowBorder,
                 InitialConditions initialConditions);

void exchangeData(std::vector<double> &yLocal,
                  std::vector<double> &yLocalPrevious,
                  std::vector<double> &yLocalPreviousUpHighBorder,
                  std::vector<double> &yLocalPreviousDownLowBorder,
                  std::vector<double> &buf1,
                  std::vector<double> &buf2,
                  int numberOfProcesses,
                  int processId,
                  InitialConditions initialConditions);

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
    InitialConditions initialConditions);
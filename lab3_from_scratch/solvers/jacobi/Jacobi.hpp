#include <vector>
#include <iostream>
#include "mpi.h"
#include "../../common/Common.hpp"
#include "../../init_conds/InitialConditions.hpp"

void jacobiV1(std::vector<double> &y,
              InitialConditions initialConditions,
              int numberOfProcesses,
              int processId);

void jacobiV2(std::vector<double> &y,
              InitialConditions initialConditions,
              int numberOfProcesses,
              int processId);

void jacobiV3(std::vector<double> &y,
              InitialConditions initialConditions,
              int numberOfProcesses,
              int processId);

void solveSystem(std::vector<double> &yLocal,
                 std::vector<double> &yLocalPrevious,
                 std::vector<double> &yLocalPreviousUpHighBorder,
                 std::vector<double> &yLocalPreviousDownLowBorder,
                 int numberOfProcesses,
                 int processId,
                 int localRows,
                 int localSize,
                 int localOffsetInRows,
                 InitialConditions initialConditions);
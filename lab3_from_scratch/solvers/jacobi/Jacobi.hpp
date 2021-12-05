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

void exchangeDataV2(std::vector<double> &yLocal,
                    std::vector<double> &yLocalPrevious,
                    std::vector<double> &yLocalPreviousUpHighBorder,
                    std::vector<double> &yLocalPreviousDownLowBorder,
                    std::vector<double> &buf1,
                    std::vector<double> &buf2,
                    int numberOfProcesses,
                    int processId,
                    InitialConditions initialConditions);

void setSourceAndDestination(const int numberOfProcesses,
                             const int processId,
                             int &higherRankProcess,
                             int &lowerRankProcess);

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
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

void solveSystemV3(std::vector<double> &yLocal,
                   std::vector<double> &yLocalPrevious,
                   std::vector<double> &yLocalPreviousUpHighBorder,
                   std::vector<double> &yLocalPreviousDownLowBorder,
                   int numberOfProcesses,
                   int processId,
                   int localRows,
                   int localSize,
                   int localOffsetInRows,
                   int &lowRequests,
                   int &highRequests,
                   InitialConditions initialConditions,
                   std::vector<double> &buf1,
                   std::vector<double> &buf2,
                   std::vector<MPI_Request> &requestsFromLowerToHigher,
                   std::vector<MPI_Request> &requestsFromHigherToLower,
                   std::vector<MPI_Status> &transactionStateFromTop,
                   std::vector<MPI_Status> &transactionStateFromBottom);

void prepareJacobiRequests(int processId,
                           int numberOfProcesses,
                           InitialConditions initialConditions,
                           int &lowRequests,
                           int &highRequests,
                           std::vector<double> &yLocalPreviousUpHighBorder,
                           std::vector<double> &yLocalPreviousDownLowBorder,
                           std::vector<double> &buf1,
                           std::vector<double> &buf2,
                           std::vector<MPI_Request> &requestsFromLowerToHigher,
                           std::vector<MPI_Request> &requestsFromHigherToLower);
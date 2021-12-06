#include <vector>
#include "../../common/Common.hpp"
#include "mpi.h"
#include <iostream>
#include "../../init_conds/InitialConditions.hpp"

void seidelV1(std::vector<double> &y,
              InitialConditions initialConditions,
              int numberOfProcesses,
              int processId);

void seidelV2(std::vector<double> &y,
              InitialConditions initialConditions,
              int numberOfProcesses,
              int processId);

void seidelV3(std::vector<double> &y,
              InitialConditions initialConditions,
              int numberOfProcesses,
              int processId);

void solveRed(std::vector<double> &yLocal,
              std::vector<double> &yLocalPrevious,
              std::vector<double> &yLocalPreviousUpHighBorder,
              std::vector<double> &yLocalPreviousDownLowBorder,
              int numberOfProcesses,
              int processId,
              int localRows,
              int localSize,
              int localOffsetInRows,
              InitialConditions initialConditions);

void solveBlack(std::vector<double> &yLocal,
                std::vector<double> &yLocalPrevious,
                std::vector<double> &yLocalPreviousUpHighBorder,
                std::vector<double> &yLocalPreviousDownLowBorder,
                int numberOfProcesses,
                int processId,
                int localRows,
                int localSize,
                int localOffsetInRows,
                InitialConditions initialConditions);

void solveBlackV3(int processId,
                  int numberOfProcesses,
                  InitialConditions initialConditions,
                  int &lowRequestsBlack,
                  int &highRequestsBlack,
                  int &localRows,
                  int &localSize,
                  int &localOffsetInRows,
                  std::vector<double> &buf1,
                  std::vector<double> &buf2,
                  std::vector<double> &yLocal,
                  std::vector<double> &yLocalPrevious,
                  std::vector<double> &yLocalPreviousUpHighBorder,
                  std::vector<double> &yLocalPreviousDownLowBorder,
                  std::vector<MPI_Status> &stateFromLowerToHigherBlack,
                  std::vector<MPI_Status> &stateFromHigherToLowerBlack,
                  std::vector<MPI_Request> &requestsFromLowerToHigherBlack,
                  std::vector<MPI_Request> &requestsFromHigherToLowerBlack);

void solveRedV3(int processId,
                int numberOfProcesses,
                InitialConditions initialConditions,
                int &lowRequestsRed,
                int &highRequestsRed,
                int &localRows,
                int &localSize,
                int &localOffsetInRows,
                std::vector<double> &buf1,
                std::vector<double> &buf2,
                std::vector<double> &yLocal,
                std::vector<double> &yLocalPrevious,
                std::vector<double> &yLocalPreviousUpHighBorder,
                std::vector<double> &yLocalPreviousDownLowBorder,
                std::vector<MPI_Status> &stateFromLowerToHigherRed,
                std::vector<MPI_Status> &stateFromHigherToLowerRed,
                std::vector<MPI_Request> &requestsFromLowerToHigherRed,
                std::vector<MPI_Request> &requestsFromHigherToLowerRed);

void prepareRequests(int processId,
                     int numberOfProcesses,
                     InitialConditions initialConditions,
                     int &lowRequestsBlack,
                     int &highRequestsBlack,
                     int &lowRequestsRed,
                     int &highRequestsRed,
                     std::vector<double> &yLocalPreviousUpHighBorder,
                     std::vector<double> &yLocalPreviousDownLowBorder,
                     std::vector<double> &buf1,
                     std::vector<double> &buf2,
                     std::vector<MPI_Request> &requestsFromLowerToHigherBlack,
                     std::vector<MPI_Request> &requestsFromHigherToLowerBlack,
                     std::vector<MPI_Request> &requestsFromLowerToHigherRed,
                     std::vector<MPI_Request> &requestsFromHigherToLowerRed);
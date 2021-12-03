#include <vector>
#include "../../common/Common.hpp"
#include "mpi.h"
#include <iostream>
#include "../../init_conds/InitialConditions.hpp"

void jacobiV1(
    std::vector<double> &y,
    double h,
    int size,
    double kSquare,
    int numberOfProcesses,
    int processId,
    double eps,
    InitialConditions initialConditions,
    std::vector<double> &solution);

void jacobiV2(
    std::vector<double> &y,
    double h,
    int size,
    double kSquare,
    int numberOfProcesses,
    int processId,
    double eps,
    InitialConditions initialConditions,
    std::vector<double> &solution);

void jacobiV3(std::vector<double> &y,
               double h,
               int size,
               double kSquare,
               int numberOfProcesses,
               int processId,
               double eps,
               InitialConditions initialConditions,
               std::vector<double> &solution);
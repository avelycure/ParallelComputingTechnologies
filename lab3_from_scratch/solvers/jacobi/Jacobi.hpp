#include <vector>
#include <iostream>
#include "mpi.h"
#include "../../common/Common.hpp"
#include "../../init_conds/InitialConditions.hpp"


void jacobiV1(
    std::vector<double> &y,
    InitialConditions initialConditions,
    int numberOfProcesses,
    int processId);
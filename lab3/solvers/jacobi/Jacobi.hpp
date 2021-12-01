#include <vector>
#include "../../matrix_funcs/Common.hpp"
#include "../../mpi.h"
#include <iostream>
#include "../../init_conds/InitialConditions.hpp"

void jacobiSendReceive(
    std::vector<double> &y,
    double h,
    int size,
    double k_square,
    int numberOfProcesses,
    int processId,
    double eps,
    double &time,
    InitialConditions initialConditions);
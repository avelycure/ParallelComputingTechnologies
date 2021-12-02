#include <vector>
#include "../../common/Common.hpp"
#include "../../mpi.h"
#include <iostream>
#include "../../init_conds/InitialConditions.hpp"

void seidelV1(std::vector<double> &y,
                        double h,
                        int size,
                        double k_square,
                        int processesNumber,
                        int processId,
                        double eps,
                        InitialConditions initialConditions);

void seidelV2(std::vector<double> &y,
                        double h,
                        int size,
                        double k_square,
                        int processesNumber,
                        int processId,
                        double eps,
                        InitialConditions initialConditions);

void seidelV3(std::vector<double> &y,
                        double h,
                        int size,
                        double k_square,
                        int processesNumber,
                        int processId,
                        double eps,
                        InitialConditions initialConditions);
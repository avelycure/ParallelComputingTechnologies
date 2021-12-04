#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "mpi.h"
#include "../init_conds/InitialConditions.hpp"

void printLog(
    int numberOfProcesses,
    int sumSize,
    double eps);

void divideResponsibilities(
    std::vector<double> &y,
    std::vector<double> &yPart,
    int numberOfProcesses,
    int processId,
    int &locationSize,
    int &displacement,
    InitialConditions initialConditions);

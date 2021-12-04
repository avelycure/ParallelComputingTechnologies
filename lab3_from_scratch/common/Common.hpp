#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "mpi.h"
#include "../init_conds/InitialConditions.hpp"

void printLog(
    int numberOfProcesses,
    int sumSize,
    double epsilon,
    bool isDebugMode);

void divideResponsibilities(
    std::vector<double> &y,
    std::vector<double> &yLocal,
    std::vector<double> &yLocalPrevious,
    int numberOfProcesses,
    int processId,
    int &localSize,
    int &localDisplacement,
    InitialConditions initialConditions);

void fillVectorWithZeros(std::vector<double> &y);

void init(
    int processId,
    std::vector<double> &y,
    std::vector<double> &yLocal,
    std::vector<double> &yLocalPrevious,
    int localSize);

void setUpLocations(std::vector<int> &processesLocationSizes,
                    std::vector<int> &processesDisplacement,
                    int numberOfProcesses,
                    InitialConditions initialConditions);

void printProcessData(
    std::vector<double> yLocal,
    std::vector<double> yLocalPrevious,
    int processId,
    int localSize,
    int localDisplacement,
    bool isDebugMode);
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
    int &localRows,
    int &localOffsetInRows,
    InitialConditions initialConditions);

void printProcessLocations(int numberOfProcesses,
                           std::vector<int> processesDisplacement,
                           std::vector<int> processesLocalRows,
                           std::vector<int> &processesOffsetInRows,
                           bool isDebugMode);

void fillVectorWithZeros(std::vector<double> &y);

void init(
    int processId,
    std::vector<double> &y,
    std::vector<double> &yLocal,
    std::vector<double> &yLocalPrevious,
    int localSize);

void setUpLocations(std::vector<int> &processesLocationSizes,
                    std::vector<int> &processesDisplacement,
                    std::vector<int> &processesLocalRows,
                    std::vector<int> &processesOffsetInRows,
                    int numberOfProcesses,
                    InitialConditions initialConditions);

void printProcessData(
    std::vector<double> yLocal,
    std::vector<double> yLocalPrevious,
    int processId,
    int localSize,
    int localDisplacement,
    int localRows,
    bool isDebugMode);

void printMethodStatistic(
    int finalNorm,
    int iterationsNumber,
    double timeStart,
    double timeEnd,
    double differenceWithAnalyticSolution,
    bool isDebugMode);

double infiniteNorm(std::vector<double> &x, std::vector<double> &y, int begin, int end);
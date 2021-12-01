#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "../../mpi.h"

void printLog(
    int numberOfProcesses,
    int sumSize,
    double khRelation,
    double eps);

void divideVectorBetweenProcesses(std::vector<double> &y,
                                  double h,
                                  int size,
                                  double kSquare,
                                  int processesNumber,
                                  int processId,
                                  std::vector<double> &yPart,
                                  std::vector<double> &yPreviousPart,
                                  std::vector<double> &partOfRightPart,
                                  std::vector<int> &processPartOfData,
                                  std::vector<int> &displacement,
                                  int &vectorPartSize,
                                  int &receiveDisplacement,
                                  int &extraSize,
                                  int &offset);

double infiniteNorm(std::vector<double> &vec1, std::vector<double> &vec2, int begin, int end);

void fillVectorWithZeros(std::vector<double> &y);

void setInteractionsScheme(int processesNumber, int processId, int &dest, int &source, int &send1, int &recv1,
                           int &send2, int &recv2);

void setSourceAndDestination(int processesNumber, int processId, int &destination, int &source);
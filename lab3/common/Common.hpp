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

void initVecLoc(std::vector<double> &y, std::vector<double> &vec_right,
                double h, int size, double k_square, int np, int myid, std::vector<double> &y_loc,
                std::vector<double> &y_loc_prev, std::vector<double> &vec_right_loc, std::vector<int> &len,
                std::vector<int> &disp, int &loc_size, int &recv_disp, int &extr_size, int &offset);

double infiniteNorm(std::vector<double> &vec1, std::vector<double> &vec2, int begin, int end);

void fillVectorWithZeros(std::vector<double> &y);

void send_recv_scheme(const int np, const int myid, int &dest, int &source, int &send1, int &recv1, \
    int &send2, int &recv2);

void setSourceAndDestination(int np, int myid, int &destination, int &source);
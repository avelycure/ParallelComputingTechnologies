#include <iostream>
#include <fstream>
#include <cmath>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <chrono>
#include "LUSolver.hpp"
#include "DemmelLuSolver.hpp"
#include "LUBlocked.hpp"

using namespace std;

void readFromFile(double *a, string fileName);
void fillMatrixRandom(double* a, int sizeM, int sizeN);
void fillVectorRandom(double* x, int sizeM);
void printVector(double *x, int size, std::string name);
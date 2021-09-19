#include <iostream>
#include <fstream>
#include <cmath>
#include <string.h>
#include <fstream>
#include <time.h>
#include "LUSolver.hpp"
#include "DemmelLuSolver.hpp"

using namespace std;

void readFromFile(double *a, string fileName);
void fillMatrixRandom(double* a, int sizeM, int sizeN);
void fillVectorRandom(double* x, int sizeM);
void printVector(double *x, int size, std::string name);
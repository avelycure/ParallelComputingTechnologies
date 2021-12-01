#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <string>
#include "mpi.h"
#include "../init_conds/InitialConditions.hpp"
#include "../solvers/jacobi/Jacobi.hpp"
#include "../solvers/seidel/Seidel.hpp"
#include "../matrix_funcs/Common.hpp"
using namespace std;
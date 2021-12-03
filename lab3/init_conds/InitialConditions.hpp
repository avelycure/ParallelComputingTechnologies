#include <vector>
#pragma once

/**
 * Data class for storing initial information about sum
 * */
class InitialConditions
{
public:
    int n = 128 * 128 - 1;
    double leftBorder = 0.0;
    double rightBorder = 1.0;
    double h = (rightBorder - leftBorder) / n;
    double eps = 1e-8;
    double kSquare = 1.0;

    const double PI = 3.14159265358979323846;
    const double PI_SQUARE = PI * PI;

    /**
     * Right part of equation
     * */
    double f(double x, double y, double kSquare)
    {
        return sin(PI * y) * (2.0 + (kSquare + PI_SQUARE) * (1.0 - x) * x);
    }

    /**
    * Analytic solution of the equation
    * */
    double analyticSolution(double x, double y)
    {
        return (1.0 - x) * x * sin(PI * y);
    }

    void computeSolution(std::vector<double> &solution)
    {
        for (int i = 0; i < n + 1; ++i)
            for (int j = 0; j < n + 1; ++j)
                solution[i * (n + 1) + j] = analyticSolution(i * h, j * h);
    }
};

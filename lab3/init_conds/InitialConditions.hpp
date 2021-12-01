#include <vector>
#pragma once

/**
 * Data class for storing initial information about ssun
 * */
class InitialConditions
{
public:
    int n = 128 * 16;
    double leftBorder = 0.0;
    double rightBorder = 1.0;
    double h = (rightBorder - leftBorder) / n;
    double eps = 1e-8;
    double kSquare = 1.0;

    const double PI = 3.1415;
    const double PI_SQUARE = PI * PI;

    /**
     * Right part of equation
     * */
    double f(double x, double y, double k_square)
    {
        return sin(PI * y) * (2.0 + (k_square + PI_SQUARE) * (1.0 - x) * x);
    }

    double fSolution(double x, double y)
    {
        return (1.0 - x) * x * sin(PI * y);
    }

    void computeSolution(std::vector<double> &solution)
    {
        for (int i = 0; i < n + 1; ++i)
            for (int j = 0; j < n + 1; ++j)
                solution[i * (n + 1) + j] = fSolution(i * h, j * h);
    }
};
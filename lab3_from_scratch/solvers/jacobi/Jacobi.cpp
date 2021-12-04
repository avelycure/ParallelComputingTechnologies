#include "Jacobi.hpp"

/**
 * Used methods: MPI_Send. MPI_Recv
 * */
void jacobiV1(
    std::vector<double> &y,
    InitialConditions initialConditions,
    int numberOfProcesses,
    int processId)
{
    //Number of nodes in the process domain after division between processes
    int locationSize;

    //Displacement in relation with beginning of the result vector
    int displacement;

    //Norm on current iteration
    double currentNorm;
    //Norm on previous iteration
    double previousNorm;

    //Part of the solution, here is situated part with which every process works
    std::vector<double> yPart;

    divideResponsibilities(
        y,
        yPart,
        numberOfProcesses,
        processId,
        locationSize,
        displacement,
        initialConditions);
}
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
    int localSize;

    //Displacement in relation with beginning of the result vector
    int localDisplacement;

    //Norm on current iteration
    double currentNorm;
    //Norm on previous iteration
    double previousNorm;

    //Part of the solution, here is situated part with which every process works
    std::vector<double> yLocal;
    std::vector<double> yLocalPrevious;

    divideResponsibilities(y,
                           yLocal,
                           yLocalPrevious,
                           numberOfProcesses,
                           processId,
                           localSize,
                           localDisplacement,
                           initialConditions);

    init(processId,
         y,
         yLocal,
         yLocalPrevious,
         localSize);

    printProcessData(yLocal,
                     yLocalPrevious,
                     processId,
                     localSize,
                     localDisplacement,
                     initialConditions.isDebugMode);
}
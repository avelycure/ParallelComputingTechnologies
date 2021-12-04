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

    int iterationsNumber = 0;

    //Difference with analytic solution
    double difference = 0.0;

    //Part of the solution, here is situated part with which every process works
    std::vector<double> yLocal;
    std::vector<double> yLocalPrevious;

    //Variables to get transaction status
    MPI_Status statU, statL;

    //Time counters
    double timeStart, timeEnd;

    //Number of rows which belongs to process
    int localRows;

    //Offset in rows from the beginning
    int localOffsetInRows;

    //Just coefficient of the equation
    double c = 1.0 / (4.0 + initialConditions.h * initialConditions.h * initialConditions.k * initialConditions.k);

    divideResponsibilities(y,
                           yLocal,
                           yLocalPrevious,
                           numberOfProcesses,
                           processId,
                           localSize,
                           localDisplacement,
                           localRows,
                           localOffsetInRows,
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
                     localRows,
                     initialConditions.isDebugMode);

    if (processId == 0)
        timeStart = MPI_Wtime();

    for (int i = 0; i < yLocalPrevious.size(); i++)
        yLocalPrevious[i] = i + processId * 10;

    /*std::cout << "Process id: " << processId << " " << yLocalPrevious.back();
    std::cout << " " << &yLocalPrevious.back();
    std::cout << " address: " << &(yLocalPrevious.at(yLocalPrevious.size() - initialConditions.n)) << " ";
    std::cout << " value:" << (*(&(yLocalPrevious.back()) - initialConditions.n + 1)) << std::endl;*/

    /*do
    {
        iterationsNumber++;

        yLocal.swap(yLocalPrevious);

        //Send data from lower rank processes to higher
        MPI_Send(&yLocalPrevious.back() - initialConditions.n + 1, initialConditions.n, MPI_DOUBLE, processId + 1, 42, MPI_COMM_WORLD);
        MPI_Recv(yLocalPrevious.data(), initialConditions.n, MPI_DOUBLE, processId - 1, 42, MPI_COMM_WORLD, &statL);

        //Send data from higher rank processes to lower
        MPI_Send(yLocalPrevious.data(), initialConditions.n, MPI_DOUBLE, processId - 1, 41, MPI_COMM_WORLD);
        MPI_Recv(&yLocalPrevious.back() - initialConditions.n + 1, initialConditions.n, MPI_DOUBLE, processId + 1, 41, MPI_COMM_WORLD, &statU);

        for (int i = 0; i < localRows; i++)
            for (int j = 0; j < initialConditions.n; j++)
                yLocal[i * initialConditions.n + j] = c * (initialConditions.h * initialConditions.h *
                                                               initialConditions.f((i +) * initialConditions.h, j * initialConditions.h) +
                                                           yLocalPrevious[(i - 1) * initialConditions.n + j] +
                                                           yLocalPrevious[(i + 1) * initialConditions.n + j] +
                                                           yLocalPrevious[i * initialConditions.n + (j - 1)] +
                                                           yLocalPrevious[i * initialConditions.n + (j + 1)]);

        norm = infiniteNorm(yPart, yPreviousPart, receiveDisplacement, locationSize - size);

    } while (currentNorm > initialConditions.epsilon);*/

    if (processId == 0)
        timeEnd = MPI_Wtime();

    if (processId == 0)
        printMethodStatistic(currentNorm,
                             iterationsNumber,
                             timeStart,
                             timeEnd,
                             difference,
                             initialConditions.isDebugMode);
}
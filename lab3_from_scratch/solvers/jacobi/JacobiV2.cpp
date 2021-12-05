#include "Jacobi.hpp"

/**
 * Used methods: MPI_Sendrecv
 * */
void jacobiV2(
    std::vector<double> &y,
    InitialConditions initialConditions,
    int numberOfProcesses,
    int processId)
{
    //Number of nodes in the process domain after division between processes
    int localSize;

    //Displacement in relation with beginning of the result vector
    int localDisplacement;

    //Norm on current iteration in process
    double localNorm;
    //Maximum norm between all processes
    double globalNorm;

    int iterationsNumber = 0;

    //Difference with analytic solution
    double difference = 0.0;

    //Part of the solution, here is situated part with which every process works
    std::vector<double> yLocal;
    std::vector<double> yLocalPrevious;

    //Time counters
    double timeStart, timeEnd;

    //Number of rows which belongs to process
    int localRows;

    //Offset in rows from the beginning
    int localOffsetInRows;

    //Just coefficient of the equation
    //double c = 1.0 / (4.0 + initialConditions.h * initialConditions.h * initialConditions.k * initialConditions.k);

    //This vectors are needed to handle data transmission
    //They store first and last rows of process part
    std::vector<double> yLocalPreviousUpHighBorder;
    std::vector<double> yLocalPreviousDownLowBorder;

    //Vectors to handle data transmission
    std::vector<double> buf1;
    std::vector<double> buf2;

    //@processesLocationSizes stores data about number of nodes which belongs to each process
    std::vector<int> processesLocationSizes;
    //@processesLocationSizes stores data about displacement from the beginning of y
    std::vector<int> processesDisplacement;

    divideResponsibilities(y,
                           yLocal,
                           yLocalPrevious,
                           numberOfProcesses,
                           processId,
                           localSize,
                           localDisplacement,
                           localRows,
                           localOffsetInRows,
                           processesLocationSizes,
                           processesDisplacement,
                           initialConditions);

    init(processId,
         y,
         yLocal,
         yLocalPrevious,
         yLocalPreviousUpHighBorder,
         yLocalPreviousDownLowBorder,
         buf1,
         buf2,
         initialConditions,
         localSize);

    printProcessData(yLocal,
                     yLocalPrevious,
                     processId,
                     localSize,
                     localDisplacement,
                     localRows,
                     localOffsetInRows,
                     initialConditions.isDebugMode);

    MPI_Barrier(MPI_COMM_WORLD);
    if (processId == 0)
        timeStart = MPI_Wtime();

    do
    {
        iterationsNumber++;

        yLocal.swap(yLocalPrevious);

        exchangeDataV2(
            yLocal,
            yLocalPrevious,
            yLocalPreviousUpHighBorder,
            yLocalPreviousDownLowBorder,
            buf1,
            buf2,
            numberOfProcesses,
            processId,
            initialConditions);

        solveSystem(
            yLocal,
            yLocalPrevious,
            yLocalPreviousUpHighBorder,
            yLocalPreviousDownLowBorder,
            numberOfProcesses,
            processId,
            localRows,
            localSize,
            localOffsetInRows,
            initialConditions);

        localNorm = infiniteNorm(yLocal, yLocalPrevious);

        //If maximum norm is decreasing that means that all norms are decreasing
        MPI_Allreduce(&localNorm, &globalNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    } while (globalNorm > initialConditions.epsilon);

    if (processId == 0)
        timeEnd = MPI_Wtime();

    //Filling answer
    MPI_Gatherv(yLocal.data(), yLocal.size(), MPI_DOUBLE, y.data(), processesLocationSizes.data(),
                processesDisplacement.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    if (processId == 0)
        printMethodStatistic(globalNorm,
                             iterationsNumber,
                             timeStart,
                             timeEnd,
                             difference,
                             true);
}

void exchangeDataV2(std::vector<double> &yLocal,
                    std::vector<double> &yLocalPrevious,
                    std::vector<double> &yLocalPreviousUpHighBorder,
                    std::vector<double> &yLocalPreviousDownLowBorder,
                    std::vector<double> &buf1,
                    std::vector<double> &buf2,
                    int numberOfProcesses,
                    int processId,
                    InitialConditions initialConditions)
{
    MPI_Status statU, statL;
    int lowerRankProcess;
    int higherRankProcess;

    setSourceAndDestination(numberOfProcesses, processId, higherRankProcess, lowerRankProcess);

    copyLastRow(yLocalPrevious, buf1, initialConditions);
    MPI_Sendrecv(buf1.data(), initialConditions.n, MPI_DOUBLE, higherRankProcess, 56,
                 yLocalPreviousUpHighBorder.data(), initialConditions.n, MPI_DOUBLE, lowerRankProcess, 56, MPI_COMM_WORLD, &statU);

    copyFirstRow(yLocalPrevious, buf2, initialConditions);
    MPI_Sendrecv(buf2.data(), initialConditions.n, MPI_DOUBLE, lowerRankProcess, 100,
                 yLocalPreviousDownLowBorder.data(), initialConditions.n, MPI_DOUBLE, higherRankProcess, 100, MPI_COMM_WORLD, &statL);
}

void setSourceAndDestination(const int numberOfProcesses,
                             const int processId,
                             int &higherRankProcess,
                             int &lowerRankProcess)
{
    higherRankProcess = processId + 1;

    lowerRankProcess = processId - 1;

    if (processId == 0)
        lowerRankProcess = numberOfProcesses - 1;
    else if (processId == numberOfProcesses - 1)
        higherRankProcess = 0;
}
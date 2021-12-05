#include "Jacobi.hpp"

/**
 * Used methods: MPI_Recv_init. MPI_Send_init
 * */
void jacobiV3(
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

    //Only processes not on the edges of domain can transfer data to top and to bottom
    std::vector<MPI_Request> req1, req2;
    if ((processId == 0) || (processId == numberOfProcesses - 1))
    {
        req1.resize(1);
        req2.resize(1);
    }
    else
    {
        req1.resize(2);
        req2.resize(2);
    }
    int requestsToTop = 0, requestsToBottom = 0;

    if (processId != numberOfProcesses - 1)
    {
        MPI_Send_init(buf1.data(), initialConditions.n, MPI_DOUBLE, processId + 1, 50, MPI_COMM_WORLD, req1.data() + requestsToTop);
        MPI_Recv_init(yLocalPreviousDownLowBorder.data(), initialConditions.n,
                      MPI_DOUBLE, processId + 1, 70, MPI_COMM_WORLD, req2.data() + requestsToBottom);
        requestsToTop++;
        requestsToBottom++;
    }

    if (processId != 0)
    {
        MPI_Recv_init(yLocalPreviousUpHighBorder.data(), initialConditions.n, MPI_DOUBLE, processId - 1,
                      50, MPI_COMM_WORLD, req1.data() + requestsToTop);
        MPI_Send_init(buf2.data(), initialConditions.n, MPI_DOUBLE, processId - 1,
                      70, MPI_COMM_WORLD, req2.data() + requestsToBottom);
        requestsToTop++;
        requestsToBottom++;
    }

    std::vector<MPI_Status> stat1(requestsToTop), stat2(requestsToBottom);

    MPI_Barrier(MPI_COMM_WORLD);
    if (processId == 0)
        timeStart = MPI_Wtime();

    do
    {
        iterationsNumber++;
        //std::cout << iterationsNumber << ::std::endl;

        yLocal.swap(yLocalPrevious);

        copyLastRow(yLocalPrevious, buf1, initialConditions);
        MPI_Startall(requestsToTop, req1.data());

        copyFirstRow(yLocalPrevious, buf2, initialConditions);
        MPI_Startall(requestsToBottom, req2.data());

        MPI_Waitall(requestsToBottom, req2.data(), stat2.data());
        MPI_Waitall(requestsToTop, req1.data(), stat1.data());

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
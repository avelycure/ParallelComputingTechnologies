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

    //Part of the solution, here is situated part with which every process works
    std::vector<double> yLocal;
    std::vector<double> yLocalPrevious;

    //Time counters
    double timeStart, timeEnd;

    //Number of rows which belongs to process
    int localRows;

    //Offset in rows from the beginning
    int localOffsetInRows;

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
    //@requestsFromLowerToHigher is for sending data from lower rank processes to higher
    //@requestsFromHigherToLower is for sending data from higher rank processes ro lower
    std::vector<MPI_Request> requestsFromLowerToHigher, requestsFromHigherToLower;

    //Variables for easy handling number of requests
    //@highRequests is for sending to the top(to lower rank processes) and to receive from the top
    //@lowRequests if for sending to the bottom(to higher rank processes) and to receive from bottom
    int highRequests = 0, lowRequests = 0;

    prepareJacobiRequests(processId,
                          numberOfProcesses,
                          initialConditions,
                          lowRequests,
                          highRequests,
                          yLocalPreviousUpHighBorder,
                          yLocalPreviousDownLowBorder,
                          buf1,
                          buf2,
                          requestsFromLowerToHigher,
                          requestsFromHigherToLower);

    std::vector<MPI_Status> transactionStateFromTop(highRequests), transactionStateFromBottom(lowRequests);

    MPI_Barrier(MPI_COMM_WORLD);
    if (processId == 0)
        timeStart = MPI_Wtime();

    do
    {
        iterationsNumber++;

        yLocal.swap(yLocalPrevious);

        solveSystemV3(yLocal,
                      yLocalPrevious,
                      yLocalPreviousUpHighBorder,
                      yLocalPreviousDownLowBorder,
                      numberOfProcesses,
                      processId,
                      localRows,
                      localSize,
                      localOffsetInRows,
                      lowRequests,
                      highRequests,
                      initialConditions,
                      buf1,
                      buf2,
                      requestsFromLowerToHigher,
                      requestsFromHigherToLower,
                      transactionStateFromTop,
                      transactionStateFromBottom);

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
        printMethodStatistic("Jacobi.MPI_Send_init. MPI_Recv_init",
                             globalNorm,
                             iterationsNumber,
                             timeStart,
                             timeEnd,
                             true);
}

/**
 * Make requests to higher rank process and lower rank process
 * */
void prepareJacobiRequests(int processId,
                           int numberOfProcesses,
                           InitialConditions initialConditions,
                           int &lowRequests,
                           int &highRequests,
                           std::vector<double> &yLocalPreviousUpHighBorder,
                           std::vector<double> &yLocalPreviousDownLowBorder,
                           std::vector<double> &buf1,
                           std::vector<double> &buf2,
                           std::vector<MPI_Request> &requestsFromLowerToHigher,
                           std::vector<MPI_Request> &requestsFromHigherToLower)
{
    //Only the first and the last process can send data to one direction others can send in twice
    if ((processId == 0) || (processId == numberOfProcesses - 1))
    {
        requestsFromLowerToHigher.resize(1);
        requestsFromHigherToLower.resize(1);
    }
    else
    {
        requestsFromLowerToHigher.resize(2);
        requestsFromHigherToLower.resize(2);
    }

    //All processes except numberOfProcesses - 1 rank, could send data forward and can receive data from higher rank process
    if (processId != numberOfProcesses - 1)
    {
        //sending last row to next process
        MPI_Send_init(buf1.data(), initialConditions.n, MPI_DOUBLE, processId + 1, 50, MPI_COMM_WORLD, requestsFromLowerToHigher.data());
        MPI_Recv_init(yLocalPreviousDownLowBorder.data(), initialConditions.n,
                      MPI_DOUBLE, processId + 1, 70, MPI_COMM_WORLD, requestsFromHigherToLower.data());
        highRequests++;
        lowRequests++;
    }

    //All processes except 0 rank can receive data from lower rank processes and can send data to them
    if (processId != 0)
    {
        MPI_Recv_init(yLocalPreviousUpHighBorder.data(), initialConditions.n, MPI_DOUBLE, processId - 1,
                      50, MPI_COMM_WORLD, requestsFromLowerToHigher.data() + highRequests);
        MPI_Send_init(buf2.data(), initialConditions.n, MPI_DOUBLE, processId - 1,
                      70, MPI_COMM_WORLD, requestsFromHigherToLower.data() + lowRequests);
        highRequests++;
        lowRequests++;
    }
}

/**
 * Solving system using MPI functions MPI_Send_init and MPI_Recv_init. This method is similar to solveSystem() but we can here
 * we can request data only when we need it. So we wait for first row in the begining of the function and we request last row at the end
 * */
void solveSystemV3(std::vector<double> &yLocal,
                   std::vector<double> &yLocalPrevious,
                   std::vector<double> &yLocalPreviousUpHighBorder,
                   std::vector<double> &yLocalPreviousDownLowBorder,
                   int numberOfProcesses,
                   int processId,
                   int localRows,
                   int localSize,
                   int localOffsetInRows,
                   int &lowRequests,
                   int &highRequests,
                   InitialConditions initialConditions,
                   std::vector<double> &buf1,
                   std::vector<double> &buf2,
                   std::vector<MPI_Request> &requestsFromLowerToHigher,
                   std::vector<MPI_Request> &requestsFromHigherToLower,
                   std::vector<MPI_Status> &transactionStateFromTop,
                   std::vector<MPI_Status> &transactionStateFromBottom)
{
    //Local renaming to increase readability
    double h = initialConditions.h;
    int n = initialConditions.n;
    //Just coefficient of the equation
    double c = 1.0 / (4.0 + initialConditions.h * initialConditions.h * initialConditions.k * initialConditions.k);

    //Begin transfering data, as we need the first row to begin solving system
    copyLastRow(yLocalPrevious, buf1, initialConditions);
    MPI_Startall(highRequests, requestsFromLowerToHigher.data());

    //We need last row so we begin transfering data in first rows to lower rank processes
    copyFirstRow(yLocalPrevious, buf2, initialConditions);
    MPI_Startall(lowRequests, requestsFromHigherToLower.data());

    //Calculate only rows that are not borders of process part
    for (int i = 1; i < localRows - 1; i++)
        for (int j = 1; j < n - 1; j++)
            yLocal[i * n + j] = c * (h * h * initialConditions.f((localOffsetInRows + i) * h, j * h) +
                                     yLocalPrevious[(i - 1) * n + j] +
                                     yLocalPrevious[(i + 1) * n + j] +
                                     yLocalPrevious[i * n + (j - 1)] +
                                     yLocalPrevious[i * n + (j + 1)]);

    //Wait until we get the row before first to begin calculation
    MPI_Waitall(highRequests, requestsFromLowerToHigher.data(), transactionStateFromTop.data());
    if (processId != 0)
        for (int j = 1; j < n - 1; j++)
            yLocal[j] = c * (h * h * initialConditions.f(localOffsetInRows * h, j * h) +
                             yLocalPreviousUpHighBorder[j] +
                             yLocalPrevious[n + j] +
                             yLocalPrevious[j - 1] +
                             yLocalPrevious[j + 1]);
    //Else if it is first process do nothing because initial conditions in the first row are zeros

    //Wait until we get the row after last to begin calculation
    MPI_Waitall(lowRequests, requestsFromHigherToLower.data(), transactionStateFromBottom.data());
    if (processId != numberOfProcesses - 1)
        for (int j = 1; j < n - 1; j++)
            yLocal[localSize - n + j] = c * (h * h * initialConditions.f((localOffsetInRows + localRows - 1) * h, j * h) +
                                             yLocalPrevious[localSize - 2 * n + j] +
                                             yLocalPreviousDownLowBorder[j] +
                                             yLocalPrevious[localSize - n + j - 1] +
                                             yLocalPrevious[localSize - n + j + 1]);
}
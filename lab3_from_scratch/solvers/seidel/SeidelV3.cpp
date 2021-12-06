#include "Seidel.hpp"

/**
 * Used methods: MPI_Send_init. MPI_Recv_init
 * */
void seidelV3(
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

    //Vectors to handle black and red requests from higher rank processes to lower and vice versa
    std::vector<MPI_Request> requestsFromLowerToHigherBlack, requestsFromHigherToLowerBlack;
    std::vector<MPI_Request> requestsFromLowerToHigherRed, requestsFromHigherToLowerRed;
    int highRequestsBlack = 0, lowRequestsBlack = 0, highRequestsRed = 0, lowRequestsRed = 0;

    prepareRequests(processId,
                    numberOfProcesses,
                    initialConditions,
                    lowRequestsBlack,
                    highRequestsBlack,
                    lowRequestsRed,
                    highRequestsRed,
                    yLocalPreviousUpHighBorder,
                    yLocalPreviousDownLowBorder,
                    buf1,
                    buf2,
                    requestsFromLowerToHigherBlack,
                    requestsFromHigherToLowerBlack,
                    requestsFromLowerToHigherRed,
                    requestsFromHigherToLowerRed);

    std::vector<MPI_Status> stateFromLowerToHigherBlack(highRequestsBlack), stateFromHigherToLowerBlack(lowRequestsBlack);
    std::vector<MPI_Status> stateFromLowerToHigherRed(highRequestsRed), stateFromHigherToLowerRed(lowRequestsRed);

    MPI_Barrier(MPI_COMM_WORLD);
    if (processId == 0)
        timeStart = MPI_Wtime();

    do
    {
        iterationsNumber++;

        yLocal.swap(yLocalPrevious);

        solveBlackV3(processId,
                     numberOfProcesses,
                     initialConditions,
                     lowRequestsBlack,
                     highRequestsBlack,
                     localRows,
                     localSize,
                     localOffsetInRows,
                     buf1,
                     buf2,
                     yLocal,
                     yLocalPrevious,
                     yLocalPreviousUpHighBorder,
                     yLocalPreviousDownLowBorder,
                     stateFromLowerToHigherBlack,
                     stateFromHigherToLowerBlack,
                     requestsFromLowerToHigherBlack,
                     requestsFromHigherToLowerBlack);

        solveRedV3(processId,
                   numberOfProcesses,
                   initialConditions,
                   lowRequestsRed,
                   highRequestsRed,
                   localRows,
                   localSize,
                   localOffsetInRows,
                   buf1,
                   buf2,
                   yLocal,
                   yLocalPrevious,
                   yLocalPreviousUpHighBorder,
                   yLocalPreviousDownLowBorder,
                   stateFromLowerToHigherRed,
                   stateFromHigherToLowerRed,
                   requestsFromLowerToHigherRed,
                   requestsFromHigherToLowerRed);

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
        printMethodStatistic("Seidel. MPI_Send_init. MPI_Recv_init",
                             globalNorm,
                             iterationsNumber,
                             timeStart,
                             timeEnd,
                             true);
}

/**
 * Make requests to higher rank process and lower rank process
 * */
void prepareRequests(int processId,
                     int numberOfProcesses,
                     InitialConditions initialConditions,
                     int &lowRequestsBlack,
                     int &highRequestsBlack,
                     int &lowRequestsRed,
                     int &highRequestsRed,
                     std::vector<double> &yLocalPreviousUpHighBorder,
                     std::vector<double> &yLocalPreviousDownLowBorder,
                     std::vector<double> &buf1,
                     std::vector<double> &buf2,
                     std::vector<MPI_Request> &requestsFromLowerToHigherBlack,
                     std::vector<MPI_Request> &requestsFromHigherToLowerBlack,
                     std::vector<MPI_Request> &requestsFromLowerToHigherRed,
                     std::vector<MPI_Request> &requestsFromHigherToLowerRed)
{
    //Only the first and the last process can send data to one direction others can send in twice
    if ((processId == 0) || (processId == numberOfProcesses - 1))
    {
        requestsFromLowerToHigherBlack.resize(1);
        requestsFromHigherToLowerBlack.resize(1);
        requestsFromLowerToHigherRed.resize(1);
        requestsFromHigherToLowerRed.resize(1);
    }
    else
    {
        requestsFromLowerToHigherBlack.resize(2);
        requestsFromHigherToLowerBlack.resize(2);
        requestsFromLowerToHigherRed.resize(2);
        requestsFromHigherToLowerRed.resize(2);
    }

    // black
    if (processId != numberOfProcesses - 1)
    {
        MPI_Send_init(buf1.data(), initialConditions.n, MPI_DOUBLE, processId + 1, 56, MPI_COMM_WORLD, requestsFromLowerToHigherBlack.data());
        MPI_Recv_init(yLocalPreviousDownLowBorder.data(), initialConditions.n, MPI_DOUBLE, processId + 1, 65, MPI_COMM_WORLD, requestsFromHigherToLowerBlack.data());
        lowRequestsBlack++;
        highRequestsBlack++;
    }

    if (processId != 0)
    {
        MPI_Recv_init(yLocalPreviousUpHighBorder.data(), initialConditions.n, MPI_DOUBLE, processId - 1, 56, MPI_COMM_WORLD, requestsFromLowerToHigherBlack.data() + highRequestsBlack);
        MPI_Send_init(buf2.data(), initialConditions.n, MPI_DOUBLE, processId - 1, 65, MPI_COMM_WORLD, requestsFromHigherToLowerBlack.data() + lowRequestsBlack);
        lowRequestsBlack++;
        highRequestsBlack++;
    }

    //red
    if (processId != numberOfProcesses - 1)
    {
        MPI_Send_init(buf1.data(), initialConditions.n, MPI_DOUBLE, processId + 1, 56, MPI_COMM_WORLD, requestsFromLowerToHigherRed.data());
        MPI_Recv_init(yLocalPreviousDownLowBorder.data(), initialConditions.n, MPI_DOUBLE, processId + 1, 65, MPI_COMM_WORLD, requestsFromHigherToLowerRed.data());
        lowRequestsRed++;
        highRequestsRed++;
    }
    if (processId != 0)
    {
        MPI_Recv_init(yLocalPreviousUpHighBorder.data(), initialConditions.n, MPI_DOUBLE, processId - 1, 56, MPI_COMM_WORLD, requestsFromLowerToHigherRed.data() + highRequestsRed);
        MPI_Send_init(buf2.data(), initialConditions.n, MPI_DOUBLE, processId - 1, 65, MPI_COMM_WORLD, requestsFromHigherToLowerRed.data() + lowRequestsRed);
        lowRequestsRed++;
        highRequestsRed++;
    }
}

/**
 * Solving system using MPI functions MPI_Send_init and MPI_Recv_init. This method is similar to solveBlack() but we can here
 * we can request data only when we need it. So we wait for first row in the begining of the function and we request last row at the end
 * */
void solveBlackV3(int processId,
                  int numberOfProcesses,
                  InitialConditions initialConditions,
                  int &lowRequestsBlack,
                  int &highRequestsBlack,
                  int &localRows,
                  int &localSize,
                  int &localOffsetInRows,
                  std::vector<double> &buf1,
                  std::vector<double> &buf2,
                  std::vector<double> &yLocal,
                  std::vector<double> &yLocalPrevious,
                  std::vector<double> &yLocalPreviousUpHighBorder,
                  std::vector<double> &yLocalPreviousDownLowBorder,
                  std::vector<MPI_Status> &stateFromLowerToHigherBlack,
                  std::vector<MPI_Status> &stateFromHigherToLowerBlack,
                  std::vector<MPI_Request> &requestsFromLowerToHigherBlack,
                  std::vector<MPI_Request> &requestsFromHigherToLowerBlack)
{
    //Local renaming to increase readability
    double h = initialConditions.h;
    int n = initialConditions.n;
    //Just coefficient of the equation
    double c = 1.0 / (4.0 + initialConditions.h * initialConditions.h * initialConditions.k * initialConditions.k);

    //Exchanging previous! values of solution
    copyLastRow(yLocalPrevious, buf1, initialConditions);
    MPI_Startall(highRequestsBlack, requestsFromLowerToHigherBlack.data());

    MPI_Waitall(highRequestsBlack, requestsFromLowerToHigherBlack.data(), stateFromLowerToHigherBlack.data());
    //In this part we solve equation with previous! values of yLocal

    if (processId != 0)
        for (int j = (localOffsetInRows % 2 == 0) ? 2 : 1; j < n - 1; j += 2)
            yLocal[j] = c * (h * h * initialConditions.f(localOffsetInRows * h, j * h) +
                             yLocalPreviousUpHighBorder[j] +
                             yLocalPrevious[n + j] +
                             yLocalPrevious[j - 1] +
                             yLocalPrevious[j + 1]);
    //Else if it is first process do nothing because initial conditions in the first row are zeros

    //Calculate only rows that are not borders of process part
    for (int i = 1; i < localRows - 1; i++)
        for (int j = ((localOffsetInRows + i) % 2 == 0) ? 2 : 1; j < n - 1; j += 2)
            yLocal[i * n + j] = c * (h * h * initialConditions.f((localOffsetInRows + i) * h, j * h) +
                                     yLocalPrevious[(i - 1) * n + j] +
                                     yLocalPrevious[(i + 1) * n + j] +
                                     yLocalPrevious[i * n + (j - 1)] +
                                     yLocalPrevious[i * n + (j + 1)]);

    //Exchanging current! values of solution
    copyFirstRow(yLocalPrevious, buf2, initialConditions);
    MPI_Startall(lowRequestsBlack, requestsFromHigherToLowerBlack.data());
    MPI_Waitall(lowRequestsBlack, requestsFromHigherToLowerBlack.data(), stateFromHigherToLowerBlack.data());

    if (processId != numberOfProcesses - 1)
        for (int j = ((localOffsetInRows + localRows - 1) % 2 == 0) ? 2 : 1; j < n - 1; j += 2)
            yLocal[localSize - n + j] = c * (h * h * initialConditions.f((localOffsetInRows + localRows - 1) * h, j * h) +
                                             yLocalPrevious[localSize - 2 * n + j] +
                                             yLocalPreviousDownLowBorder[j] +
                                             yLocalPrevious[localSize - n + j - 1] +
                                             yLocalPrevious[localSize - n + j + 1]);
}

/**
 * Solving system using MPI functions MPI_Send_init and MPI_Recv_init. This method is similar to solveBlack() but we can here
 * we can request data only when we need it. So we wait for first row in the begining of the function and we request last row at the end
 * */
void solveRedV3(int processId,
                int numberOfProcesses,
                InitialConditions initialConditions,
                int &lowRequestsRed,
                int &highRequestsRed,
                int &localRows,
                int &localSize,
                int &localOffsetInRows,
                std::vector<double> &buf1,
                std::vector<double> &buf2,
                std::vector<double> &yLocal,
                std::vector<double> &yLocalPrevious,
                std::vector<double> &yLocalPreviousUpHighBorder,
                std::vector<double> &yLocalPreviousDownLowBorder,
                std::vector<MPI_Status> &stateFromLowerToHigherRed,
                std::vector<MPI_Status> &stateFromHigherToLowerRed,
                std::vector<MPI_Request> &requestsFromLowerToHigherRed,
                std::vector<MPI_Request> &requestsFromHigherToLowerRed)
{
    //Local renaming to increase readability
    double h = initialConditions.h;
    int n = initialConditions.n;
    //Just coefficient of the equation
    double c = 1.0 / (4.0 + initialConditions.h * initialConditions.h * initialConditions.k * initialConditions.k);

    //Exchanging previous! values of solution
    copyLastRow(yLocal, buf1, initialConditions);
    MPI_Startall(highRequestsRed, requestsFromLowerToHigherRed.data());
    MPI_Waitall(highRequestsRed, requestsFromLowerToHigherRed.data(), stateFromLowerToHigherRed.data());

    //In this part we solve equation with current! values of yLocal
    if (processId != 0)
        for (int j = (localOffsetInRows % 2 == 0) ? 1 : 2; j < n - 1; j += 2)
            yLocal[j] = c * (h * h * initialConditions.f(localOffsetInRows * h, j * h) +
                             yLocalPreviousUpHighBorder[j] +
                             yLocal[n + j] +
                             yLocal[j - 1] +
                             yLocal[j + 1]);
    //Else if it is first process do nothing because initial conditions in the first row are zeros

    //Calculate only rows that are not borders of process part
    for (int i = 1; i < localRows - 1; i++)
        for (int j = ((localOffsetInRows + i) % 2 == 0) ? 1 : 2; j < n - 1; j += 2)
            yLocal[i * n + j] = c * (h * h * initialConditions.f((localOffsetInRows + i) * h, j * h) +
                                     yLocal[(i - 1) * n + j] +
                                     yLocal[(i + 1) * n + j] +
                                     yLocal[i * n + (j - 1)] +
                                     yLocal[i * n + (j + 1)]);

    //Exchanging current! values of solution
    copyFirstRow(yLocal, buf2, initialConditions);
    MPI_Startall(lowRequestsRed, requestsFromHigherToLowerRed.data());
    MPI_Waitall(lowRequestsRed, requestsFromHigherToLowerRed.data(), stateFromHigherToLowerRed.data());

    if (processId != numberOfProcesses - 1)
        for (int j = ((localOffsetInRows + localRows - 1) % 2 == 0) ? 1 : 2; j < n - 1; j += 2)
            yLocal[localSize - n + j] = c * (h * h * initialConditions.f((localOffsetInRows + localRows - 1) * h, j * h) +
                                             yLocal[localSize - 2 * n + j] +
                                             yLocalPreviousDownLowBorder[j] +
                                             yLocal[localSize - n + j - 1] +
                                             yLocal[localSize - n + j + 1]);
}
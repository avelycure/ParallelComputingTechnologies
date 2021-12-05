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

    //This vectors are needed to handle data transmission
    //They store first and last rows of process part
    std::vector<double> yLocalPreviousUpHighBorder;
    std::vector<double> yLocalPreviousDownLowBorder;

    std::vector<double> buf1;
    buf1.resize(initialConditions.n);
    std::vector<double> buf2;
    buf2.resize(initialConditions.n);

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

    //Local renaming to increase readability
    double h = initialConditions.h;
    int n = initialConditions.n;

    /*for (int i = 0; i < yLocal.size(); i++)
    {
        yLocal[i] = i + processId * 100;
        //yLocalPreviousUpHighBorder[i] = i + processId * 10;
        //yLocalPreviousDownLowBorder[i] = i + processId * 100 + 50;
    }*/

    do
    {
        iterationsNumber++;

        yLocal.swap(yLocalPrevious);

        //Send data from lower rank processes to higher
        if (processId != numberOfProcesses - 1)
        {
            //copyLastRow(yLocalPrevious, yLocalPreviousDownLowBorder, initialConditions);
            //MPI_Send(yLocalPreviousDownLowBorder.data(), n, MPI_DOUBLE, processId + 1, 42, MPI_COMM_WORLD);
            copyLastRow(yLocalPrevious, buf1, initialConditions);
            MPI_Send(buf1.data(), n, MPI_DOUBLE, processId + 1, 42, MPI_COMM_WORLD);
        }
        if (processId != 0)
            MPI_Recv(yLocalPreviousUpHighBorder.data(), n, MPI_DOUBLE, processId - 1, 42, MPI_COMM_WORLD, &statL);

        //Send data from higher rank processes to lower
        if (processId != 0)
        {
            //copyFirstRow(yLocalPrevious, yLocalPreviousUpHighBorder, initialConditions);
            //MPI_Send(yLocalPreviousUpHighBorder.data(), n, MPI_DOUBLE, processId - 1, 41, MPI_COMM_WORLD);
            copyFirstRow(yLocalPrevious, buf2, initialConditions);
            MPI_Send(buf2.data(), n, MPI_DOUBLE, processId - 1, 41, MPI_COMM_WORLD);
        }
        if (processId != numberOfProcesses - 1)
            MPI_Recv(yLocalPreviousDownLowBorder.data(), n, MPI_DOUBLE, processId + 1, 41, MPI_COMM_WORLD, &statU);

        /*for (int i = 0; i < yLocalPreviousUpHighBorder.size(); i++)
    {
        std::cout << "UProcess" << processId << " [" << i << "]=" << yLocalPreviousUpHighBorder[i] << std::endl;
    }

    for (int i = 0; i < yLocalPreviousUpHighBorder.size(); i++)
    {
        std::cout << "DProcess" << processId << " [" << i << "]=" << yLocalPreviousDownLowBorder[i] << std::endl;
    }*/

        if (processId != 0)
            for (int j = 1; j < n - 1; j++)
                yLocal[j] = c * (h * h * initialConditions.f(localOffsetInRows * h, j * h) +
                                 yLocalPreviousUpHighBorder[j] +
                                 yLocalPrevious[n + j] +
                                 yLocalPrevious[j - 1] +
                                 yLocalPrevious[j + 1]);
        //Else if it is first process do nothing because initial conditions in the first row are zeros

        //Calculate only rows that are not borders of process part
        for (int i = 1; i < localRows - 1; i++)
            for (int j = 1; j < n - 1; j++)
                yLocal[i * n + j] = c * (h * h * initialConditions.f((localOffsetInRows + i) * h, j * h) +
                                         yLocalPrevious[(i - 1) * n + j] +
                                         yLocalPrevious[(i + 1) * n + j] +
                                         yLocalPrevious[i * n + (j - 1)] +
                                         yLocalPrevious[i * n + (j + 1)]);

        if (processId != numberOfProcesses - 1)
            for (int j = 1; j < n - 1; j++)
                yLocal[localSize - n + j] = c * (h * h * initialConditions.f((localOffsetInRows + localRows - 1) * h, j * h) +
                                                 yLocalPrevious[localSize - 2 * n + j] +
                                                 yLocalPreviousDownLowBorder[j] +
                                                 yLocalPrevious[localSize - n + j - 1] +
                                                 yLocalPrevious[localSize - n + j + 1]);
        //Else if it is last process do nothing because initial conditions in the last row are zeros

        localNorm = infiniteNorm(yLocal, yLocalPrevious);

        MPI_Allreduce(&localNorm, &globalNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        //std::cout << iterationsNumber << "Process id: " << processId << " local norm: " << localNorm << " global: " << globalNorm << std::endl;

    } while (globalNorm > initialConditions.epsilon);

    if (processId == 0)
        timeEnd = MPI_Wtime();

    MPI_Gatherv(yLocal.data(), yLocal.size(), MPI_DOUBLE, y.data(), processesLocationSizes.data(),
                processesDisplacement.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    if (processId == 0)
        printMethodStatistic(globalNorm,
                             iterationsNumber,
                             timeStart,
                             timeEnd,
                             difference,
                             initialConditions.isDebugMode);
}

/**
 * Copy first row of local part of solution
 * This function is needed for transmition
 * */
void copyFirstRow(std::vector<double> &yLocal,
                  std::vector<double> &localHighBorder,
                  InitialConditions initialConditions)
{
    for (int i = 0; i < initialConditions.n; i++)
        localHighBorder[i] = yLocal[i];
}

/**
 * Copy last row of local part of solution
 * This function is needed for transmition
 * */
void copyLastRow(std::vector<double> &yLocal,
                 std::vector<double> &localLowBorder,
                 InitialConditions initialConditions)
{
    for (int i = 0; i < initialConditions.n; i++)
        localLowBorder[i] = yLocal[yLocal.size() - initialConditions.n + i];
}
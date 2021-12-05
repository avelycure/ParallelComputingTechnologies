#include "Common.hpp"

/**
 * We have to exchange two rows: the first and the last. We have two buffers to simplify readability: @buf1, @buf2.
 * Firstly we send from top of the matrix to the bottom. The schema is like this: 0->1->2->3->... where numbers are processes ranks
 * Then we send data from bottom to the top: 4->3->2->1->0. Each process(except 0 and numberOfProcesses - 1,  which are
 * special cases) has one row before(or up) it and one row after. So we have two vectors to work with: @yLocalPreviousUpHighBorder
 * and @yLocalPreviousDownLowBorder. @yLocal is vector of solution on current iteration
 * */
void exchangeDataV1(std::vector<double> &yLocal,
                    std::vector<double> &yLocalPreviousUpHighBorder,
                    std::vector<double> &yLocalPreviousDownLowBorder,
                    std::vector<double> &buf1,
                    std::vector<double> &buf2,
                    int numberOfProcesses,
                    int processId,
                    InitialConditions initialConditions)
{
    //Variables to get transaction status
    MPI_Status statU, statL;

    //Send data from lower rank processes to higher
    if (processId != numberOfProcesses - 1)
    {
        copyLastRow(yLocal, buf1, initialConditions);
        MPI_Send(buf1.data(), initialConditions.n, MPI_DOUBLE, processId + 1, 800, MPI_COMM_WORLD);
    }
    if (processId != 0)
        MPI_Recv(yLocalPreviousUpHighBorder.data(), initialConditions.n, MPI_DOUBLE, processId - 1, 800, MPI_COMM_WORLD, &statL);

    //Send data from higher rank processes to lower
    if (processId != 0)
    {
        copyFirstRow(yLocal, buf2, initialConditions);
        MPI_Send(buf2.data(), initialConditions.n, MPI_DOUBLE, processId - 1, 25, MPI_COMM_WORLD);
    }
    if (processId != numberOfProcesses - 1)
        MPI_Recv(yLocalPreviousDownLowBorder.data(), initialConditions.n, MPI_DOUBLE, processId + 1, 25, MPI_COMM_WORLD, &statU);
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
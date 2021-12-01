#include "Main.hpp"

int main(int argc, char **argv)
{
    int processId;
    int numberOfProcesses;

    //MPI_MAX PROCESSOR_NAME is a constant: length(any processorName) < it
    char processorName[MPI_MAX_PROCESSOR_NAME];
    int processorNameLength;

    vector<double> y;
    vector<double> solution;

    MPI_Init(&argc, &argv);

    //get number of processes in world communicator
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

    //get current process num
    MPI_Comm_rank(MPI_COMM_WORLD, &processId);

    //get process name and lenght of process name
    MPI_Get_processor_name(processorName, &processorNameLength);

    InitialConditions initialCondition = InitialConditions();

    if (processId == 0)
    {

        printLog(
            numberOfProcesses,
            initialCondition.n,
            initialCondition.kSquare / (initialCondition.h * initialCondition.h),
            initialCondition.eps);

        y.resize((initialCondition.n + 1) * (initialCondition.n + 1));
        solution.resize((initialCondition.n + 1) * (initialCondition.n + 1));

        initialCondition.computeSolution(solution);
    }

    //Jacobi send receceive part

    //synchronization, all processes will be blocked until end of the procedure
    double time;
    MPI_Barrier(MPI_COMM_WORLD);
    jacobiSendReceive(y,
                      initialCondition.h,
                      initialCondition.n + 1,
                      initialCondition.kSquare,
                      numberOfProcesses,
                      processId,
                      initialCondition.eps,
                      time,
                      initialCondition);

    MPI_Finalize();
}

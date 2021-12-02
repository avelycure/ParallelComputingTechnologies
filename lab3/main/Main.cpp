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

    //Jacobi part
    // Synchronization, all processes will be blocked until end of the procedure
    MPI_Barrier(MPI_COMM_WORLD);
    //Used methods: MPI_Send. MPI_Recv
    jacobiV1(y,
             initialCondition.h,
             initialCondition.n + 1,
             initialCondition.kSquare,
             numberOfProcesses,
             processId,
             initialCondition.eps,
             initialCondition,
             solution);

    MPI_Barrier(MPI_COMM_WORLD);
    //Used methods: MPI_Sendrecv
    jacobiV2(y,
             initialCondition.h,
             initialCondition.n + 1,
             initialCondition.kSquare,
             numberOfProcesses,
             processId,
             initialCondition.eps,
             initialCondition,
             solution);

    MPI_Barrier(MPI_COMM_WORLD);
    //Used methods: MPI_Send_init, MPI_Recv_init
    jacobiV3(y,
             initialCondition.h,
             initialCondition.n + 1,
             initialCondition.kSquare,
             numberOfProcesses,
             processId,
             initialCondition.eps,
             initialCondition,
             solution);

    MPI_Barrier(MPI_COMM_WORLD);
    //Used methods: MPI_Send. MPI_Recv
    seidelV1(y,
             initialCondition.h,
             initialCondition.n + 1,
             initialCondition.kSquare,
             numberOfProcesses,
             processId,
             initialCondition.eps,
             initialCondition,
             solution);

    MPI_Barrier(MPI_COMM_WORLD);
    //Used methods: MPI_Sendrecv
    seidelV2(y,
             initialCondition.h,
             initialCondition.n + 1,
             initialCondition.kSquare,
             numberOfProcesses,
             processId,
             initialCondition.eps,
             initialCondition,
             solution);

    MPI_Barrier(MPI_COMM_WORLD);
    //Used methods: MPI_Send_init, MPI_Recv_init
    seidelV3(y,
             initialCondition.h,
             initialCondition.n + 1,
             initialCondition.kSquare,
             numberOfProcesses,
             processId,
             initialCondition.eps,
             initialCondition,
             solution);

    MPI_Finalize();
}

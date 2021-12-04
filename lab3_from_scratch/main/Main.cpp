#include "Main.hpp"

int main(int argc, char **argv)
{
    int processId;
    int numberOfProcesses;

    //MPI_MAX PROCESSOR_NAME is a constant: length(any processorName) < it
    char processorName[MPI_MAX_PROCESSOR_NAME];
    int processorNameLength;

    vector<double> y;
    //Analytic solution of the equation
    vector<double> solution;

    MPI_Init(&argc, &argv);

    //get number of processes in world communicator
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

    //get current process num
    MPI_Comm_rank(MPI_COMM_WORLD, &processId);

    //get process name and lenght of process name
    MPI_Get_processor_name(processorName, &processorNameLength);

    InitialConditions initialConditions = InitialConditions();

    if (processId == 0)
    {
        printLog(
            numberOfProcesses,
            initialConditions.n,
            initialConditions.epsilon);

        y.resize((initialConditions.n) * (initialConditions.n));
        solution.resize((initialConditions.n) * (initialConditions.n));

        initialConditions.computeSolution(solution);
    }

    jacobiV1(y,
             initialConditions,
             numberOfProcesses,
             processId);

    MPI_Finalize();
}

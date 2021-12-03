#include "Jacobi.hpp"

/**
 * Jacobi method of solving system. Result of this function is in y vector.
 * Used methods: MPI_Send. MPI_Recv
 * */
void jacobiV1(
    std::vector<double> &y,
    double h,
    int size,
    double kSquare,
    int numberOfProcesses,
    int processId,
    double eps,
    InitialConditions initialConditions,
    std::vector<double> &solution)
{
    int locationSize;
    int receiveDisplacement;
    int extraSize;
    int offset;

    int iterationsNumber = 0;
    double norm;
    double finalNorm = 0.0;

    double c = 1.0 / (4.0 + kSquare);
    kSquare = kSquare / (h * h);

    std::vector<double> rightPart;
    std::vector<int> numbersOfProcessDataParts;
    std::vector<int> displacement;
    std::vector<double> yPart;
    std::vector<double> yPreviousPart;
    std::vector<double> partOfRightPart;

    divideVectorBetweenProcesses(y, h, size, kSquare, numberOfProcesses, processId,
                                 yPart, yPreviousPart, partOfRightPart, numbersOfProcessDataParts, displacement,
                                 locationSize, receiveDisplacement, extraSize, offset);

    MPI_Status statU, statL;
    double timeStart;
    double timeEnd;

    if (processId == 0)
        timeStart = MPI_Wtime();

    do
    {
        iterationsNumber++;

        yPreviousPart.swap(yPart);

        //send data 0 -> 1, receive in 1
        for (int id = 0; id < numberOfProcesses - 1; id++)
        {
            if (processId == id)
                MPI_Send(yPreviousPart.data() + locationSize - 2 * size, size, MPI_DOUBLE, processId + 1, 42, MPI_COMM_WORLD);

            if (processId == id + 1)
                MPI_Recv(yPreviousPart.data(), size, MPI_DOUBLE, processId - 1, 42, MPI_COMM_WORLD, &statL);
        }

        // receive data from 1 to 0
        for (int id = numberOfProcesses - 1; id > 0; id--)
        {
            if (processId == id)
                MPI_Send(yPreviousPart.data() + size, size, MPI_DOUBLE, processId - 1, 41, MPI_COMM_WORLD);

            if (processId == id - 1)
                MPI_Recv(yPreviousPart.data() + locationSize - size, size, MPI_DOUBLE, processId + 1, 41, MPI_COMM_WORLD, &statU);
        }

        for (int i = 1; i < locationSize / size - 1; ++i)
            for (int j = 1; j < size - 1; ++j)
                yPart[i * size + j] = c * (h * h * initialConditions.f((i + offset) * h, j * h, kSquare) +
                                           yPreviousPart[(i - 1) * size + j] +
                                           yPreviousPart[(i + 1) * size + j] +
                                           yPreviousPart[i * size + (j - 1)] +
                                           yPreviousPart[i * size + (j + 1)]);

        norm = infiniteNorm(yPart, yPreviousPart, receiveDisplacement, locationSize - size);

        // one time find max from norm in all processes
        MPI_Allreduce(&norm, &finalNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    } while (finalNorm > eps);

    if (processId == 0)
        timeEnd = MPI_Wtime();

    //page 40, gathers different number of data from array
    MPI_Gatherv(yPart.data() + receiveDisplacement, locationSize - size - extraSize, MPI_DOUBLE, y.data(),
                numbersOfProcessDataParts.data(), displacement.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (processId == 0)
    {
        std::cout << "\033[1;32mJacobi. MPI_Send. MPI_Recv\033[0m" << std::endl;
        std::cout << "Substraction norm: " << finalNorm << std::endl;
        std::cout << "Number of iterations: " << iterationsNumber << std::endl;
        std::cout << "Time: " << timeEnd - timeStart << std::endl;
        double differenceWithAccurateSolution = infiniteNorm(y, solution, 0, y.size());
        std::cout << "Diffrence: " << differenceWithAccurateSolution << std::endl;
    }
};

/**
 * Used: MPI_Sendrecv
 * */
void jacobiV2(std::vector<double> &y,
              double h,
              int size,
              double kSquare,
              int numberOfProcesses,
              int processId,
              double eps,
              InitialConditions initialConditions,
              std::vector<double> &solution)
{
    int locationSize;
    int receiveDisplacement;
    int extraSize;
    int offset;

    std::vector<double> rightPart;
    std::vector<int> numbersOfProcessDataParts;
    std::vector<int> displacement;
    std::vector<double> yPart;
    std::vector<double> yPreviousPart;
    std::vector<double> partOfRightPart;

    int iterationsNumber = 0;
    double norm;
    double finalNorm;

    double c = 1.0 / (4.0 + kSquare);
    kSquare /= (h * h);

    int destination = 0;
    int source = 0;

    int sendCount = (processId - (numberOfProcesses - 1)) ? size : 0;
    int receiveCount = processId ? size : 0;

    double timeStart;
    double timeEnd;

    divideVectorBetweenProcesses(y, h, size, kSquare, numberOfProcesses, processId, yPart, yPreviousPart, partOfRightPart, numbersOfProcessDataParts, displacement,
                                 locationSize, receiveDisplacement, extraSize, offset);

    if (numberOfProcesses > 1)
        setSourceAndDestination(numberOfProcesses, processId, destination, source);

    if (processId == 0)
        timeStart = MPI_Wtime();

    do
    {
        iterationsNumber++;

        yPreviousPart.swap(yPart);

        MPI_Status statU, statL;

        //Block sending, page 33
        MPI_Sendrecv(yPreviousPart.data() + locationSize - 2 * size, sendCount, MPI_DOUBLE, destination, 56,
                     yPreviousPart.data(), receiveCount, MPI_DOUBLE, source, 56, MPI_COMM_WORLD, &statU);

        MPI_Sendrecv(yPreviousPart.data() + size, receiveCount, MPI_DOUBLE, source, 65,
                     yPreviousPart.data() + locationSize - size, sendCount, MPI_DOUBLE, destination, 65, MPI_COMM_WORLD, &statL);

        for (int i = 1; i < locationSize / size - 1; ++i)
            for (int j = 1; j < size - 1; ++j)
                yPart[i * size + j] = c * (h * h * initialConditions.f((i + offset) * h, j * h, kSquare) + yPreviousPart[(i - 1) * size + j] +
                                           yPreviousPart[(i + 1) * size + j] +
                                           yPreviousPart[i * size + (j - 1)] +
                                           yPreviousPart[i * size + (j + 1)]);

        norm = infiniteNorm(yPart, yPreviousPart, size, locationSize - size);

        MPI_Allreduce(&norm, &finalNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    } while (finalNorm > eps);

    if (processId == 0)
        timeEnd = MPI_Wtime();

    MPI_Gatherv(yPart.data() + receiveDisplacement, locationSize - size - extraSize, MPI_DOUBLE,
                y.data(), numbersOfProcessDataParts.data(), displacement.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (processId == 0)
    {
        std::cout << "\033[1;32mJacobi. MPI_Sendrecv\033[0m" << std::endl;
        std::cout << "Substraction norm: " << finalNorm << std::endl;
        std::cout << "Number of iterations: " << iterationsNumber << std::endl;
        std::cout << "Time: " << timeEnd - timeStart << std::endl;
        double differenceWithAccurateSolution = infiniteNorm(y, solution, 0, y.size());
        std::cout << "Diffrence: " << differenceWithAccurateSolution << std::endl;
    }
};

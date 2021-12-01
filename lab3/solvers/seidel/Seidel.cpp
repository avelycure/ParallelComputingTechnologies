#include "Seidel.hpp"

void seidelV1(std::vector<double> &y,
                        double h,
                        int size,
                        double k_square,
                        int processesNumber,
                        int processId,
                        double eps,
                        double &time,
                        InitialConditions initialConditions)
{
    int locationSize;
    int receiveDisplace;
    int extaSize;
    int offset;

    std::vector<double> vec_right;
    std::vector<int> len, disp;
    std::vector<double> yPart, yPreviousPart, vec_right_loc;

    divideVectorBetweenProcesses(y, h, size, k_square, processesNumber, processId, yPart, yPreviousPart, vec_right_loc, len, disp,
                                 locationSize, receiveDisplace, extaSize, offset);

    int iterationsNumber = 0;
    double norm;
    double finalNorm;

    double c = 1.0 / (4.0 + k_square);
    k_square /= (h * h);

    int destination = 0;
    int source = 0;

    double timeStart;
    double timeEnd;

    int sendCount = (processId - (processesNumber - 1)) ? size : 0;
    int receiveCount = processId ? size : 0;

    if (processesNumber > 1)
        setSourceAndDestination(processesNumber, processId, destination, source);

    if (processId == 0)
        timeStart = MPI_Wtime();

    do
    {
        iterationsNumber++;

        yPreviousPart.swap(yPart);

        MPI_Status statU, statL;

        MPI_Sendrecv(yPreviousPart.data() + locationSize - 2 * size, sendCount, MPI_DOUBLE, destination, 56,
                     yPreviousPart.data(), receiveCount, MPI_DOUBLE, source, 56, MPI_COMM_WORLD, &statU);

        MPI_Sendrecv(yPreviousPart.data() + size, receiveCount, MPI_DOUBLE, source, 65,
                     yPreviousPart.data() + locationSize - size, sendCount, MPI_DOUBLE, destination, 65, MPI_COMM_WORLD, &statU);

        for (int i = 1; i < locationSize / size - 1; ++i)
            for (int j = ((i + offset + 1) % 2) + 1; j < size - 1; j += 2)
                yPart[i * size + j] = c * (h * h * initialConditions.f((i + offset) * h, j * h, k_square) +
                                           yPreviousPart[(i - 1) * size + j] +
                                           yPreviousPart[(i + 1) * size + j] +
                                           yPreviousPart[i * size + (j - 1)] +
                                           yPreviousPart[i * size + (j + 1)]);

        MPI_Sendrecv(yPart.data() + locationSize - 2 * size, sendCount, MPI_DOUBLE, destination, 56,
                     yPart.data(), receiveCount, MPI_DOUBLE, source, 56, MPI_COMM_WORLD, &statU);

        MPI_Sendrecv(yPart.data() + size, receiveCount, MPI_DOUBLE, source, 65,
                     yPart.data() + locationSize - size, sendCount, MPI_DOUBLE, destination, 65, MPI_COMM_WORLD, &statU);

        for (int i = 1; i < locationSize / size - 1; ++i)
            for (int j = ((i + offset) % 2) + 1; j < size - 1; j += 2)
                yPart[i * size + j] = c * (h * h * initialConditions.f((i + offset) * h, j * h, k_square) +
                                           yPart[(i - 1) * size + j] +
                                           yPart[(i + 1) * size + j] +
                                           yPart[i * size + (j - 1)] +
                                           yPart[i * size + (j + 1)]);

        norm = infiniteNorm(yPart, yPreviousPart, size, locationSize - size);

        MPI_Allreduce(&norm, &finalNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    } while (finalNorm > eps);

    if (processId == 0)
        timeEnd = MPI_Wtime();

    MPI_Gatherv(yPart.data() + receiveDisplace, locationSize - size - extaSize,
                MPI_DOUBLE, y.data(), len.data(), disp.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (processId == 0)
    {
        std::cout << "\033[1;32mSeidel. MPI_Sendrecv\033[0m" << std::endl;
        std::cout << "Substraction norm: " << finalNorm << std::endl;
        std::cout << "Number of iterations: " << iterationsNumber << std::endl;
        std::cout << "Time: " << timeEnd - timeStart << std::endl;
    }

    time = timeEnd - timeStart;
}
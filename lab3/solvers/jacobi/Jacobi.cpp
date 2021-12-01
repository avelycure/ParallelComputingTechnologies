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
    double &time,
    InitialConditions initialConditions)
{
    int vectorPartSize;
    int receiveDisplacement;
    int extraSize;
    int offset;

    int iterationsNumber = 0;
    double norm;
    double finalNorm = 0.0;

    double c = 1.0 / (4.0 + kSquare);
    kSquare = kSquare / (h * h);

    int destination = 0;
    int source = 0;
    int send1 = 0;
    int recv1 = 0;
    int send2 = 0;
    int recv2 = 0;

    std::vector<double> rightPart;
    std::vector<int> numbersOfProcessDataParts;
    std::vector<int> displacement;
    std::vector<double> yPart;
    std::vector<double> yPreviousPart;
    std::vector<double> partOfRightPart;

    divideVectorBetweenProcesses(y, h, size, kSquare, numberOfProcesses, processId,
                                 yPart, yPreviousPart, partOfRightPart, numbersOfProcessDataParts, displacement,
                                 vectorPartSize, receiveDisplacement, extraSize, offset);

    if (numberOfProcesses > 1)
        setInteractionsScheme(numberOfProcesses, processId, destination, source, send1, recv1, send2, recv2);

    //Upper Lower?
    MPI_Status statU, statL;
    double timeStart;
    double timeEnd;

    if (processId == 0)
        timeStart = MPI_Wtime();

    do
    {
        iterationsNumber++;

        yPreviousPart.swap(yPart);

        // Block sending of data, page 12,16. Sends size * send1 elements to processId + 1 destination
        MPI_Send(yPreviousPart.data() + vectorPartSize - 2 * size, size * send1, MPI_DOUBLE, destination, 56, MPI_COMM_WORLD);
        MPI_Recv(yPreviousPart.data(), size * recv1, MPI_DOUBLE, source, 56, MPI_COMM_WORLD, &statU);
        MPI_Send(yPreviousPart.data() + vectorPartSize - 2 * size, size * send2, MPI_DOUBLE, destination, 57, MPI_COMM_WORLD);
        MPI_Recv(yPreviousPart.data(), size * recv2, MPI_DOUBLE, source, 57, MPI_COMM_WORLD, &statU);

        MPI_Send(yPreviousPart.data() + size, size * recv2, MPI_DOUBLE, source, 65, MPI_COMM_WORLD);
        MPI_Recv(yPreviousPart.data() + vectorPartSize - size, size * send2, MPI_DOUBLE, destination, 65, MPI_COMM_WORLD, &statL);
        MPI_Send(yPreviousPart.data() + size, size * recv1, MPI_DOUBLE, source, 67, MPI_COMM_WORLD);
        MPI_Recv(yPreviousPart.data() + vectorPartSize - size, size * send1, MPI_DOUBLE, destination, 67, MPI_COMM_WORLD, &statL);

        for (int i = 1; i < vectorPartSize / size - 1; ++i)
            for (int j = 1; j < size - 1; ++j)
                yPart[i * size + j] = c * (h * h * initialConditions.f((i + offset) * h, j * h, kSquare) +
                                           yPreviousPart[(i - 1) * size + j] +
                                           yPreviousPart[(i + 1) * size + j] +
                                           yPreviousPart[i * size + (j - 1)] +
                                           yPreviousPart[i * size + (j + 1)]);

        norm = infiniteNorm(yPart, yPreviousPart, receiveDisplacement, vectorPartSize - size);

        // one time find max from norm in all processes
        MPI_Allreduce(&norm, &finalNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    } while (finalNorm > eps);

    if (processId == 0)
        timeEnd = MPI_Wtime();

    MPI_Gatherv(yPart.data() + receiveDisplacement, vectorPartSize - size - extraSize, MPI_DOUBLE, y.data(),
                numbersOfProcessDataParts.data(), displacement.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (processId == 0)
    {
        std::cout << "\033[1;32mJacobi. MPI_Send. MPI_Recv\033[0m" << std::endl;
        std::cout << "Substraction norm: " << finalNorm << std::endl;
        std::cout << "Number of iterations: " << iterationsNumber << std::endl;
        std::cout << "Time: " << timeEnd - timeStart << std::endl;
        //for (int i = 0; i < size; ++i) {
        //    for (int j = 0; j < size; ++j) {
        //        std::cout << y[i * size + j] << "\t";
        //    }
        //    std::cout << "\n";
        //}
    }
    time = timeEnd - timeStart;
};

/**
 * Used: MPI_Sendrecv
 * */
void jacobiV2(std::vector<double> &y,
              double h,
              int size,
              double kSquare,
              int processesNumber,
              int processId,
              double eps,
              double &time,
              InitialConditions initialConditions)
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
    kSquare = kSquare / (h * h);

    int destination = 0;
    int source = 0;

    int sendCount = (processId - (processesNumber - 1)) ? size : 0;
    int receiveCount = processId ? size : 0;

    double timeStart;
    double timeEnd;

    divideVectorBetweenProcesses(y, h, size, kSquare, processesNumber, processId, yPart, yPreviousPart, partOfRightPart, numbersOfProcessDataParts, displacement,
                                 locationSize, receiveDisplacement, extraSize, offset);

    if (processesNumber > 1)
        setSourceAndDestination(processesNumber, processId, destination, source);

    if (processId == 0)
        timeStart = MPI_Wtime();

    do
    {
        iterationsNumber++;

        yPreviousPart.swap(yPart);

        //typo in 177??? statU or statL?
        MPI_Status statU, statL;

        //Block sending, page 33
        MPI_Sendrecv(yPreviousPart.data() + locationSize - 2 * size, sendCount, MPI_DOUBLE, destination, 56,
                     yPreviousPart.data(), receiveCount, MPI_DOUBLE, source, 56, MPI_COMM_WORLD, &statU);

        MPI_Sendrecv(yPreviousPart.data() + size, receiveCount, MPI_DOUBLE, source, 65,
                     yPreviousPart.data() + locationSize - size, sendCount, MPI_DOUBLE, destination, 65, MPI_COMM_WORLD, &statL);

        for (int i = 1; i < locationSize / size - 1; ++i)
            for (int j = 1; j < size - 1; ++j)
                yPart[i * size + j] = c * (h * h * ((i + offset) * h, j * h, kSquare) + yPreviousPart[(i - 1) * size + j] +
                                           yPreviousPart[(i + 1) * size + j] +
                                           yPreviousPart[i * size + (j - 1)] +
                                           yPreviousPart[i * size + (j + 1)]);

        norm = infiniteNorm(yPart, yPreviousPart, size, locationSize - size);

        MPI_Allreduce(&norm, &finalNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    } while (finalNorm > eps);

    if (processId == 0)
        timeEnd = MPI_Wtime();

    MPI_Gatherv(yPart.data() + receiveDisplacement, locationSize - size - extraSize, MPI_DOUBLE, y.data(), numbersOfProcessDataParts.data(), displacement.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (processId == 0)
    {
        std::cout << "\033[1;32mJacobi. MPI_Sendrecv\033[0m" << std::endl;
        std::cout << "Substraction norm: " << finalNorm << std::endl;
        std::cout << "Number of iterations: " << iterationsNumber << std::endl;
        std::cout << "Time: " << timeEnd - timeStart << std::endl;
    }

    time = timeEnd - timeStart;
};

void jacobiV3(std::vector<double> &y,
              double h,
              int size,
              double kSquare,
              int processesNumber,
              int processId,
              double eps,
              double &time,
              InitialConditions initialConditions)
{
    int locationSize;
    int receiceDisplacement;
    int extraSize;
    int offset;

    std::vector<double> rightPart;
    std::vector<int> numbersOfProcessDataParts;
    std::vector<int> displace;
    std::vector<double> yPart;
    std::vector<double> yPreviousPart;
    std::vector<double> partOfRightPart;

    divideVectorBetweenProcesses(y, h, size, kSquare, processesNumber, processId, yPart, yPreviousPart, partOfRightPart, numbersOfProcessDataParts, displace,
                                 locationSize, receiceDisplacement, extraSize, offset);

    int iterationsNumber = 0;
    double norm;
    double finalNorm;

    double c = 1.0 / (4.0 + kSquare);
    kSquare /= (h * h);

    double timeStart;
    double timeEnd;
    if (processId == 0)
        timeStart = MPI_Wtime();

    std::vector<MPI_Request> reqInEven;
    std::vector<MPI_Request> reqOutEven;
    std::vector<MPI_Request> reqInOdd;
    std::vector<MPI_Request> reqOutOdd;

    if ((processId == 0) || (processId == processesNumber - 1))
    {
        reqInEven.resize(1);
        reqOutEven.resize(1);
        reqInOdd.resize(1);
        reqOutOdd.resize(1);
    }
    else
    {
        reqInEven.resize(2);
        reqOutEven.resize(2);
        reqInOdd.resize(2);
        reqOutOdd.resize(2);
    }
    int nInEven = 0;
    int nOutEven = 0;
    int nInOdd = 0;
    int nOutOdd = 0;

    if (processId != processesNumber - 1)
    {
        // page 30, non blocking sending
        MPI_Send_init(yPreviousPart.data() + locationSize - 2 * size, size, MPI_DOUBLE, processId + 1, 56, MPI_COMM_WORLD, reqOutEven.data() + nOutEven);
        MPI_Recv_init(yPreviousPart.data() + locationSize - size, size, MPI_DOUBLE, processId + 1, 65, MPI_COMM_WORLD, reqInEven.data() + nInEven);
        nOutEven++;
        nInEven++;
        MPI_Send_init(yPart.data() + locationSize - 2 * size, size, MPI_DOUBLE, processId + 1, 56, MPI_COMM_WORLD, reqOutOdd.data() + nOutOdd);
        MPI_Recv_init(yPart.data() + locationSize - size, size, MPI_DOUBLE, processId + 1, 65, MPI_COMM_WORLD, reqInOdd.data() + nInOdd);
        nOutOdd++;
        nInOdd++;
    }
    if (processId != 0)
    {
        MPI_Send_init(yPreviousPart.data() + size, size, MPI_DOUBLE, processId - 1, 65, MPI_COMM_WORLD, reqOutEven.data() + nOutEven);
        MPI_Recv_init(yPreviousPart.data(), size, MPI_DOUBLE, processId - 1, 56, MPI_COMM_WORLD, reqInEven.data() + nInEven);
        nOutEven++;
        nInEven++;
        MPI_Send_init(yPart.data() + size, size, MPI_DOUBLE, processId - 1, 65, MPI_COMM_WORLD, reqOutOdd.data() + nOutOdd);
        MPI_Recv_init(yPart.data(), size, MPI_DOUBLE, processId - 1, 56, MPI_COMM_WORLD, reqInOdd.data() + nInOdd);
        nOutOdd++;
        nInOdd++;
    }

    std::vector<MPI_Status> statInEven(nInEven);
    std::vector<MPI_Status> statOutEven(nOutEven);
    std::vector<MPI_Status> statInOdd(nInOdd);
    std::vector<MPI_Status> statOutOdd(nOutOdd);

    do
    {
        iterationsNumber++;

        yPreviousPart.swap(yPart);

        if (iterationsNumber % 2 != 0)
        {
            MPI_Startall(nOutOdd, reqOutOdd.data());
            MPI_Startall(nInOdd, reqInOdd.data());
        }
        else
        {
            MPI_Startall(nOutEven, reqOutEven.data());
            MPI_Startall(nInEven, reqInEven.data());
        }

        for (int i = 2; i < locationSize / size - 2; ++i)
            for (int j = 1; j < size - 1; ++j)
                yPart[i * size + j] = c * (h * h * initialConditions.f((i + offset) * h, j * h, kSquare) +
                                           yPreviousPart[(i - 1) * size + j] +
                                           yPreviousPart[(i + 1) * size + j] +
                                           yPreviousPart[i * size + (j - 1)] +
                                           yPreviousPart[i * size + (j + 1)]);

        if (iterationsNumber % 2 != 0)
        {
            MPI_Waitall(nInOdd, reqInOdd.data(), statInOdd.data());
            MPI_Waitall(nOutOdd, reqOutOdd.data(), statInOdd.data());
        }
        else
        {
            MPI_Waitall(nInEven, reqInEven.data(), statInEven.data());
            MPI_Waitall(nOutEven, reqOutEven.data(), statInEven.data());
        }

        int index = 1;
        for (int j = 1; j < size - 1; ++j)
            yPart[index * size + j] = c * (h * h * initialConditions.f((index + offset) * h, j * h, kSquare) +
                                           yPreviousPart[(index - 1) * size + j] +
                                           yPreviousPart[(index + 1) * size + j] +
                                           yPreviousPart[index * size + (j - 1)] +
                                           yPreviousPart[index * size + (j + 1)]);

        index = locationSize / size - 2;
        for (int j = 1; j < size - 1; ++j)
            yPart[index * size + j] = c * (h * h * initialConditions.f((index + offset) * h, j * h, kSquare) +
                                           yPreviousPart[(index - 1) * size + j] +
                                           yPreviousPart[(index + 1) * size + j] +
                                           yPreviousPart[index * size + (j - 1)] +
                                           yPreviousPart[index * size + (j + 1)]);

        norm = infiniteNorm(yPart, yPreviousPart, size, locationSize - size);

        MPI_Allreduce(&norm, &finalNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    } while (finalNorm > eps);

    if (processId == 0)
        timeEnd = MPI_Wtime();

    MPI_Gatherv(yPart.data() + receiceDisplacement, locationSize - size - extraSize, MPI_DOUBLE, y.data(), numbersOfProcessDataParts.data(), displace.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (processId == 0)
    {
        std::cout << "\033[1;32mJacobi. MPI_Send_init. MPI_Recv_init\033[0m" << std::endl;
        std::cout << "Substraction norm: " << finalNorm << std::endl;
        std::cout << "Number of iterations: " << iterationsNumber << std::endl;
        std::cout << "Time: " << timeEnd - timeStart << std::endl;
    }

    time = timeEnd - timeStart;
}
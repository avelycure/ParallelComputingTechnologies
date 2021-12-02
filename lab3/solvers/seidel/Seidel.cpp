#include "Seidel.hpp"

void seidelV1(std::vector<double> &y,
              double h,
              int size,
              double kSquare,
              int processesNumber,
              int processId,
              double eps,
              InitialConditions initialConditions)
{
    int locationSize;
    int receiveDisplace;
    int extaSize;
    int offset;

    std::vector<double> rightPart;
    std::vector<int> numbersOfProcessDataParts, displace;
    std::vector<double> yPart, yPreviousPart, partOfRightPart;

    divideVectorBetweenProcesses(y, h, size, kSquare, processesNumber, processId, yPart, yPreviousPart, partOfRightPart, numbersOfProcessDataParts, displace,
                                 locationSize, receiveDisplace, extaSize, offset);

    int iterationsNumber = 0;
    double norm;
    double finalNorm;

    double c = 1.0 / (4.0 + kSquare);
    kSquare /= (h * h);

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
                yPart[i * size + j] = c * (h * h * initialConditions.f((i + offset) * h, j * h, kSquare) +
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
                yPart[i * size + j] = c * (h * h * initialConditions.f((i + offset) * h, j * h, kSquare) +
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
                MPI_DOUBLE, y.data(), numbersOfProcessDataParts.data(), displace.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (processId == 0)
    {
        std::cout << "\033[1;32mSeidel. MPI_Sendrecv\033[0m" << std::endl;
        std::cout << "Substraction norm: " << finalNorm << std::endl;
        std::cout << "Number of iterations: " << iterationsNumber << std::endl;
        std::cout << "Time: " << timeEnd - timeStart << std::endl;
    }
}

void seidelV2(std::vector<double> &y,
              double h,
              int size,
              double kSquare,
              int processesNumber,
              int processId,
              double eps,
              InitialConditions initialConditions)
{
    int location_size;
    int receiveDisplace;
    int extraSize;
    int offset;

    std::vector<double> rightPart;
    std::vector<int> numbersOfProcessDataParts;
    std::vector<int> displace;
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
    int sendCount = (processId - (processesNumber - 1)) ? size : 0;
    int receiveCount = processId ? size : 0;

    double timeStart;
    double timeEnd;

    divideVectorBetweenProcesses(y, h, size, kSquare, processesNumber, processId, yPart, yPreviousPart, partOfRightPart, numbersOfProcessDataParts, displace,
                                 location_size, receiveDisplace, extraSize, offset);

    if (processesNumber > 1)
        setSourceAndDestination(processesNumber, processId, destination, source);

    if (processId == 0)
        timeStart = MPI_Wtime();

    do
    {
        iterationsNumber++;

        yPreviousPart.swap(yPart);

        MPI_Status statU, statL;

        MPI_Sendrecv(yPreviousPart.data() + location_size - 2 * size, sendCount, MPI_DOUBLE, destination, 56,
                     yPreviousPart.data(), receiveCount, MPI_DOUBLE, source, 56, MPI_COMM_WORLD, &statU);

        MPI_Sendrecv(yPreviousPart.data() + size, receiveCount, MPI_DOUBLE, source, 65,
                     yPreviousPart.data() + location_size - size, sendCount, MPI_DOUBLE, destination, 65, MPI_COMM_WORLD, &statU);

        for (int i = 1; i < location_size / size - 1; ++i)
            for (int j = ((i + offset + 1) % 2) + 1; j < size - 1; j += 2)
                yPart[i * size + j] = c * (h * h * initialConditions.f((i + offset) * h, j * h, kSquare) +
                                           yPreviousPart[(i - 1) * size + j] +
                                           yPreviousPart[(i + 1) * size + j] +
                                           yPreviousPart[i * size + (j - 1)] +
                                           yPreviousPart[i * size + (j + 1)]);

        MPI_Sendrecv(yPart.data() + location_size - 2 * size, sendCount, MPI_DOUBLE, destination, 56,
                     yPart.data(), receiveCount, MPI_DOUBLE, source, 56, MPI_COMM_WORLD, &statU);

        MPI_Sendrecv(yPart.data() + size, receiveCount, MPI_DOUBLE, source, 65,
                     yPart.data() + location_size - size, sendCount, MPI_DOUBLE, destination, 65, MPI_COMM_WORLD, &statU);

        for (int i = 1; i < location_size / size - 1; ++i)
            for (int j = ((i + offset) % 2) + 1; j < size - 1; j += 2)
                yPart[i * size + j] = c * (h * h * initialConditions.f((i + offset) * h, j * h, kSquare) +
                                           yPart[(i - 1) * size + j] +
                                           yPart[(i + 1) * size + j] +
                                           yPart[i * size + (j - 1)] +
                                           yPart[i * size + (j + 1)]);

        norm = infiniteNorm(yPart, yPreviousPart, size, location_size - size);

        MPI_Allreduce(&norm, &finalNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    } while (finalNorm > eps);

    if (processId == 0)
        timeEnd = MPI_Wtime();

    MPI_Gatherv(yPart.data() + receiveDisplace, location_size - size - extraSize, MPI_DOUBLE, y.data(), numbersOfProcessDataParts.data(), displace.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (processId == 0)
    {
        std::cout << "\033[1;32mSeidel. MPI_Sendrecv\033[0m" << std::endl;
        std::cout << "Substraction norm: " << finalNorm << std::endl;
        std::cout << "Number of iterations: " << iterationsNumber << std::endl;
        std::cout << "Time: " << timeEnd - timeStart << std::endl;
    }
}

void seidelV3(std::vector<double> &y,
              double h,
              int size,
              double kSquare,
              int processesNumber,
              int processId,
              double eps,
              InitialConditions initialConditions)
{
    int locationSize;
    int receiveDisplace;
    int extraSize;
    int offset;

    std::vector<double> vec_right;
    std::vector<int> numbersOfProcessDataParts, disp;
    std::vector<double> yPart, yPreviousPart, vec_right_loc;

    divideVectorBetweenProcesses(y, h, size, kSquare, processesNumber, processId, yPart, yPreviousPart, vec_right_loc, numbersOfProcessDataParts, disp,
                                 locationSize, receiveDisplace, extraSize, offset);

    int iterationsNumber = 0;
    double norm;
    double finalNorm;

    double h_square = h * h;
    double c = 1.0 / (4.0 + kSquare);
    kSquare /= h_square;

    double timeStart;
    double timeEnd;
    if (processId == 0)
        timeStart = MPI_Wtime();

    std::vector<MPI_Request> reqInBlack, reqOutBlack, reqInRed, reqOutRed;
    std::vector<MPI_Request> reqInBlack2, reqOutBlack2, reqInRed2, reqOutRed2;
    if ((processId == 0) || (processId == processesNumber - 1))
    {
        reqInBlack.resize(1);
        reqOutBlack.resize(1);
        reqInRed.resize(1);
        reqOutRed.resize(1);
        reqInBlack2.resize(1);
        reqOutBlack2.resize(1);
        reqInRed2.resize(1);
        reqOutRed2.resize(1);
    }
    else
    {
        reqInBlack.resize(2);
        reqOutBlack.resize(2);
        reqInRed.resize(2);
        reqOutRed.resize(2);
        reqInBlack2.resize(2);
        reqOutBlack2.resize(2);
        reqInRed2.resize(2);
        reqOutRed2.resize(2);
    }

    int nInBlack = 0, nOutBlack = 0, nInRed = 0, nOutRed = 0;
    int nInBlack2 = 0, nOutBlack2 = 0, nInRed2 = 0, nOutRed2 = 0;

    // black
    if (processId != processesNumber - 1)
    {
        MPI_Send_init(yPreviousPart.data() + locationSize - 2 * size, size, MPI_DOUBLE, processId + 1, 56, MPI_COMM_WORLD, reqOutBlack.data() + nOutBlack);
        MPI_Recv_init(yPreviousPart.data() + locationSize - size, size, MPI_DOUBLE, processId + 1, 65, MPI_COMM_WORLD, reqInBlack.data() + nInBlack);
        nOutBlack++;
        nInBlack++;
        MPI_Send_init(yPart.data() + locationSize - 2 * size, size, MPI_DOUBLE, processId + 1, 57, MPI_COMM_WORLD, reqOutBlack2.data() + nOutBlack2);
        MPI_Recv_init(yPart.data() + locationSize - size, size, MPI_DOUBLE, processId + 1, 67, MPI_COMM_WORLD, reqInBlack2.data() + nInBlack2);
        nOutBlack2++;
        nInBlack2++;
    }
    if (processId != 0)
    {
        MPI_Send_init(yPreviousPart.data() + size, size, MPI_DOUBLE, processId - 1, 65, MPI_COMM_WORLD, reqOutBlack.data() + nOutBlack);
        MPI_Recv_init(yPreviousPart.data(), size, MPI_DOUBLE, processId - 1, 56, MPI_COMM_WORLD, reqInBlack.data() + nInBlack);
        nOutBlack++;
        nInBlack++;
        MPI_Send_init(yPart.data() + size, size, MPI_DOUBLE, processId - 1, 67, MPI_COMM_WORLD, reqOutBlack2.data() + nOutBlack2);
        MPI_Recv_init(yPart.data(), size, MPI_DOUBLE, processId - 1, 57, MPI_COMM_WORLD, reqInBlack2.data() + nInBlack2);
        nOutBlack2++;
        nInBlack2++;
    }

    //red
    if (processId != processesNumber - 1)
    {
        MPI_Send_init(yPart.data() + locationSize - 2 * size, size, MPI_DOUBLE, processId + 1, 56, MPI_COMM_WORLD, reqOutRed.data() + nOutRed);
        MPI_Recv_init(yPart.data() + locationSize - size, size, MPI_DOUBLE, processId + 1, 65, MPI_COMM_WORLD, reqInRed.data() + nInRed);
        nOutRed++;
        nInRed++;
        MPI_Send_init(yPreviousPart.data() + locationSize - 2 * size, size, MPI_DOUBLE, processId + 1, 58, MPI_COMM_WORLD, reqOutRed2.data() + nOutRed2);
        MPI_Recv_init(yPreviousPart.data() + locationSize - size, size, MPI_DOUBLE, processId + 1, 68, MPI_COMM_WORLD, reqInRed2.data() + nInRed2);
        nOutRed2++;
        nInRed2++;
    }
    if (processId != 0)
    {
        MPI_Send_init(yPart.data() + size, size, MPI_DOUBLE, processId - 1, 65, MPI_COMM_WORLD, reqOutRed.data() + nOutRed);
        MPI_Recv_init(yPart.data(), size, MPI_DOUBLE, processId - 1, 56, MPI_COMM_WORLD, reqInRed.data() + nInRed);
        nOutRed++;
        nInRed++;
        MPI_Send_init(yPreviousPart.data() + size, size, MPI_DOUBLE, processId - 1, 68, MPI_COMM_WORLD, reqOutRed2.data() + nOutRed2);
        MPI_Recv_init(yPreviousPart.data(), size, MPI_DOUBLE, processId - 1, 58, MPI_COMM_WORLD, reqInRed2.data() + nInRed2);
        nOutRed2++;
        nInRed2++;
    }

    std::vector<MPI_Status> statInBlack(nInBlack), statOutBlack(nOutBlack), statInRed(nInRed), statOutRed(nOutRed);
    std::vector<MPI_Status> statInBlack2(nInBlack2), statOutBlack2(nOutBlack2), statInRed2(nInRed2), statOutRed2(nOutRed2);

    do
    {
        iterationsNumber++;

        yPreviousPart.swap(yPart);

        if (iterationsNumber % 2 != 0)
        {
            MPI_Startall(nOutBlack2, reqOutBlack2.data());
            MPI_Startall(nInBlack2, reqInBlack2.data());
        }
        else
        {
            MPI_Startall(nOutBlack, reqOutBlack.data());
            MPI_Startall(nInBlack, reqInBlack.data());
        }

        for (int i = 2; i < locationSize / size - 2; ++i)
            for (int j = ((i + offset + 1) % 2) + 1; j < size - 1; j += 2)
                yPart[i * size + j] = c * (h_square * initialConditions.f((i + offset) * h, j * h, kSquare) +
                                           yPreviousPart[(i - 1) * size + j] +
                                           yPreviousPart[(i + 1) * size + j] +
                                           yPreviousPart[i * size + (j - 1)] +
                                           yPreviousPart[i * size + (j + 1)]);

        if (iterationsNumber % 2 != 0)
        {
            MPI_Waitall(nOutBlack2, reqOutBlack2.data(), statOutBlack2.data());
            MPI_Waitall(nInBlack2, reqInBlack2.data(), statInBlack2.data());
        }
        else
        {
            MPI_Waitall(nInBlack, reqInBlack.data(), statInBlack.data());
            MPI_Waitall(nOutBlack, reqOutBlack.data(), statOutBlack.data());
        }

        {
            int i = 1;
            for (int j = ((i + offset + 1) % 2) + 1; j < size - 1; j += 2)
                yPart[i * size + j] = c * (h_square * initialConditions.f((i + offset) * h, j * h, kSquare) +
                                           yPreviousPart[(i - 1) * size + j] +
                                           yPreviousPart[(i + 1) * size + j] +
                                           yPreviousPart[i * size + (j - 1)] +
                                           yPreviousPart[i * size + (j + 1)]);
        }

        {
            int i = locationSize / size - 2;
            for (int j = ((i + offset + 1) % 2) + 1; j < size - 1; j += 2)
                yPart[i * size + j] = c * (h_square * initialConditions.f((i + offset) * h, j * h, kSquare) +
                                           yPreviousPart[(i - 1) * size + j] +
                                           yPreviousPart[(i + 1) * size + j] +
                                           yPreviousPart[i * size + (j - 1)] +
                                           yPreviousPart[i * size + (j + 1)]);
        }

        if (iterationsNumber % 2 != 0)
        {
            MPI_Startall(nOutRed2, reqOutRed2.data());
            MPI_Startall(nInRed2, reqInRed2.data());
        }
        else
        {
            MPI_Startall(nOutRed, reqOutRed.data());
            MPI_Startall(nInRed, reqInRed.data());
        }

        for (int i = 2; i < locationSize / size - 2; ++i)
            for (int j = ((i + offset) % 2) + 1; j < size - 1; j += 2)
                yPart[i * size + j] = c * (h_square * initialConditions.f((i + offset) * h, j * h, kSquare) +
                                           yPart[(i - 1) * size + j] +
                                           yPart[(i + 1) * size + j] +
                                           yPart[i * size + (j - 1)] +
                                           yPart[i * size + (j + 1)]);

        if (iterationsNumber % 2 != 0)
        {
            MPI_Waitall(nInRed2, reqInRed2.data(), statInRed2.data());
            MPI_Waitall(nOutRed2, reqOutRed2.data(), statOutRed2.data());
        }
        else
        {
            MPI_Waitall(nInRed, reqInRed.data(), statInRed.data());
            MPI_Waitall(nOutRed, reqOutRed.data(), statOutRed.data());
        }

        {
            int i = 1;
            for (int j = ((i + offset) % 2) + 1; j < size - 1; j += 2)
                yPart[i * size + j] = c * (h_square * initialConditions.f((i + offset) * h, j * h, kSquare) +
                                           yPart[(i - 1) * size + j] +
                                           yPart[(i + 1) * size + j] +
                                           yPart[i * size + (j - 1)] +
                                           yPart[i * size + (j + 1)]);
        }

        {
            int i = locationSize / size - 2;
            for (int j = ((i + offset) % 2) + 1; j < size - 1; j += 2)
                yPart[i * size + j] = c * (h_square * initialConditions.f((i + offset) * h, j * h, kSquare) +
                                           yPart[(i - 1) * size + j] +
                                           yPart[(i + 1) * size + j] +
                                           yPart[i * size + (j - 1)] +
                                           yPart[i * size + (j + 1)]);
        }

        norm = infiniteNorm(yPart, yPreviousPart, size, locationSize - size);

        MPI_Allreduce(&norm, &finalNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    } while (finalNorm > eps);

    if (processId == 0)
        timeEnd = MPI_Wtime();

    MPI_Gatherv(yPart.data() + receiveDisplace, locationSize - size - extraSize, MPI_DOUBLE, y.data(), numbersOfProcessDataParts.data(), disp.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (processId == 0)
    {
        std::cout << "\033[1;32mSeidel. MPI_Send_init. MPI_Recv_init\033[0m" << std::endl;
        std::cout << "Substraction norm: " << finalNorm << std::endl;
        std::cout << "Number of iterations: " << iterationsNumber << std::endl;
        std::cout << "Time: " << timeEnd - timeStart << std::endl;
    }
}
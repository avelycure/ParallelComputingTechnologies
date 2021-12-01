#include "Common.hpp"

void printLog(
    int numberOfProcesses,
    int sumSize,
    double khRelation,
    double eps)
{
    std::cout << "\033[1;33mLOG\033[0m" << std::endl;
    std::cout << "numberOfProcesses = " << numberOfProcesses << std::endl;
    std::cout << "sumSize = " << sumSize << std::endl;
    std::cout << "k^2 / h^2 = " << khRelation << std::endl;
    std::cout << "eps = " << eps << std::endl;
}

double infiniteNorm(std::vector<double> &vec1, std::vector<double> &vec2, int begin, int end)
{
    double norm = 0.0;
    double elem = 0.0;
    for (int i = begin; i < end; ++i)
    {
        elem = fabs(vec1[i] - vec2[i]);
        if (norm < elem)
            norm = elem;
    }
    return norm;
}

/**
 * Set what parts of initial vector will work with every process
 * @offest and @displacement is all about offset, but @offset is offset in rows
 * while displacement is offset in elements
 * */
void divideVectorBetweenProcesses(std::vector<double> &y,
                                  double h,
                                  int size,
                                  double kSquare,
                                  int processesNumber,
                                  int processId,
                                  std::vector<double> &yPart,
                                  std::vector<double> &yPreviousPart,
                                  std::vector<double> &partOfRightPart,
                                  //integer array of number of elements sended to each process
                                  std::vector<int> &numbersOfProcessDataParts,
                                  //integer array of offsets in relation of y beginning, page 41
                                  std::vector<int> &displacement,
                                  int &vectorPartSize,
                                  //????what for
                                  int &receiveDisplacement,
                                  int &extraSize,
                                  int &offset)
{

    std::vector<int> vecOffset(processesNumber);

    if (processId == 0)
    {
        numbersOfProcessDataParts.resize(processesNumber);
        displacement.resize(processesNumber);

        // divide initial vector on parts
        for (int i = 0; i < processesNumber; i++)
            numbersOfProcessDataParts[i] = (size / processesNumber) * size;

        // case of not integer division
        for (int i = 0; i < size - (size / processesNumber) * (processesNumber); i++)
            numbersOfProcessDataParts[i] += size;

        displacement[0] = 0;
        vecOffset[0] = -1;

        //set displacement and offset for all processes
        for (int i = 0; i < processesNumber - 1; i++)
        {
            displacement[i + 1] = displacement[i] + numbersOfProcessDataParts[i];
            vecOffset[i + 1] = vecOffset[i] + numbersOfProcessDataParts[i] / size;
        }
        vecOffset[0] = 0;
    }

    // this procedures send 1 element of MPI_INT type from led.data/vecOffset from ROOT process
    // in loc_size/offset variables of MPI_WORLD commutator
    MPI_Scatter(numbersOfProcessDataParts.data(), 1, MPI_INT, &vectorPartSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(vecOffset.data(), 1, MPI_INT, &offset, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (processId == 0)
    {
        y.resize(size * size);
        fillVectorWithZeros(y);
    }

    // may be this also handles cases when we have not equal parts of the matrix
    if (processesNumber > 1)
    {
        if ((processId == 0) || (processId == processesNumber - 1))
        {
            vectorPartSize += size;
            extraSize = 0.0;
        }
        else
        {
            vectorPartSize += 2 * size;
            extraSize = size;
        }
    }
    else
        extraSize = -size;

    yPart.resize(vectorPartSize);
    yPreviousPart.resize(vectorPartSize);
    partOfRightPart.resize(vectorPartSize);

    receiveDisplacement = (processId == 0) ? 0 : size;

    //Send data from root process
    MPI_Scatterv(y.data(),
                 numbersOfProcessDataParts.data(),
                 displacement.data(),
                 MPI_DOUBLE,
                 yPart.data() + receiveDisplacement,
                 vectorPartSize - size - extraSize,
                 MPI_DOUBLE,
                 0,
                 MPI_COMM_WORLD);
}

void fillVectorWithZeros(std::vector<double> &y)
{
    for (int i = 0; i < y.size(); i++)
        y[i] = 0.0;
}

void setInteractionsScheme(
    int processesNumber,
    const int processId,
    int &destination,
    int &source,
    int &send1,
    int &recv1,
    int &send2,
    int &recv2)
{

    setSourceAndDestination(processesNumber, processId, destination, source);

    if (processId % 2 == 0)
    {
        send1 = 1;
        recv1 = 0;
    }
    else
    {
        send1 = 0;
        recv1 = 1;
    }

    if (processId == processesNumber - 1)
        send1 = 0;

    recv2 = send1;

    if (processId == 0)
        recv2 = 0;

    if ((processId == processesNumber - 1) && (processesNumber % 2 != 0))
        recv2 = 1;

    send2 = recv1;
    if ((processId == processesNumber - 1) && (processesNumber % 2 == 0))
        send2 = 0;
}

void setSourceAndDestination(const int processesNumber, const int processId, int &destination, int &source)
{
    destination = processId + 1;

    source = processId - 1;

    if (processId == 0)
        source = processesNumber - 1;
    else if (processId == processesNumber - 1)
        destination = 0;
}
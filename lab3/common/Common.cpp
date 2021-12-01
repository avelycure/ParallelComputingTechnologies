#include "Common.hpp"

void printLog(
    int numberOfProcesses,
    int sumSize,
    double khRelation,
    double eps)
{
    std::cout << "LOG" << std::endl;
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

void initVecLoc(std::vector<double> &y,
                std::vector<double> &vec_right,
                double h,
                int size,
                double k_square,
                int processesNumber,
                int processId,
                std::vector<double> &y_loc,
                std::vector<double> &y_loc_prev,
                std::vector<double> &vec_right_loc,
                std::vector<int> &len,
                std::vector<int> &disp,
                int &loc_size,
                int &recv_disp,
                int &extr_size,
                int &offset)
{

    std::vector<int> vec_offset(processesNumber);

    if (processId == 0)
    {
        len.resize(processesNumber);
        disp.resize(processesNumber);

        for (int i = 0; i < processesNumber; ++i)
            len[i] = (size / processesNumber) * size;

        for (int i = 0; i < size - (size / processesNumber) * (processesNumber); ++i)
            len[i] += size;

        disp[0] = 0;
        vec_offset[0] = -1;
        //if (np > 1)
        for (int i = 0; i < processesNumber - 1; ++i)
        {
            disp[i + 1] = disp[i] + len[i];
            vec_offset[i + 1] = vec_offset[i] + len[i] / size;
        }
        vec_offset[0] = 0;
    }
    //MPI_Barrier(MPI_COMM_WORLD);

    // this procedure sends data, page 40
    MPI_Scatter(len.data(), 1, MPI_INT, &loc_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(vec_offset.data(), 1, MPI_INT, &offset, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (processId == 0)
    {
        y.resize(size * size);
        //vec_right.resize(size * size);
        fillVectorWithZeros(y);
        //comp_vec_right(vec_right, size, h, k_square);
    }

    if (processesNumber > 1)
    {
        if ((processId == 0) || (processId == processesNumber - 1))
        {
            loc_size += size;
            extr_size = 0.0;
        }
        else
        {
            loc_size += 2 * size;
            extr_size = size;
        }
    }
    else
        extr_size = -size;

    y_loc.resize(loc_size);
    y_loc_prev.resize(loc_size);
    vec_right_loc.resize(loc_size);

    recv_disp = (processId == 0) ? 0.0 : size;

    MPI_Scatterv(y.data(), len.data(), disp.data(), MPI_DOUBLE, y_loc.data() + recv_disp, loc_size - size - extr_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Scatterv(vec_right.data() + size, len.data(), disp.data(), MPI_DOUBLE, vec_right_loc.data() + size, loc_size - 2 * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void fillVectorWithZeros(std::vector<double> &y)
{
    for (int i = 0; i < y.size(); i++)
        y[i] = 0.0;
}

void send_recv_scheme(
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

    //fprintf(stdout, "myid = %d, source = %d, dest = %d, send1 = %d, recv1 = %d, send2 = %d, recv2 = %d", myid, source, dest, send1, recv1, send2, recv2);
    //fflush(stdout);
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
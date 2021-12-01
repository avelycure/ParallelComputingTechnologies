#include "Jacobi.hpp"

/**
 * Jacobi method of solving system. Result of this function is in y vector.
 * */
void jacobiSendReceive(
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
    int loc_size;
    int recv_disp;
    int extr_size;
    int offset;

    int iterationsNumber = 0;
    double norm;
    double finalNorm = 0.0;

    double c = 1.0 / (4.0 + kSquare);
    kSquare = kSquare / (h * h);

    int dest = 0, source = 0;
    int send1 = 0, recv1 = 0, send2 = 0, recv2 = 0;

    std::vector<double> vec_right;
    std::vector<int> len;
    std::vector<int> disp;
    std::vector<double> y_loc;
    std::vector<double> y_loc_prev;
    std::vector<double> vec_right_loc;

    initVecLoc(y, vec_right, h, size, kSquare, numberOfProcesses, processId, y_loc, y_loc_prev, vec_right_loc, len, disp,
               loc_size, recv_disp, extr_size, offset);

    if (numberOfProcesses > 1)
        send_recv_scheme(numberOfProcesses, processId, dest, source, send1, recv1, send2, recv2);

    MPI_Status statU, statL;
    double timeStart;
    double timeEnd;

    if (processId == 0)
        timeStart = MPI_Wtime();

    do
    {
        iterationsNumber++;

        y_loc_prev.swap(y_loc);

        //Send Recv
        MPI_Send(y_loc_prev.data() + loc_size - 2 * size, size * send1, MPI_DOUBLE, dest, 56, MPI_COMM_WORLD);
        MPI_Recv(y_loc_prev.data(), size * recv1, MPI_DOUBLE, source, 56, MPI_COMM_WORLD, &statU);
        MPI_Send(y_loc_prev.data() + loc_size - 2 * size, size * send2, MPI_DOUBLE, dest, 57, MPI_COMM_WORLD);
        MPI_Recv(y_loc_prev.data(), size * recv2, MPI_DOUBLE, source, 57, MPI_COMM_WORLD, &statU);

        MPI_Send(y_loc_prev.data() + size, size * recv2, MPI_DOUBLE, source, 65, MPI_COMM_WORLD);
        MPI_Recv(y_loc_prev.data() + loc_size - size, size * send2, MPI_DOUBLE, dest, 65, MPI_COMM_WORLD, &statL);
        MPI_Send(y_loc_prev.data() + size, size * recv1, MPI_DOUBLE, source, 67, MPI_COMM_WORLD);
        MPI_Recv(y_loc_prev.data() + loc_size - size, size * send1, MPI_DOUBLE, dest, 67, MPI_COMM_WORLD, &statL);

        for (int i = 1; i < loc_size / size - 1; ++i)
            for (int j = 1; j < size - 1; ++j)
            {
                //y_loc[i * size + j] = (h * h * vec_right_loc[i * size + j] \
                    //+ y_loc_prev[(i - 1) * size + j] + y_loc_prev[(i + 1) * size + j] \
                    //+ y_loc_prev[i * size + (j - 1)] + y_loc_prev[i * size + (j + 1)]) / (k_square + 4.0);
                y_loc[i * size + j] = c * (h * h * initialConditions.func_right((i + offset) * h, j * h, kSquare) + y_loc_prev[(i - 1) * size + j] + y_loc_prev[(i + 1) * size + j] + y_loc_prev[i * size + (j - 1)] + y_loc_prev[i * size + (j + 1)]);
            }

        norm = infiniteNorm(y_loc, y_loc_prev, recv_disp, loc_size - size);

        MPI_Allreduce(&norm, &finalNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    } while (finalNorm > eps);

    if (processId == 0)
        timeEnd = MPI_Wtime();

    MPI_Gatherv(y_loc.data() + recv_disp, loc_size - size - extr_size, MPI_DOUBLE, y.data(), len.data(), disp.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (processId == 0)
    {
        std::cout << "jacobi, Send Recv:\n";
        std::cout << "\tstop_crit = " << finalNorm << "\n\tNumber of iterations: " << iterationsNumber << "\n";
        //for (int i = 0; i < size; ++i) {
        //    for (int j = 0; j < size; ++j) {
        //        std::cout << y[i * size + j] << "\t";
        //    }
        //    std::cout << "\n";
        //}
    }
    time = timeEnd - timeStart;
}
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

double norm_inf(std::vector<double> &vec1, std::vector<double> &vec2, int begin, int end) {
    double norm = 0.0;
    for (int i = begin; i < end; ++i) {
        double elem = fabs(vec1[i] - vec2[i]);
        if (norm < elem) {
            {
                norm = elem;
            }
        }
    }
    return norm;
}

void initVecLoc(std::vector<double> &y, std::vector<double> &vec_right,
    double h, int size, double k_square, int np, int myid, std::vector<double> &y_loc,
    std::vector<double> &y_loc_prev, std::vector<double> &vec_right_loc, std::vector<int> &len,
    std::vector<int> &disp, int &loc_size, int &recv_disp, int &extr_size, int &offset) {
    
    std::vector<int> vec_offset(np);
    
    if (myid == 0) {
        len.resize(np);
        disp.resize(np);
        for (int i = 0; i < np; ++i) {
            len[i] = (size / np) * size;
        }
        for (int i = 0; i < size - (size / np) * (np); ++i) {
            len[i] += size;
        }

        disp[0] = 0;
        vec_offset[0] = -1;
        //if (np > 1)
        for (int i = 0; i < np - 1; ++i) {
            disp[i + 1] = disp[i] + len[i];
            vec_offset[i + 1] = vec_offset[i] + len[i] / size;
        }
        vec_offset[0] = 0;
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatter(len.data(), 1, MPI_INT, &loc_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(vec_offset.data(), 1, MPI_INT, &offset, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (myid == 0) {
        y.resize(size * size);
        //vec_right.resize(size * size);
        comp_vec_zero(y);
        //comp_vec_right(vec_right, size, h, k_square);
    }

    if (np > 1) {
        if ((myid == 0) || (myid == np - 1)) {
            loc_size += size;
            extr_size = 0.0;
        } else {
            loc_size += 2 * size;
            extr_size = size;
        }
    } else {
        extr_size = -size;
    }

    y_loc.resize(loc_size);
    y_loc_prev.resize(loc_size);
    vec_right_loc.resize(loc_size);
    
    recv_disp = (myid == 0) ? 0.0 : size;
    
    MPI_Scatterv(y.data(), len.data(), disp.data(), MPI_DOUBLE, y_loc.data() + recv_disp, loc_size - size - extr_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Scatterv(vec_right.data() + size, len.data(), disp.data(), MPI_DOUBLE, vec_right_loc.data() + size, loc_size - 2 * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void comp_vec_zero(std::vector<double> &y) {
    for (int i = 0; i < y.size(); ++i) {
        y[i] = 0.0;
    }
}

void send_recv_scheme(const int np, const int myid, int &dest, int &source, int &send1, int &recv1, \
    int &send2, int &recv2) {
    source_dest(np, myid, dest, source);
    if (myid % 2 == 0) {
        send1 = 1;
        recv1 = 0;
    } else {
        send1 = 0;
        recv1 = 1;
    }
    if (myid == np - 1) {
        send1 = 0;
    }
    recv2 = send1;
    if (myid == 0) {
        recv2 = 0;
    }
    if ((myid == np - 1) && (np % 2 != 0)) {
        recv2 = 1;
    }
    send2 = recv1;
    if ((myid == np - 1) && (np % 2 == 0)) {
        send2 = 0;
    }
    //fprintf(stdout, "myid = %d, source = %d, dest = %d, send1 = %d, recv1 = %d, send2 = %d, recv2 = %d", myid, source, dest, send1, recv1, send2, recv2);
    //fflush(stdout);
}

void source_dest(const int np, const int myid, int &dest, int &source) {
    dest = myid + 1;
    source = myid - 1;
    if (myid == 0) {
        source = np - 1;
    } else if (myid == np - 1) {
        dest = 0;
    }
}
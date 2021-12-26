#include <vector>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <random>
#include "mpi.h"
#include <omp.h>

using namespace std;

const double G = 6.67E-11;

ostream &operator<<(ostream &stream, const vector<double> &vec)
{
    //stream << "[";
    for (auto it = vec.begin(); it != --vec.end(); ++it)
    {
        stream << setw(20) << *it; //<< ", ";
    }
    stream << setw(20) << vec.back();
    //stream << "]";
    return stream;
}

vector<double> operator*(const vector<double> &lhs, const double &rhs)
{
    vector<double> result(lhs.size());
    for (int i = 0; i < lhs.size(); ++i)
        result[i] = lhs[i] * rhs;
    return result;
}

vector<double> operator*(const double &lhs, const vector<double> &rhs)
{
    return rhs * lhs;
}

vector<double> operator+(const vector<double> &lhs, const vector<double> &rhs)
{
    if (lhs.size() != rhs.size())
        throw invalid_argument("vectors have different sizes");
    vector<double> result(lhs.size());
    for (int i = 0; i < lhs.size(); ++i)
        result[i] = lhs[i] + rhs[i];
    return result;
}

vector<double> &operator+=(vector<double> &lhs, const vector<double> &rhs)
{
    if (lhs.size() != rhs.size())
        throw invalid_argument("vectors have different sizes");
    for (int i = 0; i < lhs.size(); ++i)
        lhs[i] += rhs[i];
    return lhs;
}

double distance3(double x1, double y1, double z1, double x2, double y2, double z2)
{
    // return pow(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2), 3.0 / 2);
    return pow(pow((x2 - x1), 2) + pow((y2 - y1), 2) + pow((z2 - z1), 2), 1.5);
}

vector<double> func(vector<double> &posvel, const vector<double> &masses)
{
    double eps = 1E-15;
    size_t n = masses.size();
    vector<double> result(6 * n);
    for (int i = 0; i < n; ++i)
    {
        result[3 * i] = posvel[3 * i + 3 * n]; // dr/dt=v
        result[3 * i + 1] = posvel[3 * i + 3 * n + 1];
        result[3 * i + 2] = posvel[3 * i + 3 * n + 2];
        double ax(0), ay(0), az(0);
        for (int j = 0; j < n; ++j)
        {
            double coef = -G * masses[j] / max(distance3(posvel[3 * i], posvel[3 * i + 1], posvel[3 * i + 2], posvel[3 * j], posvel[3 * j + 1], posvel[3 * j + 2]), eps);
            ax += (posvel[3 * i] - posvel[3 * j]) * coef;
            ay += (posvel[3 * i + 1] - posvel[3 * j + 1]) * coef;
            az += (posvel[3 * i + 2] - posvel[3 * j + 2]) * coef;
        }
        result[3 * i + 3 * n] = ax; //dv/dt=a
        result[3 * i + 3 * n + 1] = ay;
        result[3 * i + 3 * n + 2] = az;
    }
    return result;
}

vector<double> func_mpi(vector<double> &posvel, const vector<double> &masses, int np, int rank)
{
    double eps = 1E-15;
    size_t n = masses.size();
    vector<double> result(6 * n);
    vector<double> resultAll;
    if (rank == 0)
        resultAll.resize(6 * n);
#pragma omp parallel for
    for (int i = rank * n / np; i < (rank + 1) * n / np; ++i)
    {
        result[3 * i] = posvel[3 * i + 3 * n]; // dr/dt=v
        result[3 * i + 1] = posvel[3 * i + 3 * n + 1];
        result[3 * i + 2] = posvel[3 * i + 3 * n + 2];
        double ax(0), ay(0), az(0);
        for (int j = 0; j < n; ++j)
        {
            double coef = -G * masses[j] / max(distance3(posvel[3 * i], posvel[3 * i + 1], posvel[3 * i + 2], posvel[3 * j], posvel[3 * j + 1], posvel[3 * j + 2]), eps);
            ax += (posvel[3 * i] - posvel[3 * j]) * coef;
            ay += (posvel[3 * i + 1] - posvel[3 * j + 1]) * coef;
            az += (posvel[3 * i + 2] - posvel[3 * j + 2]) * coef;
        }
        result[3 * i + 3 * n] = ax; //dv/dt=a
        result[3 * i + 3 * n + 1] = ay;
        result[3 * i + 3 * n + 2] = az;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //  MPI_Allgather(result.data() + 3 * rank * n / np, 3 * n / np, MPI_DOUBLE, resultAll.data() + 3 * rank * n / np , 3 * n / np, MPI_DOUBLE, MPI_COMM_WORLD);
    //  printf("rank: %d \nid: %d\n \nnum: %d\nsize: %d",rank, 3 * rank * n / np + 3 * n, 3 * n / np, result.size());
    //  cout << rank << " ????????????\n" << endl;
    //  MPI_Allgather(result.data() + 3 * rank * n / np + 3 * n, 3 * n / np, MPI_DOUBLE, resultAll.data() + 3 * rank * n / np + 3 * n, 3 * n / np, MPI_DOUBLE, MPI_COMM_WORLD);
    //  cout << rank << " watafak\n" << endl;

    MPI_Gather(result.data() + 3 * rank * n / np, 3 * n / np, MPI_DOUBLE, resultAll.data() + 3 * rank * n / np, 3 * n / np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Gather(result.data() + 3 * rank * n / np + 3 * n, 3 * n / np, MPI_DOUBLE, resultAll.data() + 3 * rank * n / np + 3 * n, 3 * n / np, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank != 0)
        resultAll.swap(result);
    MPI_Bcast(resultAll.data(), 6 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return resultAll;
}

void rungekutta2(const vector<double> &init, const vector<double> &masses, double tmax, double tau, int k, const string &filename)
{
    ofstream output(filename);
    double tcurrent = 0;
    vector<double> pcurrent(init);
    vector<double> k1(init), k2(init);
    output << setprecision(10) << tcurrent << " " << pcurrent << '\n';
    int counter = 0;
    while (tcurrent < tmax)
    {
        k1 = func(pcurrent, masses);
        auto temp = pcurrent + tau / 2 * k1;
        k2 = func(temp, masses);
        pcurrent += tau * k2;
        tcurrent += tau;
        if (++counter % k == 0)
            output << setprecision(10) << tcurrent << " " << pcurrent << '\n';
        //cout << tcurrent << endl;
    }
}

void rungekutta2_mpiopenmp(const vector<double> &init, const vector<double> &masses, double tmax, double tau, int k, const string &filename, int np, int rank)
{
    ofstream output(filename);
    double tcurrent = 0;
    vector<double> pcurrent(init);
    vector<double> k1(init), k2(init);
    if (rank == 0)
        output << setprecision(10) << tcurrent << " " << pcurrent << '\n';
    int counter = 0;
    while (tcurrent < tmax)
    {
        k1 = func_mpi(pcurrent, masses, np, rank);
        auto temp = pcurrent + tau / 2 * k1;
        k2 = func_mpi(temp, masses, np, rank);
        pcurrent += tau * k2;
        tcurrent += tau;
        if (++counter % k == 0 && rank == 0)
            output << setprecision(10) << tcurrent << " " << pcurrent << '\n';
        //cout << tcurrent << endl;
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int np, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    vector<double> init = {1.0, 0.0, 0.0, 0.0, 0.9, 0.0,
                           0.0, 1.0, 0.0, -0.9, 0.0, 0.0,
                           -1.0, 0.0, 0.0, 0.0, -0.9, 0.0,
                           0.0, -1.0, 0.0, 0.9, 0.0, 0.0};
    vector<double> init2(init.size());
    for (int i = 0; i < init.size() / 6; ++i)
    {
        init2[3 * i] = init[6 * i];
        init2[3 * i + 1] = init[6 * i + 1];
        init2[3 * i + 2] = init[6 * i + 2];
        init2[3 * i + init.size() / 2] = init[6 * i + 3];
        init2[3 * i + init.size() / 2 + 1] = init[6 * i + 4];
        init2[3 * i + init.size() / 2 + 2] = init[6 * i + 5];
    }
    vector<double> masses = {8810324116.227, 8810324116.227, 8810324116.227, 8810324116.227};
    //rungekutta2(init2, masses, 20, 0.025, 4, "outex.txt");
    //rungekutta2_mpiopenmp(init2, masses, 20, 0.05, 2, "outMPI.txt",np,rank);

    // std::random_device seeder;
    // const auto seed = seeder.entropy() ? seeder() : time(nullptr);
    // std::mt19937 eng(static_cast<std::mt19937::result_type>(seed));
    // std::uniform_real_distribution<double> distr1(-10, 10);
    // std::uniform_real_distribution<double> distr2(-1, 1);
    // std::uniform_real_distribution<double> distr3(0.5, 2.0);
    // auto rndpos = std::bind(distr1, eng);
    // auto rndvel = std::bind(distr2, eng);
    // auto rndmas = std::bind(distr2, eng);

    int N = 4;

    // init.resize(6 * N);

    // if (rank == 0)
    // {
    //     for (int i = 0; i < 3 * N; ++i)
    //         init[i] = rndpos();
    //     for (int i = 3 * N; i < 6 * N; ++i)
    //         init[i] = rndvel();
    // }
    MPI_Bcast(init.data(), 6 * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    masses.resize(N);

    auto t1 = chrono::high_resolution_clock::now();
    if (rank == 0)
        rungekutta2(init, masses, 0.1, 0.05, 1, "out.txt");
    auto t2 = chrono::high_resolution_clock::now();
    if (rank == 0)
        cout << "Elapsed seq time = " << chrono::duration_cast<chrono::seconds>(t2 - t1).count() << " seconds \n";

    t1 = chrono::high_resolution_clock::now();
    rungekutta2_mpiopenmp(init, masses, 0.1, 0.05, 1, "out3.txt", np, rank);
    t2 = chrono::high_resolution_clock::now();
    if (rank == 0)
        cout << "Elapsed mpi+openmp time = " << chrono::duration_cast<chrono::seconds>(t2 - t1).count() << " seconds \n";

    MPI_Finalize();

    return 0;
}

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <vector>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <random>

using namespace std;

__device__ const double G = 6.67E-11;

ostream &operator<<(ostream &stream, const vector<double> &vec)
{
    for (auto it = vec.begin(); it != --vec.end(); ++it)
    {
        stream << setw(20) << *it;
    }
    stream << setw(20) << vec.back();

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
    vector<double> result(lhs.size());
    for (int i = 0; i < lhs.size(); ++i)
        result[i] = lhs[i] + rhs[i];
    return result;
}

vector<double> &operator+=(vector<double> &lhs, const vector<double> &rhs)
{
    for (int i = 0; i < lhs.size(); ++i)
        lhs[i] += rhs[i];
    return lhs;
}

__host__ __device__ double distance3(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return pow(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2), 3.0 / 2);
}

vector<double> func(vector<double> &posvel, const vector<double> &masses)
{
    double eps = 1E-15;
    size_t n = masses.size();
    vector<double> result(6 * n);
    for (int i = 0; i < n; ++i)
    {
        result[3 * i] = posvel[3 * i + 3 * n];
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
        result[3 * i + 3 * n] = ax;
        result[3 * i + 3 * n + 1] = ay;
        result[3 * i + 3 * n + 2] = az;
    }
    return result;
}

__global__ void funcpar(double *result, double *posvel, double *masses, int n)
{
    double eps = 1E-15;
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < n)
    {
        result[3 * idx] = posvel[3 * idx + 3 * n];
        result[3 * idx + 1] = posvel[3 * idx + 3 * n + 1];
        result[3 * idx + 2] = posvel[3 * idx + 3 * n + 2];
        double ax = 0, ay = 0, az = 0;
        for (int j = 0; j < n; ++j)
        {
            double coef = -G * masses[j] / max(distance3(posvel[3 * idx], posvel[3 * idx + 1], posvel[3 * idx + 2], posvel[3 * j], posvel[3 * j + 1], posvel[3 * j + 2]), eps);
            ax += (posvel[3 * idx] - posvel[3 * j]) * coef;
            ay += (posvel[3 * idx + 1] - posvel[3 * j + 1]) * coef;
            az += (posvel[3 * idx + 2] - posvel[3 * j + 2]) * coef;
        }
        result[3 * idx + 3 * n] = ax;
        result[3 * idx + 3 * n + 1] = ay;
        result[3 * idx + 3 * n + 2] = az;
    }
}

__global__ void plustimesequal(double *result, double *rhs, double tau, int n)
{
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < n)
    {
        result[idx] += tau * rhs[idx];
    }
}

__global__ void plustimes(double *result, double *lhs, double *rhs, double tau, int n)
{
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < n)
    {
        result[idx] = lhs[idx] + tau * rhs[idx];
    }
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
    }
}

void rungekutta2par(const vector<double> &init, const vector<double> &masses, double tmax, double tau, int k, const string &filename)
{
    ofstream output(filename);
    double tcurrent = 0;
    vector<double> pcurrent(init);
    vector<double> k1(init), k2(init);
    output << setprecision(10) << tcurrent << " " << pcurrent << '\n';

    int counter = 0;
    int size = init.size();
    int msize = masses.size();
    double *pcurrent_dev;
    double *k1_dev;
    double *k2_dev;
    double *masses_dev;
    double *temp_dev;

    cudaMalloc(&pcurrent_dev, size * sizeof(double));
    cudaMalloc(&k1_dev, size * sizeof(double));
    cudaMalloc(&k2_dev, size * sizeof(double));
    cudaMalloc(&masses_dev, msize * sizeof(double));
    cudaMalloc(&temp_dev, size * sizeof(double));

    cudaMemcpy(masses_dev, masses.data(), msize * sizeof(double), cudaMemcpyHostToDevice);

    const int bs = 32;

    while (tcurrent < tmax)
    {
        dim3 blocks((msize + bs) / bs, 1, 1);
        dim3 threads(bs, 1, 1);

        cudaMemcpy(pcurrent_dev, pcurrent.data(), size * sizeof(double), cudaMemcpyHostToDevice);
        funcpar<<<blocks, threads>>>(k1_dev, pcurrent_dev, masses_dev, msize);// calculate k1

        cudaMemcpy(k1.data(), k1_dev, size * sizeof(double), cudaMemcpyDeviceToHost);
        plustimes<<<(size + bs) / bs, bs>>>(temp_dev, pcurrent_dev, k1_dev, tau / 2, size);// step with k1 on tau/2

        funcpar<<<blocks, threads>>>(k2_dev, temp_dev, masses_dev, msize);//calculate k2
        cudaMemcpy(k2.data(), k2_dev, size * sizeof(double), cudaMemcpyDeviceToHost);

        cudaMemcpy(pcurrent_dev, pcurrent.data(), size * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(k2_dev, k2.data(), size * sizeof(double), cudaMemcpyHostToDevice);
        plustimesequal<<<(size + bs) / bs, bs>>>(pcurrent_dev, k2_dev, tau, size);//finally, steps

        cudaMemcpy(pcurrent.data(), pcurrent_dev, size * sizeof(double), cudaMemcpyDeviceToHost);

        tcurrent += tau;
        if (++counter % k == 0)
            output << setprecision(10) << tcurrent << " " << pcurrent << '\n';
    }
    
    cudaFree(pcurrent_dev);
    cudaFree(k1_dev);
    cudaFree(k2_dev);
    cudaFree(masses_dev);
    cudaFree(temp_dev);
}

int main(int argc, char **argv)
{
    vector<double> init;
    vector<double> masses;
    
    size_t N = 5000;
    bool seq = False;
 
    init.resize(6 * N);
    masses.resize(N);

    std::random_device seeder;
    const auto seed = seeder.entropy() ? seeder() : time(nullptr);
    std::mt19937 eng(static_cast<std::mt19937::result_type>(seed));
    std::uniform_real_distribution<double> distr1(-10, 10);
    std::uniform_real_distribution<double> distr2(-1, 1);
    std::uniform_real_distribution<double> distr3(0.5, 2.0);
    auto rndpos = std::bind(distr1, eng);
    auto rndvel = std::bind(distr2, eng);
    auto rndmas = std::bind(distr2, eng);

    for (int i = 0; i < 3 * N; ++i)
        init[i] = rndpos();
    for (int i = 3 * N; i < 6 * N; ++i)
        init[i] = rndvel();
    for (int i = 0; i < N; ++i)
        masses[i] = rndmas();

    auto t1 = chrono::high_resolution_clock::now();
    if(seq)
        rungekutta2(init, masses, 0.1, 0.05, 1, "out2");
    auto t2 = chrono::high_resolution_clock::now();

    auto seqtime = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
    cout << "Seq time " << seqtime << " ms\n";

    t1 = chrono::high_resolution_clock::now();
    rungekutta2par(init, masses, 0.1, 0.05, 1, "out3");
    t2 = chrono::high_resolution_clock::now();

    auto partime = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
    cout << "Cuda time(?) " << partime << " ms\n";

    cout << "acceleration: " << seqtime * 1.0 / partime << "\n";

    return 0;
}

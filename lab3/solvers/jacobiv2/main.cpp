#include "iostream"
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "mpi.h"
#include "../../common/Common.hpp"
#include "../../init_conds/InitialConditions.hpp"
const double PI = 3.14159265358979323846;
const double EPS = 1.e-8;

void divideVectorBetweenProcesses(int size, int &locationSize, std::vector<double> &y, std::vector<double> &yPart,
                                  std::vector<double> &yPreviousPart, std::vector<double> &f_loc,
                                  double h, double k, int numberOfProcesses, int processId,
                                  std::vector<int> &len, std::vector<int> &displace, std::vector<int> &i_glob,
                                  int &locationShift, int &extraLocationSize, InitialConditions initialConditions);

void jacobiV3(int size, std::vector<double> &y,
              double h, double k, double coef1,
              double coef2, int numberOfProcesses, int processId,
              InitialConditions initialConditions);

//void filing_sol(std::string filename, int N, double h, std::vector<double> &y);
double normV1(const std::vector<double> &vec1, const std::vector<double> &vec2);
double normV2(const std::vector<double> &vec1, const std::vector<double> &vec2, int itbegin, int itend);

double rightPart(double x, double y, double k)
{
        return 2.0 * sin(PI * y) + k * k * (1.0 - x) * x * sin(PI * y) + PI * PI * (1.0 - x) * x * sin(PI * y);
}

double analyticSolution(double x, double y)
{
        return x * (1 - x) * sin(PI * y);
}

void fillYWithZeros(std::vector<double> &y)
{
        for (int i = 0; i < y.size(); ++i)
                y[i] = 0.0;
}

void computeRealSolution(const int size, std::vector<double> &u, double h)
{
        u.resize(size * size);
        for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                        u[i * size + j] = analyticSolution(j * h, i * h);
}

int main(int argc, char **argv)
{
        int processId, numberOfProcesses;
        int namelen;
        char processor_name[MPI_MAX_PROCESSOR_NAME];
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);
        MPI_Comm_rank(MPI_COMM_WORLD, &processId);
        MPI_Get_processor_name(processor_name, &namelen);

        InitialConditions initialConditions = InitialConditions();
        int size = initialConditions.n + 1;
        double coef1, coef2, k;
        k = 1 / initialConditions.h;
        coef1 = 1.0 / (4.0 + initialConditions.h * initialConditions.h * k * k);
        coef2 = initialConditions.h * initialConditions.h * coef1;

        std::vector<double> y, u;
        if (processId == 0)
                computeRealSolution(size, u, initialConditions.h);

        MPI_Barrier(MPI_COMM_WORLD);
        jacobiV3(size, y, initialConditions.h, k, coef1, coef2, numberOfProcesses, processId, initialConditions);
        if (processId == 0)
                fprintf(stdout, "Difference: %e \n", normV1(u, y));

        MPI_Finalize();
        return 0;
}

void divideVectorBetweenProcesses(int size, int &locationSize, std::vector<double> &y, std::vector<double> &yPart,
                                  std::vector<double> &yPreviousPart, std::vector<double> &f_loc,
                                  double h, double k, int numberOfProcesses, int processId,
                                  std::vector<int> &len, std::vector<int> &displace, std::vector<int> &i_glob,
                                  int &locationShift, int &extraLocationSize, InitialConditions initialConditions)
{

        i_glob.resize(numberOfProcesses);

        if (processId == 0)
        {
                len.resize(numberOfProcesses);
                displace.resize(numberOfProcesses);
                for (int i = 0; i < numberOfProcesses; ++i)
                        len[i] = size / numberOfProcesses * size;

                for (int i = 0; i < size - (size / numberOfProcesses) * numberOfProcesses; ++i)
                        len[i] += size;

                displace[0] = 0;
                for (int i = 0; i < numberOfProcesses - 1; ++i)
                        displace[i + 1] = displace[i] + len[i];

                i_glob[0] = 0;
                for (int i = 0; i < numberOfProcesses - 1; ++i)
                        i_glob[i + 1] = i_glob[i] + len[i] / size;
        }

        MPI_Barrier;
        MPI_Scatter(len.data(), 1, MPI_INT, &locationSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(i_glob.data(), i_glob.size(), MPI_INT, 0, MPI_COMM_WORLD);

        f_loc.resize(locationSize);
        for (int i = 0; i < f_loc.size() / size; i++)
                for (int j = 0; j < size; j++)
                        f_loc[j + i * size] = rightPart(j * h, (i + i_glob[processId]) * h, k);

        extraLocationSize = 0;
        if (numberOfProcesses > 1)
        {
                if (processId == 0 or processId == numberOfProcesses - 1)
                {
                        extraLocationSize = size;
                        locationSize += extraLocationSize;
                }
                else
                {
                        extraLocationSize = 2 * size;
                        locationSize += extraLocationSize;
                }
        }

        locationShift = (processId == 0) ? 0 : 1;

        if (processId == 0)
        {
                y.resize(size * size);
                fillYWithZeros(y);
        }

        yPart.resize(locationSize);
        yPreviousPart.resize(locationSize);
        MPI_Scatterv(y.data(),
                     len.data(),
                     displace.data(),
                     MPI_DOUBLE,
                     yPart.data() + locationShift * size,
                     locationSize - extraLocationSize,
                     MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void jacobiV3(int size, std::vector<double> &y, double h, double k,
              double coef1, double coef2,
              int numberOfProcesses, int processId, InitialConditions initialConditions)
{
        int locationSize, locationShift, extraLocationSize;
        std::vector<double> yPart, yPreviousPart, f_loc;
        std::vector<int> len, displace, i_glob;

        divideVectorBetweenProcesses(size, locationSize, y,
                                     yPart, yPreviousPart, f_loc, h, k,
                                     numberOfProcesses, processId, len, displace, i_glob,
                                     locationShift, extraLocationSize, initialConditions);

        double norm, finalNorm;
        int iterationsNumber = 0;

        std::vector<MPI_Request> reqB, reqA;
        if ((processId == 0) || (processId == numberOfProcesses - 1))
        {
                reqB.resize(1);
                reqA.resize(1);
        }
        else
        {
                reqB.resize(2);
                reqA.resize(2);
        }
        int nB = 0, nA = 0;
        double timeStart, timeEnd;

        if (processId != numberOfProcesses - 1)
        {
                MPI_Send_init(yPreviousPart.data() + locationSize - 2 * size, size, MPI_DOUBLE, processId + 1, 42, MPI_COMM_WORLD, reqA.data() + nA);
                MPI_Recv_init(yPreviousPart.data() + locationSize - size, size, MPI_DOUBLE, processId + 1, 41, MPI_COMM_WORLD, reqB.data() + nB);
                nA++;
                nB++;
        }
        if (processId != 0)
        {
                MPI_Recv_init(yPreviousPart.data(), size, MPI_DOUBLE, processId - 1, 42, MPI_COMM_WORLD, reqA.data() + nA);
                MPI_Send_init(yPreviousPart.data() + size, size, MPI_DOUBLE, processId - 1, 41, MPI_COMM_WORLD, reqB.data() + nB);
                nA++;
                nB++;
        }

        std::vector<MPI_Status> statB(nB), statA(nA);

        if (processId == 0)
                timeStart = MPI_Wtime();
        do
        {
                iterationsNumber++;

                memcpy(yPreviousPart.data(), yPart.data(), locationSize * sizeof(double));

                MPI_Startall(nA, reqA.data());
                MPI_Startall(nB, reqB.data());

                for (int i = 2; i < locationSize / size - 2; i++)
                        for (int j = 1; j < size - 1; j++)
                                yPart[i * size + j] = coef1 * (yPreviousPart[(i + 1) * size + j] +
                                                               yPreviousPart[(i - 1) * size + j] +
                                                               yPreviousPart[i * size + (j + 1)] +
                                                               yPreviousPart[i * size + (j - 1)]) +
                                                      coef2 * f_loc[j + (i - locationShift) * size];

                MPI_Waitall(nB, reqB.data(), statB.data());
                MPI_Waitall(nA, reqA.data(), statA.data());

                int i = 1;
                for (int j = 1; j < size - 1; j++)
                        yPart[i * size + j] = coef1 * (yPreviousPart[(i + 1) * size + j] +
                                                       yPreviousPart[(i - 1) * size + j] +
                                                       yPreviousPart[i * size + (j + 1)] +
                                                       yPreviousPart[i * size + (j - 1)]) +
                                              coef2 * f_loc[j + (i - locationShift) * size];

                i = locationSize / size - 2;
                for (int j = 1; j < size - 1; j++)
                        yPart[i * size + j] = coef1 * (yPreviousPart[(i + 1) * size + j] +
                                                       yPreviousPart[(i - 1) * size + j] +
                                                       yPreviousPart[i * size + (j + 1)] +
                                                       yPreviousPart[i * size + (j - 1)]) +
                                              coef2 * f_loc[j + (i - locationShift) * size];

                norm = normV2(yPart, yPreviousPart, locationShift * size, locationSize - extraLocationSize);
                MPI_Allreduce(&norm, &finalNorm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        } while (finalNorm > EPS);
        if (processId == 0)
                timeEnd = MPI_Wtime();

        MPI_Gatherv(yPart.data() + locationShift * size, locationSize - extraLocationSize, MPI_DOUBLE,
                    y.data(), len.data(), displace.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (processId == 0)
        {
                std::cout << "\033[1;32mJacobi. MPI_Send_init. MPI_Recv_init\033[0m" << std::endl;
                std::cout << "Substraction norm: " << finalNorm << std::endl;
                std::cout << "Number of iterationsNumber: " << iterationsNumber << std::endl;
                std::cout << "Time: " << timeEnd - timeStart << std::endl;
        }
}

double normV1(const std::vector<double> &vec1, const std::vector<double> &vec2)
{
        double res;
        double norm = 0.0;
        for (int i = 0; i < vec1.size(); i++)
        {
                res = fabs(vec1[i] - vec2[i]);
                if (res > norm)
                        norm = res;
        }
        return (norm);
}

double normV2(const std::vector<double> &vec1, const std::vector<double> &vec2, int it_begin, int it_end)
{
        double res;
        double norm = 0.0;
        for (int i = it_begin; i < it_end; i++)
        {
                res = fabs(vec1[i] - vec2[i]);
                if (res > norm)
                        norm = res;
        }
        return (norm);
}

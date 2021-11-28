#include <iostream>
#include <fstream>
// #include "ofstream"
#include <cmath>
#include <vector> //maybe array, but size may be large
#include <string>
#include "cassert"


const double EPS = 1e-3;

class Solver
{
public:
    Solver(size_t n_, std::string method);
    ~Solver();

    void solve();

    size_t n; //matrix n*n, vectors n
    // Ax=b
    std::vector<std::vector<double>> A; //left part matrix
    std::vector<double> x; //vector of unknowns
    std::vector<double> prev_x; //vector of unknowns(for stop function)
    std::vector<double> b; //right part vector

    std::string method; // Z - Zeidel, J - Jacobi, M - Red-Black mixed algorithm

// private:
    void step();
    void init();
};

Solver::Solver(size_t n_=8, std::string method_="J")
{
    n = n_;
    method = method_;

    prev_x = std::vector<double>(n, 0);
    x = std::vector<double>(n, 0);
    b = std::vector<double>(n, 0);

    A = std::vector<std::vector<double>>(n, std::vector<double>(n, 0));
}

Solver::~Solver()
{
    //nothing?
}

void Solver::solve()
{
    for(int i = 0; i < 10; ++i)
    {
        step();
    }
}

void Solver::init()
{
    size_t n_ = 0;
    std::ifstream file; //read matrix
    file.open("matrix.txt");
    file >> n_;
    assert(n == n_);
    for (size_t i = 0; i < n; ++i)
    {
        for(size_t j = 0; j < n; ++j)
        {
            file >> A[i][j];
        }
    }
    file.close();

    file.open("vector.txt"); // read vector
    file >> n_;
    assert(n == n_);
    for(int i = 0; i < n; ++i)
    {
        file >> b[i];
    }
    file.close();
}

void Solver::step()
{
    if (method == "J")
    {
        for(size_t i = 0; i < n; ++i)
        {
            x[i] = b[i];
            for(size_t j = 0; j < n; ++j)
            {
                if((i != j) & (prev_x[j] > EPS))
                {
                    x[i] -= A[i][j] * prev_x[j];
                }
            }
            x[i] /= A[i][i];
        }
    }
    else if (method == "Z")
    {
        for(size_t i = 0; i < n; ++i)
        {
            x[i] = b[i];
            for(size_t j = 0; j < i; ++j)
            {
                x[i] -= A[i][j] * x[j];
            }
            for(size_t j = i + 1; j < n; ++j)
            {
                x[i] -= A[i][j] * prev_x[j];
            }
            x[i] /= A[i][i];
        }
    }
    // std::cout << "iter \n";
    // for(int i = 0; i < n; ++i)
    // {
    //     std::cout << x[i] << "\n";
    // }
    prev_x = x;
    // for(int i = 0; i < n; ++i)
    // {
    //     std::cout << x[i] << "\n";
    // }
    // std::cout << "\n";
    // for(int i = 0; i < n; ++i)
    // {
    //     std::cout << x[i] << "\n";
    // }
}

int main()
{
    size_t n = 8;
    Solver A = Solver(n, "Z");
    A.init();

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            std::cout << A.A[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    for(int i = 0; i < n; ++i)
    {
        std::cout << A.b[i] << " ";
    }
    std::cout << "\n";

    A.solve();
    for(int i = 0; i < n; ++i)
    {
        std::cout << A.x[i] << "\n";
    }
}

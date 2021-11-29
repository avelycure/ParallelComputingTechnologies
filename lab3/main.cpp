#include "header.hpp"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    cout << "Hello" << endl;
    MPI_Finalize();
}

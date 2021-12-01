all:
	mpiCC main/main.cpp \
	solvers/seidel/seidel.cpp \
	solvers/jacobi/jacobi.cpp \
	matrix_funcs/common.cpp \
	-o compiled/lab3.out
	mpirun -np 4 compiled/lab3.out
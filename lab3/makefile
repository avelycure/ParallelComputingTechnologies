all:
	mpiCC main.cpp \
	src/solvers/black_red_seidel.cpp \
	src/solvers/jacobi.cpp \
	src/solvers/seidel.cpp \
	src/matrix_funcs/common.cpp

	-o lab3.out
	mpirun -np 4 lab3.out
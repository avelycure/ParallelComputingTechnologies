all:
	mpiCC main/Main.cpp \
	solvers/seidel/Seidel.cpp \
	solvers/jacobi/Jacobi.cpp \
	matrix_funcs/Common.cpp \
	-o compiled/lab3.out
	mpirun -np 4 compiled/lab3.out
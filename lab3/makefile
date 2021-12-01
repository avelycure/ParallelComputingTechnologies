all:
	mpiCC main/Main.cpp \
	solvers/seidel/Seidel.cpp \
	solvers/jacobi/Jacobi.cpp \
	common/Common.cpp \
	-o compiled/lab3.out
	mpirun -np 3 compiled/lab3.out
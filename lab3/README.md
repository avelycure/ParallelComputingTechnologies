# Laboratory work 3. Helmholtz equation.

## Lab assignment

Develop a program for solving the two-dimensional Helmholtz equation

<p>
  <img src="https://github.com/avelycure/avelycure/blob/master/assets/parallel_technologies/helmholtz/helmholts.png" width="300" />
</p>

in a square domain (x,y) in [0,1]x[0,1] with border conditions: u(x,0)=u(x,1)=u(0,y)=u(1,y)=0, and right part

<p>
  <img src="https://github.com/avelycure/avelycure/blob/master/assets/parallel_technologies/helmholtz/right_part.png" width="300" />
</p>

To numerically solve the equation on a rectangular uniform grid use the finite-difference scheme "cross" of the second order of approximation in
both independent variables.

Solve the resulting system of linear algebraic equations by iterative by the Jacobi and Seidel method.

## Features

The following library functions were used:
* MPI_Send. MPI_Recv
* MPI_Sendrecv
* MPI_Send_init. MPI_Recv_init

## Results

Number of nodes in each direction: 8193
Accuracy: 1e-8

### Jacobi. MPI_Send. MPI_Recv
|Nodes\Cores|1|2|3|4|
|:----------:|:----------:|:----------:|:----------:|:----------:|
|1|35.37|0|0|41.23|
|2|61.87|0|0|16.58|
|3|0|0|0|11.06|

### Jacobi. MPI_Sendrecv
|Nodes\Cores|1|2|3|4|
|:----------:|:----------:|:----------:|:----------:|:----------:|
|1|38.93|0|0|41.44|
|2|62.24|0|0|18.10|
|3|0|0|0|12.64|

### Jacobi. MPI_Send_init. MPI_Recv_init
|Nodes\Cores|1|2|3|4|
|:----------:|:----------:|:----------:|:----------:|:----------:|
|1|41.27|0|0|41.26|
|2|62.17|0|0|19.43|
|3|0|0|0|14.52|

### Seidel. MPI_Send. MPI_Recv
|Nodes\Cores|1|2|3|4|
|:----------:|:----------:|:----------:|:----------:|:----------:|
|1|22.98|0|0|22.37|
|2|33.74|0|0|10.93|
|3|0|0|0|8.91|

### Seidel. MPI_Sendrecv
|Nodes\Cores|1|2|3|4|
|:----------:|:----------:|:----------:|:----------:|:----------:|
|1|23.41|0|0|22.38|
|2|33.71|0|0|10.73|
|3|0|0|0|7.51|

### Seidel. MPI_Send_init. MPI_Recv_init
|Nodes\Cores|1|2|3|4|
|:----------:|:----------:|:----------:|:----------:|:----------:|
|1|26.58|0|0|25.22|
|2|37.94|0|0|12.38|
|3|0|0|0|8.82|

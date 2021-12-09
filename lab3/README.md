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

Number of nodes in each direction: 128
Coefficient k: 0.0078125
Accuracy: 1e-8
Time: milliseconds

### Jacobi. MPI_Send. MPI_Recv
|Nodes\Cores|1|2|3|4|
|:----------:|:----------:|:----------:|:----------:|:----------:|
|1|0|0|0|0|
|2|46|32|31|26|
|3|0|0|0|0|

### Jacobi. MPI_Sendrecv
|Nodes\Cores|1|2|3|4|
|:----------:|:----------:|:----------:|:----------:|:----------:|
|1|0|0|0|0|
|2|48|41|34|37|
|3|0|0|0|0|

### Jacobi. MPI_Send_init. MPI_Recv_init
|Nodes\Cores|1|2|3|4|
|:----------:|:----------:|:----------:|:----------:|:----------:|
|1|0|0|0|0|
|2|34|26|31|21|
|3|0|0|0|0|

### Seidel. MPI_Send. MPI_Recv
|Nodes\Cores|1|2|3|4|
|:----------:|:----------:|:----------:|:----------:|:----------:|
|1|0|0|0|0|
|2|37|34|27|30|
|3|0|0|0|0|

### Seidel. MPI_Sendrecv
|Nodes\Cores|1|2|3|4|
|:----------:|:----------:|:----------:|:----------:|:----------:|
|1|0|0|0|0|
|2|42|37|34|38|
|3|0|0|0|0|

### Seidel. MPI_Send_init. MPI_Recv_init
|Nodes\Cores|1|2|3|4|
|:----------:|:----------:|:----------:|:----------:|:----------:|
|1|0|0|0|0|
|2|37|33|27|28|
|3|0|0|0|0|

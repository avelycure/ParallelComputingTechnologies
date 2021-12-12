# Binary tree
Creating binary tree with C++ using OMP library.

## Time

### Building tree

|Number of kernels|Time(ms)|Acceleration|
|:----------:|:----------:|:----------:|
|1|11187||
|2|6321|1.77x|
|4|3526|3.17x|

### Finding sum

|Number of kernels|Time(ms)|Acceleration|
|:----------:|:----------:|:----------:|
|1|4326||
|2|2279|1.90x|
|4|1152|3.76x|

### Overall

|Number of kernels|Time(ms)|Acceleration|
|:----------:|:----------:|:----------:|
|1|15513||
|2|8600|1.8x|
|4|4678|3.31x|

### Link to array feature

On [gwassel](https://github.com/gwassel) notebook:

### Building tree

|Number of kernels|Time(ms)|Acceleration|
|:----------:|:----------:|:----------:|
|1|16952||
|2|8833|1.9x|

More acceleration due to less time spent in the critical section

### Finding sum

|Number of kernels|Time(ms)|Acceleration|
|:----------:|:----------:|:----------:|
|1|5762||
|2|3062|1.8x|

This method is faster than the initial one due to the fact that memory is allocated all at once, but for a real task, where you will need the option to add an element to the tree (and delete),
this method will be less flexible.

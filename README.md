# ABCD Algorithm for Tridiagonal solvers

For mid-evaluation:
- [X] Get bayer matrix
- [X] Uniform Partitioning of matrix A (12/16)
- [X] Augment matrix A
- [X] Construct P
- [X] Construct S from Y and P
- [X] Make test samples from random inputs
- [X] Test the working


### Using cmake
```sh
mkdir build
cd build
cmake ..
make
./abcd
```


### To run files using armadillo

```sh
g++ <filename> -O2 -larmadillo
./a.out
```

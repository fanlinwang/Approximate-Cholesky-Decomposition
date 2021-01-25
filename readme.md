# Approximate Cholesky Decomposition

Project repository for [Advanced Systems Lab, ETH Zurich, Spring 2020](https://acl.inf.ethz.ch/teaching/fastcode/2020/)

In this project, we optimized the implementation of the [Approximate Cholesky Decomposition algorithm](https://arxiv.org/abs/1605.02353) from an engineering perspective. [The original implementation](https://danspielman.github.io/Laplacians.jl/latest/usingSolvers/) is in Julia.

The final report can be found [here](https://github.com/fanlinwang/Approximate-Cholesky-Decomposition/blob/master/report.pdf).

Team member: Daoye Wang, Fanlin Wang, Dexin Yang, Kaifeng Zhao.

## Julia Code

- [x] approxCholTypes/LLord: a node in the linked list (i.e. a column of the sparse matrix)
- [x] approxCholTypes/LLMatOrd: a sparse matrix
- [x] approxCholTypes/LDLinv: decomposition of the matrix
- [x] approxCholTypes/LLcol: output of get_ll_col method
- [x] approxChol/LDLinv(): constructor of LDLinv
- [x] approxChol/LLMatOrd(a::SparseMatrixCSC{Tval,Tind}): constructor of LLMatOrd. SparseMatrixCSC: a class of sparse matrix. For now: write some LLMatOrd manually. (maybe a simple 3 by 3 matrix)
- [x] approxChol/print_ll_col: print a column of LLMatOrd
- [x] approxChol/get_ll_col: take out a vector of LLcol from LLMatOrd
- [x] approxChol/compressCol: somehow sort the vector of LLcol *
- [x] approxChol/approxChol: input LLMatOrd, apply approximate cholosky decomposition and output LDLinv
- [x] approxChol/LDLsolver: solve Lx = b given the decomposition
- [x] approxChol/forward!: forward substitution, which is a step in LDLsolver
- [x] approxChol/backward!: backward substitution, which is a step in LDLsolver
- [x] approxChol/approxchol_lapGiven: an API over 

## Build and Test
- Generate build files for all:
``` bash
mkdir build && cd build
cmake ..
make
```
- Run:
``` bash
./src/main [vertice] [edges]
./test/test_LDLsolver
```

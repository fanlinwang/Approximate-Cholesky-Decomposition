# Approximate Cholesky Decomposition

## Julia Code

approxCholTypes/LLord: a node in the linked list (i.e. a column of the sparse matrix)
approxCholTypes/LLMatOrd: a sparse matrix
approxCholTypes/LDLinv: decomposition of the matrix
approxCholTypes/LLcol: output of get_ll_col method

approxChol/LDLinv(): constructor of LDLinv
approxChol/LLMatOrd(a::SparseMatrixCSC{Tval,Tind}): constructor of LLMatOrd 
SparseMatrixCSC: a class of sparse matrix
For now: write some LLMatOrd manually. (maybe a simple 3 by 3 matrix)
approxChol/print_ll_col: print a column of LLMatOrd
approxChol/get_ll_col: take out a vector of LLcol from LLMatOrd
approxChol/compressCol: somehow sort the vector of LLcol *
approxChol/approxChol: input LLMatOrd, apply approximate cholosky decomposition and output LDLinv
approxChol/LDLsolver: solve Lx = b given the decomposition
approxChol/forward!: forward substitution, which is a step in LDLsolver
approxChol/backward!: backward substitution, which is a step in LDLsolver
approxChol/approxchol_lapGiven: an API over 
pcg?!!! send email to ask

## 
- Generate build files:
``` bash
mkdir build && cd build
cmake ..
```
- Compile:
``` bash
make
```

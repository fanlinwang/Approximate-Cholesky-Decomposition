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

## Project Structure
``` text
.
├── CMakeLists.txt
├── src
│   └── main.cpp
│   ├── approxChol.hpp
│   └── approxChol.cpp
│   └── approxCholTypes.h
│   └── common.h
│   └── tsc_x86.h
│   └── CMakeLists.txt
└── tests
    ├── test_LDLsolver.cpp
    ├── test.h
    └── CMakeLists.txt
```
- To add a new test file `test_xxx.cpp` under the `test/`, append the following two lines to `test/CMakeLists.txt`: 
``` cmake
add_executable(test_xxx test/test_xxx.cpp ${PROJECT_SOURCES})   # Name of executable.
target_link_libraries(test_xxx PRIVATE approxChol)  # Link the executable to the library built from src/*.cpp
```

## Build and Test
- Generate build files for all:
``` bash
cd team008 && mkdir build && cd build
cmake ..
```
- Compile:
``` bash
make
```
- Run:
``` bash
./src/main
./test/test_LDLsolver
```

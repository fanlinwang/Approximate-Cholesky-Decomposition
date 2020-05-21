# Approximate Cholesky Decomposition

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
└── test
    ├── test_LDLsolver.cpp
    ├── test.h
    └── CMakeLists.txt
```
- To add a new test file `test_xxx.cpp` under the `test/`, append the following two lines to `test/CMakeLists.txt`: 
``` cmake
add_executable(test_xxx test_xxx.cpp)   # Name of executable.
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
./src/main [vertice] [edges]
./test/test_LDLsolver
```


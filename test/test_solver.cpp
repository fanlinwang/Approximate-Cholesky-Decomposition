#include <iostream>
#include "approxChol.hpp"
#include "test.h"
#include "matrix.hpp"

int main(int argc, char **argv){

    int v = atoi(argv[1]);
    int e = atoi(argv[2]);

    SparseMatrix A(v, e);

    LLMatOrd llmat = LLMatOrd(A);
    std::cout << "Created random sparse matrix A.\n";

    // std::cout << "Run approxChol once.\n"
    // approxChol(llmat);

    return 0;

}
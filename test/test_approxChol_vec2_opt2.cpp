#include <iostream>
#include "approxChol.hpp"
#include "test.h"
#include "matrix.hpp"


int main(int argc, char **argv){

    int v = atoi(argv[1]);
    int e = atoi(argv[2]);

    SparseMatrix A(v, e);
    LLMatOrd_vector2 llmat = LLMatOrd_vector2(A);

    approxChol_vector2_opt2(llmat);

    return 0;

}

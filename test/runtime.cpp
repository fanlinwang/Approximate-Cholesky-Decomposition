#include "approxChol.hpp"
#include "test.h"
#include "matrix.hpp"
#include <algorithm>
#include <iostream>

int main(int argc, char **argv){

    int v = atoi(argv[1]);
    int e = atoi(argv[2]);
    int r = atoi(argv[3]);

    SparseMatrix A(v, e);
    LLMatOrd llmat = LLMatOrd(A);
    std::cout << v << " " <<  e << " ";

    elapsed_time(approxChol, llmat, "Orig", r);
    elapsed_time(approxChol_opt, llmat, "approxChol_opt", r);

    LLMatOrd_vector2 llmat2 = LLMatOrd_vector2(A);
    elapsed_time(approxChol_vector2, llmat2, "approxChol_vector2", r);
    elapsed_time(approxChol_vector2_merge, llmat2, "approxChol_vector2_merge", r);
    elapsed_time(approxChol_vector2_opt2, llmat2, "approxChol_vector2_opt2", r);
    elapsed_time(approxChol_vector2_mergerand, llmat2, "approxChol_vector2_mergerand", r);
    std::cout << "\n";
    
    return 0;

}

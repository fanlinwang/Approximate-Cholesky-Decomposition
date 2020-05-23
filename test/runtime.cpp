#include "approxChol.hpp"
#include "test.h"
#include "matrix.hpp"
#include <algorithm>
#include <iostream>
#include <vector>

int main(int argc, char **argv){

    int v = atoi(argv[1]);
    int e = atoi(argv[2]);

    SparseMatrix A(v, e);
    LLMatOrd llmat = LLMatOrd(A);

    std::vector<LLMatOrd> llmats(2, llmat);
    elapsed_time(approxChol, llmats[0], "Orig");
    elapsed_time(approxChol_opt, llmats[1], "approxChol_opt");

    LLMatOrd_vector2 llmat2 = LLMatOrd_vector2(A);
    std::vector<LLMatOrd_vector2> llmats2(4, llmat2);
    elapsed_time(approxChol_vector2, llmats2[0], "approxChol_vector2");
    elapsed_time(approxChol_vector2_merge, llmats2[1], "approxChol_vector2_merge");
    elapsed_time(approxChol_vector2_opt2, llmats2[2], "approxChol_vector2_opt2");
    elapsed_time(approxChol_vector2_mergerand, llmats2[3], "approxChol_vector2_mergerand");
    
    return 0;

}

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
    std::cout << v << " " << e << " ";

    std::vector<LLMatOrd> llmats(4, llmat);
    elapsed_time(&approxChol, llmats[0], "baseline", r);
    elapsed_time(&approxChol_opt, llmats[1], "inline", r);
    elapsed_time(&approxChol_opt2, llmats[2], "inline SIMD", r);

    LLMatOrd_vector2 llmat2 = LLMatOrd_vector2(A);
    std::vector<LLMatOrd_vector2> llmats2(9, llmat2);
    elapsed_time(&approxChol_vector2, llmats2[0], "vector2",r);
    elapsed_time(&approxChol_vector2_merge, llmats2[1], "vec2_merge",r);
    elapsed_time(&approxChol_vector2_opt, llmats2[2], "vec2_opt", r);
    elapsed_time(&approxChol_vector2_opt2, llmats2[3], "vec2_opt2", r);
    elapsed_time(&approxChol_vector2_opt3, llmats2[4], "vec2_opt3", r);
    elapsed_time(&approxChol_vector2_opt4, llmats2[5], "vec2_opt4", r);
    elapsed_time(&approxChol_vector2_mergerand, llmats2[6], "vec2_merge_rand",r);
    elapsed_time(&approxChol_vector2_mergerand_simd, llmats2[7], "vec2_merge_rand_simd", r);

    std::cout << "\n";
    return 0;

}

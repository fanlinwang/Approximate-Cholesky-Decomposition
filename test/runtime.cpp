#include "approxChol.hpp"
#include "test.h"
#include "matrix.hpp"
#include <algorithm>
#include <iostream>
#include <vector>

int main(int argc, char **argv){

    int v = atoi(argv[1]);
    int e = atoi(argv[2]);
    int r = atoi(argv[3]);

    SparseMatrix A(v, e);
    LLMatOrd llmat = LLMatOrd(A);
    std::cout << v << " " << e << " ";

    std::vector<LLMatOrd> llmats(4, llmat);
    elapsed_time(&approxChol, llmats[0], "Baseline", r);
    elapsed_time(&approxChol_opt, llmats[1], "Inline", r);
    elapsed_time(&approxChol_opt2, llmats[2], "Inline+simd", r);

    LLMatOrd_vector2 llmat2 = LLMatOrd_vector2(A);
    std::vector<LLMatOrd_vector2> llmats2(9, llmat2);
    elapsed_time(&approxChol_vector2, llmats2[0], "Vector",r);
    elapsed_time(&approxChol_vector2_merge, llmats2[1], "Vector+Merge",r);
    elapsed_time(&approxChol_vector2_opt, llmats2[2], "Vector+Merge+simd(jkswap)", r);
    elapsed_time(&approxChol_vector2_opt2, llmats2[3], "Vector+Merge+simd(jkswap+sampling)", r);
    elapsed_time(&approxChol_vector2_opt3, llmats2[4], "Vector+Merge+simd(aligned+jkswap+sampling)", r);
    elapsed_time(&approxChol_vector2_opt4, llmats2[5], "Vector+Merge+simd(aligned+jkswap+sampling)+precompute_csum", r);
    elapsed_time(&approxChol_vector2_mergerand, llmats2[6], "Vector+Merge+pcgrand",r);
    elapsed_time(&approxChol_vector2_mergerand_simd, llmats2[7], "Vector+Merge+pcgrand+simd(jkswap)", r);

    std::cout << "\n";

    return 0;

}

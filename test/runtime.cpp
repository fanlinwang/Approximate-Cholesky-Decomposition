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

    elapsed_time(approxChol, llmat, "baseline", r);
    elapsed_time(approxChol_opt, llmat, "inline", r);
    elapsed_time(approxChol_opt, llmat, "inline SIMD", r);

    LLMatOrd_vector2 llmat2 = LLMatOrd_vector2(A);
    elapsed_time(approxChol_vector2, llmat2, "vector2",r);
    elapsed_time(approxChol_vector2_merge, llmat2, "vec2_merge",r);
    elapsed_time(approxChol_vector2_opt, llmat2, "vec2 merge simd(jkswap)", r);
    elapsed_time(approxChol_vector2_opt2, llmat2, "vect2 opt2", r);
    elapsed_time(&approxChol_vector2_opt3, llmat2, "approxChol vec2 merge simd(aligned jkswap sampling)", r);
    elapsed_time(&approxChol_vector2_opt4, llmat2, "approxChol vec2 merge simd(aligned jkswap sampling) precompute csum", r);
    elapsed_time(approxChol_vector2_mergerand, llmat2, "approxChol_vector2_mergerand",r);
    elapsed_time(&approxChol_vector2_mergerand_simd, llmat2, "approxChol vec2 merge rand simd", r);

    return 0;

}

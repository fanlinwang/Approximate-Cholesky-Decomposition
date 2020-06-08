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

    // std::vector<LLMatOrd> llmats(4, llmat);
    // elapsed_time(&approxChol, llmats[0], "Baseline", r);
    // elapsed_time(&approxChol_opt, llmats[1], "Inline", r);
    // elapsed_time(&approxChol_opt2, llmats[2], "Inline+simd", r);

    LLMatOrd_vector2 llmat2 = LLMatOrd_vector2(A);
    std::vector<LLMatOrd_vector2> llmats2(9, llmat2);
    elapsed_time(&approxChol_vector2, llmats2[0], "Vec",r);
    elapsed_time(&approxChol_vector2_merge, llmats2[1], "VecMg",r);
    // elapsed_time(&approxChol_vector2_opt, llmats2[2], "VecMgSIMD1", r);
    // elapsed_time(&approxChol_vector2_opt2, llmats2[3], "VecMgSIMD2", r);
    elapsed_time(&approxChol_vector2_opt3, llmats2[4], "VecMgSIMD3", r);
    elapsed_time(&approxChol_vector2_opt4, llmats2[5], "VecMgSIMD3csum", r);
    elapsed_time(&approxChol_vector2_mergerand, llmats2[6], "VecMgRand",r);
    elapsed_time(&approxChol_vector2_mergerand_simd, llmats2[7], "VecMgRandSIMD", r);
    elapsed_time(&approxChol_vector2_merge_search,llmats2[7], "VecMgSearch", r);
    elapsed_time(&approxChol_vector2_merge_search_opt, llmats2[7], "VecMgSearchSIMD", r);
    elapsed_time(&approxChol_vector2_merge_search_opt2, llmats2[7], "VecMgSearchSIMD2", r);
    
    LLMatOrd_vector2_struct llmat3 = LLMatOrd_vector2_struct(A);
    std::vector<LLMatOrd_vector2_struct> llmats3(3, llmat3);
    elapsed_time(&approxChol_vector2_struct_merge, llmats3[0], "VecStructMg", r);
    elapsed_time(&approxChol_vector2_struct_merge_simd, llmats3[1], "VecStructMgSIMD", r);
    elapsed_time(&approxChol_vector2_struct_merge_search, llmats3[2], "VecStructMgSearch", r);

    std::cout << "\n";

    return 0;

}

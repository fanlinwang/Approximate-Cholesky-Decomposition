#include "approxChol.hpp"
#include "test.h"
#include "matrix.hpp"
#include <algorithm>
#include <iostream>
#include <vector>
#include "search.hpp"
#include "random"

float get_random()
{
    static std::default_random_engine e;
    static std::normal_distribution<> dis(0, 1); 
    return dis(e);
}

int main(int argc, char **argv){

    int v = atoi(argv[1]);
    int e = atoi(argv[2]);
    int r = atoi(argv[3]);

    const size_t bigger_than_cachesize = 10 * 30720 * 1024;
    long *p = new long[bigger_than_cachesize];
    for (int i = 0; i < bigger_than_cachesize; i++)
        p[i] = get_random();

    SparseMatrix A(v, e);
    LLMatOrd llmat = LLMatOrd(A);
    std::cout << v << " " << e << " ";

    std::vector<LLMatOrd> llmats(3, llmat);
    elapsed_time(&approxChol, llmats[0], "Baseline", r);
    elapsed_time(&approxChol_opt, llmats[1], "Inline", r);
    elapsed_time(&approxChol_opt2, llmats[2], "Inline+simd", r);

    LLMatOrd_vector2 llmat2 = LLMatOrd_vector2(A);
    std::vector<LLMatOrd_vector2> llmats2(11, llmat2);
    elapsed_time(&approxChol_vector2, llmats2[0], "Vec",r);
    elapsed_time(&approxChol_vector2_merge, llmats2[1], "VecMg",r);
    // elapsed_time(&approxChol_vector2_opt, llmats2[2], "VecMgSIMD1", r);
    // elapsed_time(&approxChol_vector2_opt2, llmats2[3], "VecMgSIMD2", r);
    elapsed_time(&approxChol_vector2_opt3, llmats2[2], "VecMgSIMD3", r);
    // elapsed_time(&approxChol_vector2_opt4, llmats2[3], "VecMgSIMD3csum", r);
    // elapsed_time(&approxChol_vector2_mergerand, llmats2[4], "VecMgRand",r);
    // elapsed_time(&approxChol_vector2_mergerand_simd, llmats2[5], "VecMgRandSIMD", r);
    elapsed_time(&approxChol_vector2_merge_search,llmats2[6], "VecMgSearch", r);
    elapsed_time(&approxChol_vector2_merge_search_opt2, llmats2[8], "VecMgSearchSIMD", r);

    
    
    LLMatOrd_vector2_struct llmat3 = LLMatOrd_vector2_struct(A);
    std::vector<LLMatOrd_vector2_struct> llmats3(5, llmat3);
    elapsed_time(&approxChol_vector2_struct_merge, llmats3[0], "VecStructMg", r);
    elapsed_time(&approxChol_vector2_struct_merge_simd, llmats3[1], "VecStructMgSIMD", r);
    elapsed_time(&approxChol_vector2_struct_merge_search, llmats3[2], "VecStructMgSearch", r);
    elapsed_time(&approxChol_vector2_struct_merge_search_simd, llmats3[3], "VecStructMgSearchSIMD", r);
    // elapsed_time(&approxChol_vector2_struct_merge_search_nofree, llmats3[4], "VecStructMgSearch_nofree", r);
    std::cout << "\n";

   
    return 0;

}

#include "approxChol.hpp"
#include "test.h"
#include "matrix.hpp"
#include <algorithm>
#include <iostream>
#include <vector>
#include <random>

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

    SparseMatrix A(v, e);
    // int v = 4;
    // int e = 4;
    // std::vector<std::vector<double> > matrix{
    //     {0,1,1,0},
    //     {1,0,0,1},
    //     {1,0,0,1},
    //     {0,1,1,0}};
    // SparseMatrix A(matrix);
    SparseMatrix lap_A;
    laplacian(A, lap_A);
    std::vector<std::vector<Tval>> sol(5, std::vector<Tval>(v, 0.0));

    SolverParameter para;
    para.maxits = 1000;
    para.maxtime = 10000;
    para.verbose = true;
    para.tolerance = 1e-15;

    std::cout << v << " " << e << " ";
    std::vector<Tval> b(v, 0.0);
    for (auto& el: b)
        el = get_random();
    Tval m = mean(b);
    for (auto& el : b)
        el -= m;

    LLMatOrd_vector2 llmat2 = LLMatOrd_vector2(A);
    std::vector<LLMatOrd_vector2> llmats2(5, llmat2);

    std::cout << "\nVec: ";
    elapsed_time(&approxChol_vector2, llmats2[0], "Vec", r);
    LDLinv ldli1 = approxChol_vector2(llmats2[1]);

    auto start = std::chrono::steady_clock::now();
    sol[1] = LDLsolver(ldli1, b);
    sol[1] = pcg(lap_A, b, LDLsolver, sol[1], ldli1, para);
    auto end = std::chrono::steady_clock::now();
    std::cout <<  "pcg time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()/(float)(r)/1e9 << " ";

    


    std::cout << "\nVec_nosort: ";
    elapsed_time(&approxChol_nosort, llmats2[2], "Vec", r);
    LDLinv ldli2 = approxChol_nosort(llmats2[2]);

    start = std::chrono::steady_clock::now();
    sol[2] = LDLsolver(ldli2, b);
    sol[2] = pcg(lap_A, b, LDLsolver, sol[2], ldli2, para);
    end = std::chrono::steady_clock::now();
    std::cout << "pcg time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()/(float)(r)/1e9 << " ";

    LDLinv ldli3 = approxChol_vector2_opt3(llmats2[3]);
    sol[3] = LDLsolver(ldli3, b);
    sol[3] = pcg(lap_A, b, LDLsolver, sol[3], ldli3, para);

    // TODO: not correct
    // LDLinv ldli4 = approxChol_vector2_opt4(llmats2[4]);
    // sol[4] = LDLsolver(ldli4, b);
    // sol[4] = pcg(lap_A, b, LDLsolver, sol[4], ldli4, para);


    // LLMatOrd_vector2_struct llmat3 = LLMatOrd_vector2_struct(A);
    // std::vector<LLMatOrd_vector2_struct> llmats3(2, llmat3);

    // std::cout << "\nVecStructMgSearchSIMD: ";
    // elapsed_time(&approxChol_vector2_struct_merge_search_simd, llmats3[0], "Vec", r);
    // LDLinv ldli3 = approxChol_vector2_struct_merge_search_simd(llmats3[0]);

    // std::cout << " pcg: \n";
    // start = std::chrono::steady_clock::now();
    // sol[3] = LDLsolver(ldli3, b);
    // sol[3] = pcg(lap_A, b, LDLsolver, sol[3], ldli3, para);
    // end = std::chrono::steady_clock::now();
    // std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()/(float)(r)/1e9 << " ";

    std::cout << std::endl;
    return 0;

}

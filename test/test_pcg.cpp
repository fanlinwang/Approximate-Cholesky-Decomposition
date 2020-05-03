#include <iostream>
#include "approxChol.hpp"
#include "test.h"


int main(){

    // std::vector<std::vector<double> > matrix{
    //     {0,1,1,0},
    //     {1,0,0,1},
    //     {1,0,0,1},
    //     {0,1,1,0}};

    std::cout << "construct A" << std::endl;
    SparseMatrix A(4, 4);
    std::cout << A;
    // b is a normalized vector
    std::vector<Tval> b = {-0.5270249460796002, 0.7282324157897637, 
                           0.3227082222861411, -0.5239156919963046};

    // test sparse matrix * a vector:
    // should be 1.050, -1.050, -1.050, 1.050
    std::cout << "test sparse matrix multiplication: \n";
    std::vector<Tval> mul = A*b;
    for (auto m: mul)
        std::cout << m << " ";
    std::cout << "\n";

    // laplancian of A. 
    // std::vector<std::vector<double> > la{
    //     {2,-1,-1,0},
    //     {-1,2,0,-1},
    //     {-1,0,2,-1},
    //     {0,-1,-1,2}};

    SparseMatrix lap_A;
    laplacian(A, lap_A);
    std::cout << "laplacian of A" <<std::endl;
    std::cout << lap_A;

    LLMatOrd llmat = LLMatOrd(A);
    std::cout << "LLmat of A" <<std::endl;
    std::cout << llmat << std::endl;
    LDLinv ldli(llmat);
    ldli.col = {0,1,2};
    ldli.colptr = {0,2,4,5};
    ldli.rowval = {1,2,2,3,3};
    ldli.fval = {0.5, 1.0, 0.3333333333333333, 1.0, 1.0};
    ldli.d = {0.5, 0.6666666666666667, 1.3333333333333335, 0.0};

    std::vector<Tval> sol = LDLsolver(ldli, b);
    SolverParameter para;

    std::vector<Tval> iterative_sol = pcg(lap_A, b, LDLsolver, sol, ldli, para);
    for (int i = 0; i < 4; i++)
        std::cout << sol[i] <<  "," << iterative_sol[i] << "\n"; 
    return 0;
}

#include <iostream>
#include "approxChol.hpp"
#include "test.h"

#define TOLERANCE 0.01

int main(){

    std::vector<std::vector<double> > matrix{
        {0,1,1,0},
        {1,0,0,1},
        {1,0,0,1},
        {0,1,1,0}};

    SparseMatrix A(matrix);
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
    std::vector<std::vector<double> > la{
        {2,-1,-1,0},
        {-1,2,0,-1},
        {-1,0,2,-1},
        {0,-1,-1,2}};

    SparseMatrix lap_A(la);

    LLMatOrd llmat = LLMatOrd(A);
    LDLinv ldli(llmat);
    ldli.col = {0,1,2};
    ldli.colptr = {0,2,4,5};
    ldli.rowval = {1,2,2,3,3};
    ldli.fval = {0.5, 1.0, 0.3333333333333333, 1.0, 1.0};
    ldli.d = {0.5, 0.6666666666666667, 1.3333333333333335, 0.0};

    std::vector<Tval> sol = LDLsolver(ldli, b);
    std::vector<Tind> pcgIts;
    std::vector<Tval> iterative_sol = pcg(lap_A, b, LDLsolver, sol, ldli,      
                                          1e-10, 1000, 1e5, 1, pcgIts, 0);
    for (int i = 0; i < 4; i++)
        std::cout << sol[i] <<  "," << iterative_sol[i] << "\n"; 
    return 0;
}

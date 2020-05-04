#include <iostream>
#include "approxChol.hpp"
#include "test.h"
#include "matrix.hpp"

#define TOLERANCE 0.01

int main(){
    SolverParameter para;
    para.verbose = true;
    para.tolerance = 1e-8;
    
    std::vector<std::vector<double> > matrix{
        {0,1,1,0},
        {1,0,0,1},
        {1,0,0,1},
        {0,1,1,0}};
    
    SparseMatrix A(matrix);
    SparseMatrix lap_A(matrix);
    // b is a normalized vector
    std::vector<Tval> b = {-0.5270249460796002, 0.7282324157897637, 
                           0.3227082222861411, -0.5239156919963046};

    std::vector<Tval> sol;
    approxchol_lapGiven(A, lap_A, b, sol, para);
    return 0;

}
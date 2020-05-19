#include <iostream>
#include "approxChol.hpp"
#include "test.h"
#include "matrix.hpp"

#define TOLERANCE 0.01

int main(){

    std::vector<std::vector<double> > matrix{
        {0,1,1,0},
        {1,0,0,1},
        {1,0,0,1},
        {0,1,1,0}};

    SparseMatrix A(matrix);

    // b is a normalized vector
    std::vector<Tval> b = {0.3777397416887483, -1.012775866312444, 
                           1.017599128056379, -0.38256300343268407};
    LLMatOrd llmat = LLMatOrd(A);

    LDLinv ldli(llmat);
    ldli.col = {0,1,2};
    ldli.colptr = {0,2,4,5};
    ldli.rowval = {1,2,2,3,3};
    ldli.fval = {0.5, 1.0, 0.3333333333333333, 1.0, 1.0};
    ldli.d = {0.5, 0.6666666666666667, 1.3333333333333335, 0.0};
 
    std::vector<Tval> julia_sol = {0.1894727785623661, -0.5069908408742141, 
                                   0.508196656310198, -0.19067859399834997};

    std::vector<Tval> sol = LDLsolver(ldli, b);
    std::vector<std::vector<double> > la{
        {2,-1,-1,0},
        {-1,2,0,-1},
        {-1,0,2,-1},
        {0,-1,-1,2}};
    SparseMatrix lap_A(la);
    SolverParameter para;
    std::vector<Tval> iterative_sol = pcg(lap_A, b, LDLsolver, sol, ldli,      
                                          para);
    for (int i = 0; i < 4; i++)
        std::cout << sol[i] <<  "," << iterative_sol[i] << "\n"; 
    return 0;
}

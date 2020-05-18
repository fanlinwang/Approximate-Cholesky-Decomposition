#include <iostream>
#include "approxChol.hpp"
#include "test.h"
#include "conditionNumber.hpp"

#define TOLERANCE 0.01

int main(){

    std::vector<std::vector<double> > matrix{
        {0,1,1,0},
        {1,0,0,1},
        {1,0,0,1},
        {0,1,1,0}};
    std::vector<std::vector<double> > matrix2{
        {0,1,1,1,0},
        {1,0,1,1,1},
        {1,1,0,0,0},
        {1,1,0,0,0},
        {0,1,0,0,0}};
    // SparseMatrix A(matrix);
    SparseMatrix A(1000, 5000, 1);
    // std::cout << A;
    SparseMatrix la = adj2lap(A);
    std::cout << la;
    LLMatOrd llmat = LLMatOrd(A);

    LDLinv ldli = approxChol(llmat);
    std::cout << ldli << std::endl;

    // std::vector<Tval> b = {0.3777397416887483, -1.012775866312444, 
    //                        1.017599128056379, -0.38256300343268407};
    // std::vector<Tval> b = {0.3, 0.2, -0.1, -0.4};
    // std::vector<Tval> b = {0.3, 0.2, 0.2, -0.1, -0.6};
    // std::vector<Tval> b = {1.0, 1.0, 1.0, 1.0, 1.0};
    // std::vector<Tval> sol = LDLsolver(ldli, b);
    std::vector<Tval> b = {-0.46397, 0.0910944, -0.572403, 0.504353, 0.440926};


    Tval ratio = ApproxCholValidation(A, ldli, 1e-4);

    return 0;
}
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
    SparseMatrix A(100, 500, 1);
    // LLMatOrd llmat = LLMatOrd(A);

    // LDLinv ldli = approxChol(llmat);
    LLMatOrd_vector2_struct llmat(A);

    LDLinv ldli = approxChol_vector2_struct_merge(llmat);
    // std::cout << ldli << std::endl;

    Tval ratio = ApproxCholValidation(A, ldli, 1e-4);

    return 0;
}
#include <iostream>
#include "approxChol.hpp"
#include "test.h"
#include "matrix.cpp"

#define TOLERANCE 0.01

int main(){

    std::vector<std::vector<double> > matrix{
        {0,1,1,0},
        {1,0,0,1},
        {1,0,0,1},
        {0,1,1,0}};

    SparseMatrix A(matrix);

    // b is a normalized vector
    std::vector<Tval> b = {0.35177873890348577, 0.873679423502436, 
                           -1.3446358004528314, 0.11917763804690962};
    LLMatOrd llmat = LLMatOrd(A);
    LDLinv ldli = approxChol(llmat);

    std::vector<Tval> julia_sol = {0.11701982233294347, 0.4957092588700174, 
                                    -0.6134483531076161, 0.0007192719046552608};
    std::vector<Tval> sol = LDLsolver(ldli, b);

    Tval error = 0;
    for (int i = 0; i < 4; i++){
        error += julia_sol[i]-sol[i];
    }
    test(error < TOLERANCE);

    return 0;
}

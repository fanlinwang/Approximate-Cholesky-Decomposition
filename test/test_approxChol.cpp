#include <iostream>
#include "approxChol.hpp"
// #include "common.h"
// #include "approxChol.cpp"
#include "test.h"
#include "matrix.hpp" 

int main(){
    // test approxChol
    LLMatOrd llmat = simple5x5();
    llmat.n = 5;
    Tind m = 6; // m non-zero values in total
    llmat.cols = std::vector<Tind>(llmat.n);
    llmat.cols[0] = 2;
    llmat.cols[1] = 5;
    llmat.cols[2] = -1;
    llmat.cols[3] = -1;
    llmat.cols[4] = -1;
    std::vector<LLord> lles = std::vector<LLord>(m);
    lles[0].row = 1; lles[0].next = -1; lles[0].val = 1.0;
    lles[1].row = 2; lles[1].next = 0; lles[1].val = 1.0;
    lles[2].row = 3; lles[2].next = 1; lles[2].val = 1.0;
    lles[3].row = 2; lles[3].next = -1; lles[3].val = 1.0;
    lles[4].row = 3; lles[4].next = 3; lles[4].val = 1.0;
    lles[5].row = 4; lles[5].next = 4; lles[5].val = 1.0;
    llmat.lles = lles;

    LDLinv ldli = approxChol(llmat);
    std::cout << ldli << std::endl;
    std::cout << "ApproxChol tested!" << std::endl;
    return 0;
}
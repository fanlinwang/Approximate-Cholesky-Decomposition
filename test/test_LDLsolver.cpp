#include <iostream>
#include "approxChol.hpp"
#include "test.h"
#include "matrix.cpp"


int main(){
    LLMatOrd llmat = simple5x5();
    
    // TODO:
    // std::vector<Tval> b = {};
    // std::vector<Tval> solution = {};

    // const LDLinv ldli = LDLinv(llmat); 
    // std::vector<Tval> s = LDLsolver(ldli, b);
    // test(solution == s);

    return 0;
}

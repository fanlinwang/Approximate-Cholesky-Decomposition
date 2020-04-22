#include <iostream>
#include "approxChol.hpp"
// #include "common.h"
// #include "approxChol.cpp"
#include "test.h"
#include "matrix.cpp"

int main(){
    LLMatOrd llmat = simple5x5();
    print_ll_col(llmat, 1);
    return 0;
}
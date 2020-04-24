#include <iostream>
#include "approxChol.hpp"
// #include "common.h"
// #include "approxChol.cpp"
#include "test.h"
#include "matrix.hpp" 

int main(){
    // test print_ll_col, get_ll_col, compresCol
    LLMatOrd llmat = simple5x5();
    print_ll_col(llmat, 1);
    
    std::vector<LLcol> colspace;
    int len = get_ll_col(llmat, 1, colspace);
    LLcol e1; e1.row = 3; e1.ptr = 8; e1.cval = -0.5; // add multi-edges
    LLcol e2; e2.row = 2; e2.ptr = 9; e2.cval = -0.3;
    colspace.push_back(e1); colspace.push_back(e2); len += 2;

    std::cout << "Before being compressed: \n";
    for (int i = 0; i < len; i++){
        std::cout << "row = " << colspace[i].row 
            << "\t ptr = " << colspace[i].ptr
            << "\t val = " << colspace[i].cval << "\n";
    }

    std::cout << "After being compressed: \n";
    len = compressCol(colspace, len);
    for (int i = 0; i < len; i++){
        std::cout << "row = " << colspace[i].row 
            << "\t ptr =" << colspace[i].ptr
            << "\t val = " << colspace[i].cval << "\n";
    }

    return 0;
}
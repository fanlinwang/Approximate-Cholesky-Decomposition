#include "approxCholTypes.h"
#include <iostream>

/*
  Print a column in an LLMatOrd matrix.
  This is here for diagnostics.
*/
void print_ll_col(LLMatOrd llmat, int i) {
    Tind ptr = llmat.cols[i];
    while (ptr != 0)
    {
        LLord ll = llmat.lles[ptr];
        std::cout << "col " << i << ", row " << ll.row << " : " << ll.val << std::endl;

        ptr = ll.next;
    }
}


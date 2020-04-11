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


int get_ll_col(LLMatOrd llmat, int i, std::vector<LLcol> &colspace) {
    
    Tind ptr = llmat.cols[i];
    int len = 0;

    while (ptr != 0) {
        LLcol item = {llmat.lles[ptr].row, ptr, llmat.lles[ptr].val};

        if (len >= colspace.size()) {
            colspace.push_back(item);
        } else {
            colspace[len] = item;
        }

        len += 1;
        ptr = llmat.lles[ptr].next;
    }

    return len;
}


bool cmp_row(const LLcol &a, const LLcol &b) {
    return a.row < b.row;
}


bool cmp_val(const LLcol &a, const LLcol &b) {
    return a.cval < b.cval;
}


Tind compressCol(std::vector<LLcol> &colspace, int len) {

    // sort colspace ?
    std::sort(colspace.begin(), colspace.end(), cmp_row);

    std::vector<LLcol> c = colspace;

    Tind ptr = -1;
    Tind currow = c[0].row;     // julia index start from 1?
    Tval curval = c[0].cval;
    Tind curptr = c[0].ptr;

    // not consider inbounds
    for (int i = 1; i < len; i++) {

        if (c[i].row != currow) {

            ptr += 1;
            c[ptr].row = currow;
            c[ptr].cval = curval;
            c[ptr].ptr = curptr;

            currow = c[0].row;
            curval = c[0].cval;
            curptr = c[0].ptr;

        } else {
            
            curval += c[i].val;

        }
    }

    // emit the last row
    ptr += 1;
    c[ptr].row = currow;
    c[ptr].cval = curval;
    c[ptr].ptr = curptr;

    std::sort(colspace.begin(), colspace.end(), cmp_val);

    return ptr+1; // if the returned value is used for vector length, then it should return ptr+1
}

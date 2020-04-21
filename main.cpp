#include "common.h"
// #include "approxCholTypes.h"
#include "approxChol.cpp"

LLMatOrd simple5x5(){
    /*
    0 - 1 - 4
    | X |
    2   3
    Laplacian:
    3   -1  -1  -1  0
    -1  4   -1  -1  -1
    -1  -1  2   0   0
    -1  -1  0   2   0
    0   -1  0   0   1
    Notice that only the lower half triangular of the above laplacian will be stored
    */
   LLMatOrd llmat;
   llmat.n = 5;
   Tind m = 11; // m non-zero values in total
   llmat.cols = std::vector<Tind>(llmat.n);
   llmat.cols[0] = 0;
   llmat.cols[1] = 4;
   llmat.cols[2] = 8;
   llmat.cols[3] = 9;
   llmat.cols[4] = 10;
   std::vector<LLord> lles = std::vector<LLord>(m);
   lles[0].row = 0; lles[0].next = 1; lles[0].val = 1;
   lles[1].row = 1; lles[1].next = 2; lles[1].val = 1;
   lles[2].row = 2; lles[2].next = 3; lles[2].val = 1;
   lles[3].row = 3; lles[3].next = 0; lles[3].val = 1;
   lles[4].row = 1; lles[4].next = 5; lles[4].val = 1;
   lles[5].row = 2; lles[5].next = 6; lles[5].val = 1;
   lles[6].row = 3; lles[6].next = 7; lles[6].val = 1;
   lles[7].row = 4; lles[7].next = 0; lles[7].val = 1;
   lles[8].row = 2; lles[8].next = 0; lles[8].val = 1;
   lles[9].row = 3; lles[9].next = 0; lles[9].val = 1;
   lles[10].row = 4; lles[10].next = 0; lles[10].val = 1;
   llmat.lles = lles;
   return llmat;
}

int main(){
    LLMatOrd llmat = simple5x5();
    print_ll_col(llmat, 1);
}
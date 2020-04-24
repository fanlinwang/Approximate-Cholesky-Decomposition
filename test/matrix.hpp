#include "approxCholTypes.h"

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
   Tind m = 6; // m non-zero values in total
   llmat.cols = std::vector<Tind>(llmat.n);
   llmat.cols[0] = 0;
   llmat.cols[1] = 3;
   llmat.cols[2] = -1;
   llmat.cols[3] = -1;
   llmat.cols[4] = -1;
   std::vector<LLord> lles = std::vector<LLord>(m);
   lles[0].row = 1; lles[0].next = 1; lles[0].val = -1.0;
   lles[1].row = 2; lles[1].next = 2; lles[1].val = -1.0;
   lles[2].row = 3; lles[2].next = -1; lles[2].val = -1.0;
   lles[3].row = 2; lles[3].next = 4; lles[3].val = -1.0;
   lles[4].row = 3; lles[4].next = 5; lles[4].val = -1.0;
   lles[5].row = 4; lles[5].next = -1; lles[5].val = -1.0;
   llmat.lles = lles;
   return llmat;
}
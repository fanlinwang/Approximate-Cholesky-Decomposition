#include <iostream>
#include "approxChol.hpp"
#include "test.h"
#include "matrix.hpp"

double a[10000];

int main(){

    int rep = 100;
    for (int j = 0; j < rep; ++j) {
        for (int i = 0; i < 10000; ++i){
            a[i] = i*0.001;
            a[i] *= a[i];
        }
    }
/*
    int v = atoi(argv[1]);
    int e= atoi(argv[2]);

    SparseMatrix A(v, e);
    LLMatOrd llmat = LLMatOrd(A);
    std::cout << "Created random sparse matrix A.\n";

    flops = approxChol_count(llmat);
*/
    return 0;

}
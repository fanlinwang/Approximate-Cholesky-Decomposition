#include "approxChol.hpp"
#include "common.h"
#include <immintrin.h>
#include <iostream>


void forward(const LDLinv& ldli, std::vector<Tval>& y) {
    int len = ldli.col.size(); 
    for (int ii = 0; ii < len; ii++){ 
        Tind i = ldli.col[ii]; 

        Tind j0 = ldli.colptr[ii]; 
        Tind j1 = ldli.colptr[ii+1]-Tind(1);

        Tval yi = y[i];

        for (Tind jj = j0; jj <= j1-(Tind)1; jj++){
            Tval j = ldli.rowval[jj];
            y[j] += ldli.fval[jj] * yi;
            yi *= ((Tval)(1)-ldli.fval[jj]);
        }
        Tval j = ldli.rowval[j1];
        y[j] += yi;
        y[i] = yi;
    }
} 

void backward(const LDLinv& ldli, std::vector<Tval>& y) {
    int len = ldli.col.size(); 
    for (int ii = len-1; ii >= 0; ii--) {
        Tind i = ldli.col[ii];

        Tind j0 = ldli.colptr[ii];
        Tind j1 = ldli.colptr[ii+1]-Tind(1);

        Tval j = ldli.rowval[j1];
        Tval yi = y[i];
        yi = yi + y[j];

        for (Tind jj = j1-(Tind)(1); jj >= j0; jj--) {
            j = ldli.rowval[jj];
            yi = ((Tind)1-ldli.fval[jj])*yi + ldli.fval[jj]*y[j];
        }
        y[i] = yi;
    }
} 

Tval mean(const std::vector<Tval>& y){
    Tval sum = 0.0;
    for (auto elem: y){
        sum += elem;
    }
    return sum/y.size();
}

std::vector<Tval> LDLsolver(const LDLinv& ldli, const std::vector<Tval>& b){
    std::vector<Tval> y(b.begin(), b.end());

    forward(ldli, y);

    for (int i = 0; i < ldli.d.size(); i++){
        if (ldli.d[i] != 0)
            y[i] /= ldli.d[i];
    }

    backward(ldli, y);

    Tval mu = mean(y);
    for (auto elem: y){
        elem = elem - mu;
    }

    return y;
}

void approxchol_lapGiven(const SparseMatrix& A, const
                                      std::vector<Tval>& b,  std::vector<Tval>& sol, bool verbose) 
{
    // la = lap(a); ??

    // LLMatOrd llmat = LLMatOrd(a);
    // LDLinv ldli = approxChol(llmat);
    // sol = LDLsolver(ldli, b);

    if (verbose){
        // std::cout << "Ratio of operator edges to original edges: " << 2 * ldli.fval.size() / nnz(a) << "\n";
        // std::cout << "ratio of max to min diagonal of laplacian : " << maximum(diag(la))/minimum(diag(la))) << "\n";
    }
}

/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions()
{
    add_function(&approxchol_lapGiven, "Base Solver", 1);
    // add_function(&approxchol_lapGiven1, "Fast WHT", 1);
}

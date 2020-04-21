#include "approxChol.hpp"
#include "common.h"
#include <immintrin.h>
#include <iostream>
#include <random>

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
            
            curval += c[i].cval;

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

LDLinv approxChol(LLMatOrd a) {
    auto n = a.n;

    // need to make custom one without col info later
    LDLinv ldli(a);
    Tind ldli_row_ptr = 0;

    std::vector<Tval> d(n, 0);

    std::vector<LLcol> colspace(n);
    std::vector<Tval> cumspace(n);

    // random engine and distribution
    std::default_random_engine engine;
    std::uniform_real_distribution<Tval> u(0.0, 1.0);

    for (long i = 0; i <= n-2; i++) {

        ldli.col[i] = i;
        ldli.colptr[i] = ldli_row_ptr;

        int len = get_ll_col(a, i, colspace);
        len = compressCol(colspace, len);

        Tval csum = 0;
        for (int ii = 0; ii < len; ii++) {
            csum += colspace[ii].cval;
            cumspace[ii] = csum;
        }
        Tval wdeg = csum;

        Tval colScale = 1;

        for (int joffset = 0; joffset <= len-2; joffset++) {

            LLcol llcol = colspace[joffset];
            Tval w = llcol.cval * colScale;
            Tind j = llcol.row;

            Tval f = w/wdeg;

            Tval r = u(engine);
            r = r * (csum - cumspace[joffset]) + cumspace[joffset];

            auto cumspace_last = cumspace.begin();
            std::advance(cumspace_last, len);
            int koff = std::distance(cumspace.begin(), std::lower_bound(cumspace.begin(), cumspace_last, r));

            Tind k = colspace[koff].row;

            Tval newEdgeVal = w*(1-f);

            // create edge (j,k) with newEdgeVal
            // do it by reassigning ll
            if (j < k) {    // put it in col j
                Tind jhead = a.cols[j];
                a.lles[llcol.ptr].row = k;
                a.lles[llcol.ptr].next = jhead;
                a.lles[llcol.ptr].val = newEdgeVal;
                a.cols[j] = llcol.ptr;
            } else {        // put it in col k
                Tind khead = a.cols[k];
                a.lles[llcol.ptr].row = j;
                a.lles[llcol.ptr].next = khead;
                a.lles[llcol.ptr].val = newEdgeVal;
                a.cols[k] = llcol.ptr;
            }

            colScale *= 1 - f;
            wdeg = wdeg - 2*w + w*w/wdeg;

            ldli.rowval.push_back(j);
            ldli.fval.push_back(f);
            ldli_row_ptr += 1;
        }

        LLcol llcol = colspace[len];
        Tval w = llcol.cval * colScale;
        Tind j = llcol.row;

        ldli.rowval.push_back(j);
        ldli.fval.push_back(1);
        ldli_row_ptr += 1;

        d[i] = w;
    }

    ldli.colptr[n-1] = ldli_row_ptr;
    ldli.d = d;

    return ldli;
}

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

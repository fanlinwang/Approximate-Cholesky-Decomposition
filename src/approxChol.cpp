#include "approxChol.hpp"
#include "common.h"
#include <immintrin.h>
#include <iostream>
#include <random>
#include <cmath>
#include <ctime>
#include <limits>


void print_ll_col(LLMatOrd llmat, int i) {
    Tind ptr = llmat.cols[i];
    while (ptr != -1)
    {
        LLord ll = llmat.lles[ptr];
        std::cout << "col " << i << ", row " << ll.row << " : " << ll.val << std::endl;

        ptr = ll.next;
    }
}

int get_ll_col(LLMatOrd llmat, int i, std::vector<LLcol> &colspace) {
    
    Tind ptr = llmat.cols[i];
    int len = 0;

    while (ptr != -1) {
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

    // sort colspace ? DY: just checked, yes! 
    std::sort(colspace.begin(), colspace.begin() + len, cmp_row);

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

            currow = c[i].row;
            curval = c[i].cval;
            curptr = c[i].ptr;

        } else {
            
            curval += c[i].cval;

        }
    }

    // emit the last row
    ptr += 1;
    c[ptr].row = currow;
    c[ptr].cval = curval;
    c[ptr].ptr = curptr;

    std::sort(c.begin(), c.begin() + ptr + 1, cmp_val);
    
    colspace = c;
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
        // flop count: sort LLcol len*log(len)?

        Tval csum = 0;
        for (int ii = 0; ii < len; ii++) {
            csum += colspace[ii].cval;
            // flop count: 1 add
            cumspace[ii] = csum;
        }
        Tval wdeg = csum;

        Tval colScale = 1;

        for (int joffset = 0; joffset <= len-2; joffset++) {

            LLcol llcol = colspace[joffset];
            Tval w = llcol.cval * colScale;
            Tind j = llcol.row;

            Tval f = w/wdeg;
            // flop count: 1 mul

            Tval r = u(engine);
            r = r * (csum - cumspace[joffset]) + cumspace[joffset];
            // flop count: 1 mul 2 add

            auto cumspace_last = cumspace.begin();
            std::advance(cumspace_last, len);
            int koff = std::distance(cumspace.begin(), std::lower_bound(cumspace.begin(), cumspace_last, r));
            // flop count: len?

            Tind k = colspace[koff].row;

            Tval newEdgeVal = w*(1-f);
            // flop count: 1 mul 1 add


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
            wdeg = wdeg - 2*w + w*f;
            // flop count: 3 mul 3 add

            ldli.rowval.push_back(j);
            ldli.fval.push_back(f);
            ldli_row_ptr += 1;
        }

        LLcol llcol = colspace[len-1];
        Tval w = llcol.cval * colScale;
        // flop count: 1 mul
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

int approxChol_count(LLMatOrd a) {
    auto n = a.n;
    int flops_count = 0;

    // need to make custom one without col info later
    LDLinv ldli(a);
    Tind ldli_row_ptr = 0;

    std::vector<Tval> d(n, 0);

    std::vector<LLcol> colspace(n);
    std::vector<Tval> cumspace(n);

    // random engine and distribution
    std::default_random_engine engine;
    std::uniform_real_distribution<Tval> u(0.0, 1.0);

    //Tval randnums[5] = {0.7527533347596496, 0.5805790891279785, 0.14277294752566538, 0.41620584940537597, 0.4559296013211689};
    //int randptr = 0;

    for (long i = 0; i <= n-2; i++) {

        ldli.col[i] = i;
        ldli.colptr[i] = ldli_row_ptr;

        int len = get_ll_col(a, i, colspace);
        len = compressCol(colspace, len);
        // flop count: sort LLcol len*log(len)?

        Tval csum = 0;
        for (int ii = 0; ii < len; ii++) {
            csum += colspace[ii].cval;
            flops_count++;
            cumspace[ii] = csum;
        }
        Tval wdeg = csum;

        Tval colScale = 1;

        for (int joffset = 0; joffset <= len-2; joffset++) {

            LLcol llcol = colspace[joffset];
            Tval w = llcol.cval * colScale;
            Tind j = llcol.row;

            Tval f = w/wdeg;
            flops_count++;
            // flop count: 1 mul

            Tval r = u(engine);
            //Tval r = randnums[randptr++];
            r = r * (csum - cumspace[joffset]) + cumspace[joffset];
            flops_count += 3;
            // flop count: 1 mul 2 add

            auto cumspace_last = cumspace.begin();
            std::advance(cumspace_last, len);
            int koff = std::distance(cumspace.begin(), std::lower_bound(cumspace.begin(), cumspace_last, r));
            // flop count: len?

            Tind k = colspace[koff].row;

            Tval newEdgeVal = w*(1-f);
            flops_count += 2;
            // flop count: 1 mul 1 add

            /*{
                std::cout << r << "\t" << cumspace[0] << " " << cumspace[1] << " " << cumspace[2] << " " << cumspace[3] << " " << cumspace[4] << std::endl;
                std::cout << "j: " << j << "\tk: " << k << "\tkoff: " << koff << std::endl;
            }*/

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
            wdeg = wdeg - 2*w + w*f;
            flops_count += 6;
            // flop count: 3 mul 3 add

            ldli.rowval.push_back(j);
            ldli.fval.push_back(f);
            ldli_row_ptr += 1;
        }

        LLcol llcol = colspace[len-1];
        Tval w = llcol.cval * colScale;
        flops_count += 1;
        // flop count: 1 mul
        Tind j = llcol.row;

        ldli.rowval.push_back(j);
        ldli.fval.push_back(1);
        ldli_row_ptr += 1;

        d[i] = w;
    }

    ldli.colptr[n-1] = ldli_row_ptr;
    ldli.d = d;

    std::cout << "#flops: " << flops_count << std::endl;
    return flops_count;
}

void forward(const LDLinv& ldli, std::vector<Tval>& y) {
    int len = ldli.col.size(); 
    for (int ii = 0; ii < len; ii++){ 
        Tind i = ldli.col[ii]; 

        Tind j0 = ldli.colptr[ii]; 
        Tind j1 = ldli.colptr[ii+1]-1;

        Tval yi = y[i];

        for (Tind jj = j0; jj <= j1-1; jj++){
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
        Tind j1 = ldli.colptr[ii+1]-1;

        Tval j = ldli.rowval[j1];
        Tval yi = y[i];
        yi = yi + y[j];

        for (Tind jj = j1-1; jj >= j0; jj--) {
            j = ldli.rowval[jj];
            yi = ((Tind)1-ldli.fval[jj])*yi + ldli.fval[jj]*y[j];
        }
        y[i] = yi;
    }
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
    for (auto& elem: y){
        elem = elem - mu;
    }

    return y;
}

void approxchol_lapGiven(SparseMatrix& A, SparseMatrix& lap_A,
                         const std::vector<Tval>& b,          
                         std::vector<Tval>& sol, SolverParameter paras) 
{
    LLMatOrd llmat = LLMatOrd(A);
    LDLinv ldli = approxChol(llmat);
    sol = LDLsolver(ldli, b);

    Tval error = norm(lap_A*sol-b)/norm(b); 

    if (paras.verbose){
        // std::cout << "Ratio of operator edges to original edges: " << 2 * ldli.fval.size() / nnz(a) << "\n";
        // std::cout << "ratio of max to min diagonal of laplacian : " << maximum(diag(la))/minimum(diag(la))) << "\n";
        std::cout << "Error after first iteration: " << error << "\n";
    }

    // iteratively solve the equation
    sol = pcg(lap_A, b, LDLsolver, sol, ldli, paras);
}


void axpy2(const Tval al, const std::vector<Tval>& p, std::vector<Tval>& x){
    int n = x.size();
    for (int i = 0; i < n; ++i)
        x[i] = x[i] + al*p[i];
}

void bzbeta(const Tval beta, const std::vector<Tval>& z, std::vector<Tval>& p){
    int n = p.size();
    for (int i = 0; i < n; ++i)
        p[i] = z[i] + beta*p[i];
}

Tval mean(const std::vector<Tval>& y){
    Tval sum = 0.0;
    for (auto& elem: y){
        sum += elem;
    }
    return sum/y.size();
}


Tval norm(const std::vector<Tval>& y){
    Tval sum = 0.0;
    for (auto& elem: y){
        sum += elem*elem;
    }
    return std::sqrt(sum);
}

Tval dot(const std::vector<Tval>& x, const std::vector<Tval>& y){
    Tval sum = 0.0;
    for (int i = 0; i < x.size(); ++i){
        sum += x[i]*y[i];
    }
    return sum;
}


std::vector<Tval> pcg(const SparseMatrix& la, const std::vector<Tval>& b,       
                      solver ldlsolver, const std::vector<Tval>& sol, 
                      const LDLinv& ldli, SolverParameter& paras)
{ 
    Tind n = b.size();
    Tval nb = norm(b);

    // If input vector is zero, quit
    if (nb == 0)
      return sol;

    std::vector<Tval> x(n, 0.0);
    std::vector<Tval> bestx(n, 0.0);
    Tval bestnr = 1.0;

    std::vector<Tval> r(b.begin(), b.end());
    std::vector<Tval> z(sol.begin(), sol.end());
    std::vector<Tval> p(sol.begin(), sol.end());
    std::vector<Tval> q(n, 0.0);

    Tval rho = dot(r, sol);
    Tval best_rho = rho;
    Tind stag_count = 0;
    time_t t1 = time(NULL);

    Tind itcnt = 0;
    Tval oldrho = rho; 
    while (itcnt < paras.maxits){
        itcnt = itcnt+1;

        q = la*p; 
        Tval pq = dot(p,q);
        // if (paras.verbose) 
        //     std::cout << "Current iteration: " << itcnt << ", pq = " 
        //     << pq << ", rho = " << rho << "\n";
        
        if (pq < EPS || pq > INF){ 
          if (paras.verbose)
            std::cout << "PCG Stopped due to small or large pq.\n";
          break;
        }

        Tval al = rho/pq;
        // if (paras.verbose) 
        //     std::cout << "al: " << al << ", nb = " 
        //     << nb << "\n";

        // the following line could cause slowdown
        if (al*norm(p) < EPS*norm(x)){
            if (paras.verbose) 
                std::cout << "PCG: Stopped due to stagnation.\n";
            break;
        }
        
        axpy2(al,p,x);
        axpy2(-al,q,r);

        Tval nr = norm(r)/nb;
        if (nr < bestnr){
            bestnr = nr;
            for (int i = 0; i < n; ++i)
                bestx[i] = x[i];
        }
                    
        if (nr < paras.tolerance) 
            break;
        z = ldlsolver(ldli, r);

        oldrho = rho;
        rho = dot(z, r); 

        if (paras.stag_test != 0) {
            if (rho < best_rho*(1-1/paras.stag_test)){
                best_rho = rho;
                stag_count = 0;
            } else {
                if (paras.stag_test > 0) {
                    if (best_rho > (1-1/paras.stag_test)*rho) {
                        stag_count += 1;
                        if (stag_count > paras.stag_test) {
                            std::cout << "PCG Stopped by stagnation test " << paras.stag_test << "\n";
                            break;
                        }
                    }
                }
            }
        }

        if (rho < EPS || rho > INF){
          if (paras.verbose)
            std::cout << "PCG Stopped due to small or large rho.\n"; 
          break;
        }


        Tval beta = rho/oldrho; 
        if (beta < EPS || beta > INF) {
          if (paras.verbose)
            std::cout << "PCG Stopped due to small or large beta.\n"; 
          break;
        }

        bzbeta(beta,p,z);

        if ((time(NULL) - t1) > paras.maxtime){
            if (paras.verbose) 
                std::cout << "PCG New stopped at maxtime.\n"; 
            break;
        }

    }

    if (paras.verbose)
        std::cout << "PCG stopped after: " << std::round((time(NULL) - t1)) << " seconds and " << itcnt << " iterations with relative error " << (norm(la*bestx-b)/norm(b)) << ".\n";

    if (paras.pcgIts.size() > 0)
        paras.pcgIts[0] = itcnt; 

    return bestx;
}

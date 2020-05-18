#include "conditionNumber.hpp"
#include "approxChol.hpp"
#include "approxCholTypes.h"
#include <random>
#include <cmath>
#include <iostream>

SparseMatrix adj2lap(SparseMatrix& a){
    SparseMatrix la;
    la.colnum = a.colnum;
    la.elems = a.elems + a.colnum;
    la.colptrs = std::vector<Tind>(la.colnum + 1);
    for (int colidx = 0; colidx < la.colnum; colidx++){
        la.colptrs.at(colidx) = a.colptrs.at(colidx) + colidx;
    }
    la.colptrs.at(la.colnum) = la.elems;
    la.rows = std::vector<Tind>(la.elems);
    la.vals = std::vector<Tval>(la.elems);
    for (int colidx = 0; colidx < la.colnum; colidx++){
        int i = a.colptrs.at(colidx);
        int j = la.colptrs.at(colidx);
        Tval deg = 0.0;
        bool above_diag = true;
        Tind diag_place = la.elems-1;
        while (i < a.elems && i < a.colptrs.at(colidx+1)){
            if (above_diag){
                if (a.rows.at(i) > colidx){
                    above_diag = false;
                    diag_place = j;
                    j++; // leave the place for degree, which should be on the diagonal
                }
            }
            deg += a.vals.at(i);
            la.vals.at(j) = -a.vals.at(i);
            la.rows.at(j) = a.rows.at(i);
            i++; j++;
        }
        la.rows.at(diag_place) = colidx;
        la.vals.at(diag_place) = deg;
    }
    return la;
}

void op1(SparseMatrix& a, const LDLinv& ldli, Tval u, std::vector<Tval>& y){
    SparseMatrix la; 
    laplacian(a, la);
    Tval mu = mean(y);
    for (auto& elem: y){
        elem = elem - mu;
    }
    std::cout << "after projection: \t";
    for (auto elem : y ){std::cout << elem << " ";} 
    std::cout << "\n";

    // std::vector<Tval> sol = LDLsolver(ldli, y);
    std::vector<Tval> ldlsol = LDLsolver(ldli, y);
    SolverParameter para;
    para.verbose = false;
    para.tolerance = 1e-6;
    para.maxits = 1000;
    std::vector<Tval> sol = pcg(la, y, LDLsolver, ldlsol, ldli, para);
    // SolverParameter para;
    // para.verbose = true;
    // para.tolerance = 1e-8;
    // std::vector<Tval> sol(y.size());
    // approxchol_lapGiven(a_copy, la, y, sol, para);

    std::cout << "after M^{-1}: \t \t";
    for (auto elem : sol ){std::cout << elem << " ";} 
    std::cout << "\n";

    y = la * sol;

    std::cout << "L M^{-1}: \t \t";
    for (auto elem : y ){std::cout << elem << " ";} 
    std::cout << "\n";
}

void op2(SparseMatrix& a, const LDLinv& ldli, Tval u, std::vector<Tval>& y){
    // SparseMatrix la = adj2lap(a);
    SparseMatrix la; 
    laplacian(a, la);
    Tval mu = mean(y);
    for (auto& elem: y){
        elem = elem - mu;
    }
    std::cout << "after projection: \t";
    for (auto elem : y ){std::cout << elem << " ";} 
    std::cout << "\n";

    // std::vector<Tval> sol = LDLsolver(ldli, y);
    std::vector<Tval> ldlsol = LDLsolver(ldli, y);
    SolverParameter para;
    para.verbose = true;
    para.tolerance = 1e-6;
    para.maxits = 1000;
    std::vector<Tval> sol = pcg(la, y, LDLsolver, ldlsol, ldli, para);
    // SolverParameter para;
    // para.verbose = true;
    // para.tolerance = 1e-8;
    // std::vector<Tval> sol(y.size());
    // approxchol_lapGiven(a_copy, la, y, sol, para);

    std::cout << "after M^{-1}: \t \t";
    for (auto elem : sol ){std::cout << elem << " ";} 
    std::cout << "\n";

    sol = la * sol;

    std::cout << "L M^{-1}: \t \t";
    for (auto elem : sol ){std::cout << elem << " ";} 
    std::cout << "\n";

    for (auto& elem : y){
        elem *= u;
    }
    int n = y.size();
    for (int i = 0; i < n; i++){
        y[i] -= sol[i];
    }

    std::cout << "out from the operator: \t";
    for (auto elem : y ){std::cout << elem << " ";} 
    std::cout << "\n";
}

Tval eigs(LinearOperator op, SparseMatrix& la, const LDLinv& ldli, Tval u,
    std::vector<Tval>& y, Tind nev, Tval tol){
    Tind n = y.size();
    Tval error = 1.0;
    Tval eigenvalue;
    Tval ynorm = 0;
    Tval diffnorm = 1;

    // normalize y
    Tval normy = norm(y);
    for (auto& elem : y){
        elem /= normy;
    }

    while (diffnorm > tol){
        
        std::vector<Tval> x = y;

        // y <- op(y)
        op(la, ldli, u, y);
        std::vector<Tval> y_unnormalized = y;
        ynorm = norm(y);
        for (auto& elem : y){
            elem /= ynorm;
        }

        // project?
        if (nev == 2){
            Tval mu = mean(y);
            for (auto& elem: y){
                elem = elem - mu;
            }
        }

        // calculate eigenvalue
        eigenvalue = 0.0;
        Tval numerator = 0.0;
        Tval denominator = 0.0;
        // std::cout << "calculating the eigenvalue:\n";
        // for (auto elem : x ){std::cout << elem << " ";} 
        // std::cout << "\n";
        // for (auto elem : y_unnormalized ){std::cout << elem << " ";} 
        // std::cout << "\n";
        for (int i = 0; i < n; i++){
            numerator += x.at(i) * y_unnormalized.at(i);
            denominator += x.at(i) * x.at(i);
        }
        eigenvalue = numerator / denominator;

        std::vector<Tval> diff(n);
        for (int i = 0; i < n; i++){
            diff[i] = eigenvalue * x[i] - y_unnormalized[i];
        }
        diffnorm = norm(diff);

        //debug
        for (int i = 0; i < n; i++){std::cout << x[i] << " ";} 
        std::cout << "\n";
        for (int i = 0; i < n; i++){std::cout << y_unnormalized[i] << " ";}
        std::cout << "\n";
        for (int i = 0; i < n; i++){std::cout << y[i] << " ";}
        std::cout << "\n eigenvalue: \t\t" << eigenvalue << "\n";
        std::cout << "Diffnorm: \t\t" << diffnorm << "\n";
    }

    return eigenvalue;
}

Tval ApproxCholValidation(SparseMatrix& a, const LDLinv& ldli, Tval tol){
    Tind n = a.colnum;

    // randomly initialize eigenvector (projection will be done inside the op function)
    std::default_random_engine engine;
    std::uniform_real_distribution<Tval> uniform(0.0, 1.0);
    std::vector<Tval> y(n);
    for (int i = 0; i < n; i++){
        y.at(i) = uniform(engine);
    }

    // Compute max eigenvalue of A * L^{-1}
    Tval upper = eigs(&op1, a, ldli, 0.0, y, 1, tol);
    // Tval upper = 1.0;

    // initialize eigenvector again for the second operator
    for (int i = 0; i < n; i++){
        y.at(i) = uniform(engine);
    }
    
    Tval lower = upper - eigs(&op2, a, ldli, upper, y, 2, tol);
    
    std::cout << "biggest eigenvalue: " << upper << "\n \n";
    std::cout << "second smallest eigenvalue: " << lower << "\n";
    // Tval lower = 1;
    return upper / lower;
}
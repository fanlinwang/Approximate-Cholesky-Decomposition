#include "conditionNumber.hpp"
#include "approxChol.hpp"
#include "approxCholTypes.h"
#include <random>
#include <cmath>
#include <iostream>

void op1(SparseMatrix& a, const LDLinv& ldli, Tval u, std::vector<Tval>& y){
    SparseMatrix la; 
    laplacian(a, la);
    Tval mu = mean(y);
    for (auto& elem: y){
        elem = elem - mu;
    }

    std::vector<Tval> ldlsol = LDLsolver(ldli, y);
    SolverParameter para;
    para.verbose = false;
    para.tolerance = 1e-6;
    para.maxits = 1000;
    std::vector<Tval> sol = pcg(la, y, LDLsolver, ldlsol, ldli, para);

    y = la * sol;
}

void op2(SparseMatrix& a, const LDLinv& ldli, Tval u, std::vector<Tval>& y){
    SparseMatrix la; 
    laplacian(a, la);
    Tval mu = mean(y);
    for (auto& elem: y){
        elem = elem - mu;
    }

    std::vector<Tval> ldlsol = LDLsolver(ldli, y);
    SolverParameter para;
    para.verbose = true;
    para.tolerance = 1e-6;
    para.maxits = 1000;
    std::vector<Tval> sol = pcg(la, y, LDLsolver, ldlsol, ldli, para);

    sol = la * sol;

    for (auto& elem : y){
        elem *= u;
    }
    int n = y.size();
    for (int i = 0; i < n; i++){
        y[i] -= sol[i];
    }
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

    // initialize eigenvector again for the second operator
    for (int i = 0; i < n; i++){
        y.at(i) = uniform(engine);
    }
    
    Tval lower = upper - eigs(&op2, a, ldli, upper, y, 2, tol);
    
    std::cout << "biggest eigenvalue: " << upper << "\n \n";
    std::cout << "second smallest eigenvalue: " << lower << "\n";
    
    return upper / lower;
}
#ifndef __CONDITION_NUMBER__
#define __CONDITION_NUMBER__

#include "approxChol.hpp"
#include "approxCholTypes.h"
#include "common.h"

typedef void(*LinearOperator)(SparseMatrix& a, const LDLinv& ldli, Tval u, std::vector<Tval>& y);

/*
from adjacency matrix to laplacian matrix
*/
SparseMatrix adj2lap(SparseMatrix& a);

/*
In op1: la * L^{-1} * y
in op2: (uI - la * L^{-1}) * y
Parameter u in op1 is just for the sake of a consistent api...
*/
void op1(SparseMatrix& la, const LDLinv& ldli, Tval u, std::vector<Tval>& y);
void op2(SparseMatrix& la, const LDLinv& ldli, Tval u, std::vector<Tval>& y);

/*
\param{op}: suppose I want to compute the eigenvalues of matrix A, then op is the operator form of A
\param{y}: randomly initialized vector (eigenvector)
\param{nev}: number of Eigenvalue. nev = 1 means biggest eigenvalue, nev = 2 means second biggest eigenvalue
\param{tol}: stopping criterion of power iteration
return the eigenvalue wanted
Notice that this function is hardcored for our problem...
*/
Tval eigs(LinearOperator op, SparseMatrix& la, const LDLinv& ldli, Tval u,
    std::vector<Tval>& y, Tind nev, Tval tol);

/*
*/
Tval ApproxCholValidation(SparseMatrix& a, const LDLinv& ldli, Tval tol);


#endif
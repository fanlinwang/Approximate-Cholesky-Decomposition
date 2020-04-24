#pragma once

#include "approxCholTypes.h"
#include "common.h"

/*
  Print a column in an LLMatOrd matrix.
  This is here for diagnostics.
*/
void print_ll_col(LLMatOrd llmat, int i);
int get_ll_col(LLMatOrd llmat, int i, std::vector<LLcol> &colspace);
bool cmp_row(const LLcol &a, const LLcol &b);
bool cmp_val(const LLcol &a, const LLcol &b);
Tind compressCol(std::vector<LLcol> &colspace, int len);
LDLinv approxChol(LLMatOrd a);
void forward(const LDLinv& ldli, std::vector<Tval>& y);
void backward(const LDLinv& ldli, std::vector<Tval>& y); 
Tval mean(const std::vector<Tval>& y);
std::vector<Tval> LDLsolver(const LDLinv& ldli, const std::vector<Tval>& b);
void approxchol_lapGiven(const SparseMatrix& A, const std::vector<Tval>& b,   
                         std::vector<Tval>& sol, Tval error, bool verbose);


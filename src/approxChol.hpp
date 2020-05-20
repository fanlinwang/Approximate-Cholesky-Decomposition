#ifndef __APPROX_CHOL__
#define __APPROX_CHOL__

#include "approxCholTypes.h"
#include "common.h"
#include <limits>
#define EPS std::numeric_limits<Tval>::epsilon()
#define INF std::numeric_limits<Tval>::max()

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
LDLinv approxChol_vector2(LLMatOrd_vector2 a);
//LDLinv approxChol_vector3(LLMatOrd_vector2 a);
LDLinv approxChol_vector2_opt(LLMatOrd_vector2 a);
LDLinv approxChol_opt(LLMatOrd a);
void approxChol_count(LLMatOrd a, int &flops_count, int &flcomp_count, int &intops_count, int &intcomp_count);
void forward(const LDLinv& ldli, std::vector<Tval>& y);
void backward(const LDLinv& ldli, std::vector<Tval>& y); 

// TODO: faster library to do the simple matrix operations?
Tval mean(const std::vector<Tval>& y);
Tval norm(const std::vector<Tval>& y);
Tval dot(const std::vector<Tval>& x, const std::vector<Tval>& y);
void axpy2(const Tval al, const std::vector<Tval>& p, std::vector<Tval>& x);
void bzbeta(const Tval beta, const std::vector<Tval>& z, std::vector<Tval>& p);
/*
  pcg: solves a symmetric linear system using preconditioner `pre`.
  `pre` can be a function or a matrix.  If a matrix, a function to solve it is created with cholFact. 
  `tol` is set to 1e-6 by default,
  `maxits` defaults to Inf
  `maxtime` defaults to Inf.  It measures seconds.
  `verbose` defaults to false
  `pcgIts` is an array for returning the number of pcgIterations.  Default is length 0, in which case nothing is returned.
  `stag_test=k` stops the code if rho[it] > (1-1/k) rho[it-k].  Set to 0 to deactivate.
*/
typedef std::vector<Tval>(*solver)(const LDLinv& ldli, 
                                   const std::vector<Tval>& b);
      
struct SolverParameter {
    Tval tolerance;
    Tind maxits;
    time_t maxtime;
    bool verbose;
    std::vector<Tind> pcgIts;
    Tind stag_test;
    SolverParameter(){
        tolerance = 1e-5;
        maxits = 10;
        maxtime = 10000;
        verbose = true;
        stag_test = 0;
    }
};

std::vector<Tval> pcg(const SparseMatrix& la, const std::vector<Tval>& b,       
                      solver ldlsolver, const std::vector<Tval>& init_sol, 
                      const LDLinv& ldli, SolverParameter& paras);

// solver and its wrapper
std::vector<Tval> LDLsolver(const LDLinv& ldli, const std::vector<Tval>& b);
void approxchol_lapGiven(SparseMatrix& A, SparseMatrix& lap_A,
                         const std::vector<Tval>& b,          
                         std::vector<Tval>& sol, SolverParameter paras);


#endif
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

namespace std {

namespace sort_helper {

template <typename _Data, typename _Order>
struct value_reference_t;

template <typename _Data, typename _Order>
struct value_t {
    _Data data;
    _Order val;
    inline value_t(_Data _data, _Order _val) : data(_data), val(_val) {}
    inline value_t(const value_reference_t<_Data,_Order>& rhs);
};

template <typename _Data, typename _Order>
struct value_reference_t {
    _Data* pdata;
    _Order* pval;
    value_reference_t(_Data* _itData, _Order* _itVal) : pdata(_itData), pval(_itVal) {}
    inline value_reference_t& operator = (const value_reference_t& rhs) { *pdata = *rhs.pdata; *pval = *rhs.pval; return *this; }
    inline value_reference_t& operator = (const value_t<_Data,_Order>& rhs) { *pdata = rhs.data; *pval = rhs.val; return *this; }
    inline bool operator < (const value_reference_t& rhs) { return *pval < *rhs.pval; }
};

template <typename _Data, typename _Order>
struct value_iterator_t :
    iterator< random_access_iterator_tag, value_t<_Data,_Order>, ptrdiff_t, value_t<_Data,_Order>*, value_reference_t<_Data,_Order> >
{
    _Data* itData;
    _Order* itVal;
    value_iterator_t(_Data* _itData, _Order* _itVal) : itData(_itData), itVal(_itVal) {}
    inline ptrdiff_t operator - (const value_iterator_t& rhs) const { return itVal - rhs.itVal; }
    inline value_iterator_t operator + (ptrdiff_t off) const { return value_iterator_t(itData + off, itVal + off); }
    inline value_iterator_t operator - (ptrdiff_t off) const { return value_iterator_t(itData - off, itVal - off); }
    inline value_iterator_t& operator ++ () { ++itData; ++itVal; return *this; }
    inline value_iterator_t& operator -- () { --itData; --itVal; return *this; }
    inline value_iterator_t operator ++ (int) { return value_iterator_t(itData++, itVal++); }
    inline value_iterator_t operator -- (int) { return value_iterator_t(itData--, itVal--); }
    inline value_t<_Data,_Order> operator * () const { return value_t<_Data,_Order>(*itData, *itVal); }
    inline value_reference_t<_Data,_Order> operator * () { return value_reference_t<_Data,_Order>(itData, itVal); }
    inline bool operator  < (const value_iterator_t& rhs) const { return itVal  < rhs.itVal; }
    inline bool operator == (const value_iterator_t& rhs) const { return itVal == rhs.itVal; }
    inline bool operator != (const value_iterator_t& rhs) const { return itVal != rhs.itVal; }
};

template <typename _Data, typename _Order>
inline value_t<_Data,_Order>::value_t(const value_reference_t<_Data,_Order>& rhs)
    : data(*rhs.pdata), val(*rhs.pval) {}

template <typename _Data, typename _Order>
bool operator < (const value_t<_Data,_Order>& lhs, const value_reference_t<_Data,_Order>& rhs) {
    return lhs.val < *rhs.pval; }

template <typename _Data, typename _Order>
bool operator < (const value_reference_t<_Data,_Order>& lhs, const value_t<_Data,_Order>& rhs) {
    return *lhs.pval < rhs.val; }

template <typename _Data, typename _Order>
void swap(value_reference_t<_Data,_Order> lhs, value_reference_t<_Data,_Order> rhs) {
    std::swap(*lhs.pdata, *rhs.pdata);
    std::swap(*lhs.pval, *rhs.pval); }


} // namespace sort_helper

} // namespace std

#endif
#include "common.h"
typedef int Tind;
typedef double Tval;

struct Edge{
    Tind x, y;
    Tval w;
    Edge(Tind x, Tind y):x(x), y(y){}
};

struct Elem{
    Tind x, y;
    Tval w;
    Elem(Tind x, Tind y, Tval w):x(x), y(y), w(w){}
};

//not weighted
struct SparseMatrix{
    std::vector<std::vector<Elem> > matrix;
    // long cols;
    // long rows;
    SparseMatrix(std::vector<Edge> edges, long rows, long cols)//suppose input edges to be sorted and x > y
    {
        matrix = std::vector<std::vector<Elem> >(cols);
        for (int i = 0; i < edges.size(); i++)
        {
            Tind x = edges[i].x;
            Tind y = edges[i].y;
            Tval w = edges[i].w;
            //Todo
            //construct sparse matrix from edges
        }
    }
};

/*
"""
  LLp elements are all in the same column.
  row tells us the row, and val is the entry.
  val is set to zero for some edges that we should remove.
  next gives the next in the column.  It poTinds to itself to terminate.
  reverse is the index Tindo lles of the other copy of this edge,
  since every edge is stored twice as we do not know the order of elimination in advance.
"""
*/
struct LLp{
    Tind row;
    Tval val;
    LLp *next, *reverse;
    LLp():row(0), val(0), next(this), reverse(this){}
    LLp(Tind row, Tval val):row(row), val(val), next(this), reverse(this){}
    LLp(Tind row, Tval val, LLp *next):row(row), val(val), next(next), reverse(this){}
};

/*
"""
  LLmatp is the data structure used to maintain the matrix during elimination.
  It stores the elements in each column in a singly linked list (only next ptrs)
  Each element is an LLp (linked list pointer).
  The head of each column is pointed to by cols.

  We probably can get rid of degs - as it is only used to store initial degrees.
"""
*/
struct LLmatp{
    long n;
    std::vector<Tind> degs;
    std::vector<LLp> cols;
    std::vector<LLp> lles;
    // LLmatp(SparseMatrix a)
};

// these are the types we use with a fixed ordering
struct LLord{
    Tind row;
    Tind next;
    Tval val;
};

struct LLMatOrd{
    long n;
    std::vector<Tind> cols;
    std::vector<LLord> lles;
    // LLmatp(SparseMatrix a){}
    // LLmatp(SparseMatrix a, std::vector<Tind> perm)
};

struct LLcol{
    Tind row;
    Tind ptr;
    Tval cval;
};


// LDLinv

// =============================================================#

// """
//   LDLinv contains the information needed to solve the Laplacian systems.
//   It does it by applying Linv, then Dinv, then Linv (transpose).
//   But, it is specially constructed for this particular solver.
//   It does not explicitly make the matrix triangular.
//   Rather, col[i] is the name of the ith col to be eliminated
// """

struct LDLinv{
    std::vector<Tind> col;
    std::vector<Tind> colptr;
    std::vector<Tind> rowcal;
    std::vector<Tval> fval;
    std::vector<Tval> d;
    LDLinv(LLmatp a):col(std::vector<Tind>(a.n-1, 0)),
                    colptr(std::vector<Tind>(a.n, 0)),
                    rowcal(std::vector<Tind>()),
                    fval(std::vector<Tval>()),
                    d(std::vector<Tval>(a.n, 0)) {}
    LDLinv(LLMatOrd a):col(std::vector<Tind>(a.n-1, 0)),
                    colptr(std::vector<Tind>(a.n, 0)),
                    rowcal(std::vector<Tind>()),
                    fval(std::vector<Tval>()),
                    d(std::vector<Tval>(a.n, 0)) {}
    // LDLinv(SparseMatrix a);
};
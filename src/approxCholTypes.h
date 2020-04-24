#pragma once
#ifndef __CHOL_DATA_STRUCT__
#define __CHOL_DATA_STRUCT__

#include "common.h"
typedef int Tind;
typedef double Tval;

struct Edge{
    Tind r, c;
    Tval w;
    Edge(Tind r, Tind c, Tval w = 1):r(r), c(c), w(w)
    {
        
        // if (r < c)
        // {
        //     int tmp = r;
        //     r = c;
        //     c = tmp;
        // }
    }

    bool operator < (Edge &e)
    {
        return (this->c < e.c) || ((this->c == e.c) && this->r < e.r);
    }
};

bool cmp_edge(const Edge& a, const Edge& b);

//not weighted
struct SparseMatrix{
    std::vector<Tval> vals;
    std::vector<Tind> rows;
    std::vector<Tind> colptrs;
    std::vector<Tval> degs;
    long colnum;
    long elems;
    // long cols;
    // long rows;
    SparseMatrix(std::vector<Edge> edges, long colnum)//suppose input edges to be sorted and x > y
    {
        this->colnum = colnum;
        std::sort(edges.begin(), edges.end(), cmp_edge);
        long len = edges.size();
        vals.resize(len + colnum);
        rows.resize(len + colnum);
        colptrs.resize(colnum + 1);
        degs = std::vector<Tval>(colnum, 0);

        long idx = 0;
        // for (int i = 0; i < colnum; i++)
        // {
        //     while (idx < len && edges[idx].c == i)
        //     {
        //         degs[i] += edges[idx++].w;
        //     }
        // }

        idx = 0;
        long idx_elem = 0;
        for (int i = 0; i < colnum; i++)
        {
            colptrs[i] = idx_elem;
            while (idx < len && edges[idx].c == i)
            {
                vals[idx_elem] = edges[idx].w;
                rows[idx_elem] = edges[idx].r;
                idx++;
                idx_elem++;
            }
            // if (degs[i] > 0)
            // {
            //     vals[idx_elem] = degs[i];
            //     rows[idx_elem] = i;
            //     idx_elem++;
            // }
        }
        colptrs[colnum] = idx_elem;
        this->elems = idx_elem;
    }

    SparseMatrix(std::vector<std::vector<Tval> > &matrix)
    {
        std::vector<Edge> edges;
        int n = matrix.size();
        int m = matrix[0].size();
        std::cout << n << ' ' << m << std::endl;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                if (matrix[i][j])
                    edges.push_back(Edge(i, j, matrix[i][j]));
        new (this)SparseMatrix(edges, m);//placement new
        // std::cout << "in\n" << this->colnum << ' ' << this->elems;
    }
};

std::ostream& operator << (std::ostream &out, SparseMatrix &sparse);

// these are the types we use with a fixed ordering
struct LLord{
    Tind row;
    Tind next;
    Tval val;
    LLord(){}
    LLord(Tind row, Tind next, Tval val):row(row), next(next), val(val){}
};

std::vector<Tind> invperm(std::vector<Tind> &perm);

struct LLMatOrd{
    long n;
    std::vector<Tind> cols;
    std::vector<LLord> lles;
    LLMatOrd(){}
    LLMatOrd(SparseMatrix &a)
    {
        long n = a.colnum;
        this->n = n;
        long m = a.elems;
        // std::cout << n << ' ' <<  m;
        cols.resize(n + 1);

        long ptr = 0;
        for (int i = 0; i < n; i++)
        {
            long next = -1;
            for (int ind = a.colptrs[i]; ind < a.colptrs[i + 1]; ind++)
            {
                Tind j = a.rows[ind];
                if (i < j)
                {
                    Tval v = a.vals[ind];
                    lles.push_back(LLord(j, next, v));
                    next = ptr;
                    ptr++;
                }
            }
            cols[i] = next;
        }
        
    }

    LLMatOrd(SparseMatrix &a, std::vector<Tind> &perm)
    {
        long n = a.colnum;
        this->n = n;
        long m = a.elems;
        cols.resize(n + 1);
        std::vector<Tind> invp = invperm(perm);//why invperm

        long ptr = 0;
        for (int i0 = 0; i0 < n; i0++)
        {
            int i = invp[i0];
            long next = -1;
            for (int ind = a.colptrs[i0]; ind < a.colptrs[i0 + 1]; ind++)
            {
                Tind j = invp[a.rows[ind]];
                if (i < j)
                {
                    Tval v = a.vals[ind];
                    lles.push_back(LLord(j, next, v));
                    next = ptr;
                    ptr++;
                }
            }
            cols[i] = next;
        }
    }
};
std::ostream& operator << (std::ostream &out, LLMatOrd& mat);
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
    std::vector<Tind> rowval;
    std::vector<Tval> fval;
    std::vector<Tval> d;
    LDLinv(LLmatp a):col(std::vector<Tind>(a.n-1, 0)),
                    colptr(std::vector<Tind>(a.n, 0)),
                    rowval(std::vector<Tind>()),
                    fval(std::vector<Tval>()),
                    d(std::vector<Tval>(a.n, 0)) {}
    LDLinv(LLMatOrd a):col(std::vector<Tind>(a.n-1, 0)),
                    colptr(std::vector<Tind>(a.n, 0)),
                    rowval(std::vector<Tind>()),
                    fval(std::vector<Tval>()),
                    d(std::vector<Tval>(a.n, 0)) {}
    // LDLinv(SparseMatrix a);
};

#endif
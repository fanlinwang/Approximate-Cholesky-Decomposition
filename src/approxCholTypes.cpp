#include "approxCholTypes.h"

bool cmp_edge(const Edge& a, const Edge& b){
    return (a.c < b.c) || (a.c == b.c && a.r < b.r);
}

std::ostream& operator << (std::ostream &out, SparseMatrix &sparse)
{
    out << "num of columns:" << sparse.colnum << std::endl;
    out << "num of elems:" << sparse.elems << std::endl;
    out << "colptrs:" << std::endl;
    for (int i = 0; i < sparse.colnum; i++)
    {
        out << sparse.colptrs[i] << ' ';
    }
    out << std::endl;
    out << "(row, val):" << std::endl;
    for (int i = 0; i < sparse.elems; i++)
        out << "(" << sparse.rows[i] << ", " << sparse.vals[i] << ")\n";
    return out;
}

Tind find_col(const SparseMatrix &sparse, Tind i){
    auto pos = std::upper_bound(sparse.colptrs.begin(), 
                                 sparse.colptrs.end(), i);
    if (pos == sparse.colptrs.end()) // last column
        return sparse.colnum-1;
    return std::distance(sparse.colptrs.begin(), pos)-1;
}

std::vector<Tval> operator * (const SparseMatrix& A,
                               const std::vector<Tval>& x)
{
    int n = x.size(), cr = 0;
    std::vector<Tval> b(n, 0.0);
    for (int i = 0; i < A.elems; ++i){
        int cr = A.rows[i];
        Tind col = find_col(A, i);
        b[cr] += A.vals[i] * x[col];        
    }
    return b;
}

std::vector<Tind> invperm(std::vector<Tind> &perm)
{
    std::vector<Tind> invp;
    invp.resize(perm.size());
    for (int i = 0; i < perm.size(); i++)
        invp[perm[i]] = i;
    return invp;
}

std::ostream& operator << (std::ostream &out, LLMatOrd& mat)
{
    out << "num of columns:" << mat.n << std::endl;
    out << "num of elems:" << mat.lles.size() << std::endl;
    out << "colptrs:" << std::endl;
    for (int i = 0; i < mat.n; i++)
    {
        out << mat.cols[i] << ' ';
    }
    out << std::endl;
    out << "(row, next, val):" << std::endl;
    for (int i = 0; i < mat.lles.size(); i++)
        out << "(" << mat.lles[i].row << ", " << mat.lles[i].next << ", " << mat.lles[i].val << ")\n";
    return out;
}
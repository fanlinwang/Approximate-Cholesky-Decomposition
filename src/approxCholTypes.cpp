#include "approxCholTypes.h"

bool cmp_edge(const Edge& a, const Edge& b){
    return (a.c < b.c) || (a.c == b.c && a.r < b.r);
}

std::ostream& operator << (std::ostream &out, SparseMatrix &sparse)
{
    out << "SparseMatrix:" << std::endl;
    out << "num of columns:" << sparse.colnum << std::endl;
    out << "num of elems:" << sparse.elems << std::endl;
    out << "colptrs:" << std::endl;
    for (int i = 0; i < sparse.colnum; i++)
    {
        out << sparse.colptrs[i] << ' ';
    }
    out << std::endl;
    out << "(col, row, val):" << std::endl;
    int col = 0;
    for (int i = 0; i < sparse.elems; i++)
    {
        while (col < sparse.colnum && i >= sparse.colptrs[col + 1])
            col++;
        out << "(" << col << ", " << sparse.rows[i] << ", " << sparse.vals[i] << ")\n";
    }
    return out;
}

Tind find_col(const SparseMatrix &sparse, Tind i){
    auto pos = std::upper_bound(sparse.colptrs.begin(), 
                                 sparse.colptrs.end(), i);
    if (pos == sparse.colptrs.end()) // last column
        return sparse.colnum-1;
    return std::distance(sparse.colptrs.begin(), pos)-1;
}

std::vector<Tval> operator - (const std::vector<Tval>& x,
                               const std::vector<Tval>& y)
{
    assert(x.size() == y.size());
    std::vector<Tval> res(x.size());
    for (int i = 0; i < x.size(); ++i)
        res[i] = x[i] - y[i];
    return res; 
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
    out << "LLMatOrd:" << std::endl;
    out << "num of columns:" << mat.n << std::endl;
    out << "num of elems:" << mat.lles.size() << std::endl;
    out << "colptrs:" << std::endl;
    for (int i = 0; i < mat.n; i++)
    {
        out << mat.cols[i] << ' ';
    }
    out << std::endl;
    out << "(col, row, next, val):" << std::endl;
    int col = 0;
    int i = 0;
    for (int col = 0; col < mat.n; col++)
    {
        while (i <= mat.cols[col])
        {
            out << "(" << col << ", " << mat.lles[i].row << ", " << mat.lles[i].next << ", " << mat.lles[i].val << ")\n";
            i++;
        }
    }
    return out;
}

void laplacian(const SparseMatrix& A, SparseMatrix& L)
{
    std::vector<Edge> edges_list;
    std::vector<Tval> weights_sum(A.colnum, 0);
    int col = 0;
    for (int i = 0; i < A.elems; i++)
    {
        while (col < A.colnum && i >= A.colptrs[col + 1])
            col++;
        edges_list.push_back(Edge(A.rows[i], col, -A.vals[i]));
        weights_sum[col] += A.vals[i];
    }
    for (int i = 0; i < A.colnum; i++)
        edges_list.push_back(Edge(i, i, weights_sum[i]));
    L = SparseMatrix(edges_list, A.colnum);
    // new(&L) SparseMatrix(edges_list, A.colnum);
}

std::ostream& operator << (std::ostream &out, LDLinv &ldli)
{
    out << "col:" << std::endl;
    for (int i = 0; i < ldli.col.size(); i++)
    {
        out << ldli.col[i] << ' ';
    }
    out << std::endl;
    out << "colptr:" << std::endl;
    for (int i = 0; i < ldli.colptr.size(); i++)
    {
        out << ldli.colptr[i] << ' ';
    }
    out << std::endl;
    out << "rowval:" << std::endl;
    for (int i = 0; i < ldli.rowval.size(); i++)
    {
        out << ldli.rowval[i] << ' ';
    }
    out << std::endl;
    out << "fval:" << std::endl;
    for (int i = 0; i < ldli.fval.size(); i++)
    {
        out << ldli.fval[i] << ' ';
    }
    out << std::endl;
    out << "d:" << std::endl;
    for (int i = 0; i < ldli.d.size(); i++)
    {
        out << ldli.d[i] << ' ';
    }
    out << std::endl;
    return out;
}

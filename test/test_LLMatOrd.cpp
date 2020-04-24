#include <common.h>
#include <approxCholTypes.h>

// a=
// [2, 1]  =  1.0
//   [3, 1]  =  1.0
//   [1, 2]  =  1.0
//   [4, 2]  =  1.0
//   [1, 3]  =  1.0
//   [4, 3]  =  1.0
//   [2, 4]  =  1.0
//   [3, 4]  =  1.0
// LLMatOrd(a)
// LLMatOrd{Int64,Float64}(4, [2, 3, 4, 0], Laplacians.LLord{Int64,Float64}[Laplacians.LLord{Int64,Float64}(2, 0, 1.0), Laplacians.LLord{Int64,Float64}(3, 1, 1.0), Laplacians.LLord{Int64,Float64}(4, 0, 1.0), Laplacians.LLord{Int64,Float64}(4, 0, 1.0)
void test_LLMatOrd()
{
    std::vector<std::vector<double> > matrix{
        {0,1,1,0},
        {1,0,0,1},
        {1,0,0,1},
        {0,1,1,0}};

    SparseMatrix sparse(matrix);
    std::cout << sparse;
    LLMatOrd llmat(sparse);
    std::cout << llmat;
}

int main()
{
    test_LLMatOrd();
    return 0;
}
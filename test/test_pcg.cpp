#include <iostream>
#include "approxChol.hpp"
#include "approxCholTypes.h"
#include "test.h"
#include "random"

float get_random()
{
    static std::default_random_engine e;
    static std::normal_distribution<> dis(0, 1); 
    return dis(e);
}

int main(int argc, char* argv[]){

    // std::cout << "construct A" << std::endl;
    const int vertices = atoi(argv[1]);
    const int edges = atoi(argv[2]);
    SparseMatrix A(vertices, edges);
    // std::cout << A;

    std::vector<Tval> b(vertices, 0.0);
    for (auto& e: b)
        e = get_random();
    Tval m = mean(b);
    std::cout << "random b:\n";
    for (auto&e : b){
        e -= m;
        // std::cout << e << " ";
    }
    std::cout << "\ntest sparse matrix multiplication: \n";
    std::vector<Tval> mul = A*b;
    // for (auto m: mul)
    //     std::cout << m << " ";
    // std::cout << "\n";
        
    SparseMatrix lap_A;
    laplacian(A, lap_A);

    SolverParameter para;
    para.maxits = 100;
    para.maxtime = 1000;
    para.verbose = true;
    para.tolerance = 1e-10;
    std::vector<Tval> sol(vertices, 0.0);
    approxchol_lapGiven(A, lap_A, b, sol, para); 

    return 0;
}

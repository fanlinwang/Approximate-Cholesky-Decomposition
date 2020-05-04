// adapted from course website. 

#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "common.h"
#include "tsc_x86.h"
#include "approxCholTypes.h"
#include "approxChol.hpp"

using namespace std;

#define CYCLES_REQUIRED 1e8
#define REP 5
#define MAX_FUNCS 32
// #define FLOPS (4.*n)

typedef LDLinv(*comp_func)(LLMatOrd a);
// typedef void(*comp_func)(SparseMatrix& A, SparseMatrix& lap_A,
//                          const std::vector<Tval>& b,          
//                          std::vector<Tval>& sol, SolverParameter paras);

void add_function(comp_func f, std::string name, int flop);

/* prototype of the function you need to optimize */
double get_perf_score(comp_func f);
void register_functions();
double perf_test(comp_func f, string desc, int flops);
void build_x(double ** m, unsigned n1, unsigned n2);
/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions()
{
    // add_function(&approxchol_lapGiven, "Base iterative solver", 1);
    add_function(&approxChol, "Base approxChol", 1);
}

/* Global vars, used to keep track of student functions */
vector<comp_func> userFuncs;
vector<string> funcNames;
vector<int> funcFlops;
int numFuncs = 0;

float get_random()
{
    static std::random_device e;
    static std::normal_distribution<> dis(-1, 1); 
    return dis(e);
}

void init_vec(std::vector<Tval>& b)
{
    for (auto& e: b)
        e = get_random();
    Tval m = mean(b);
    std::cout << "Generating random vector b...\n";
    for (auto&e : b){
        e -= m;
    }
    std::cout << "b generated.\n";
}

/*
* Compare the solution of the functions 
*/
double nrm_sqr_diff(double *x, double *y, int n) {
    double nrm_sqr = 0.0;
    for(int i = 0; i < n; i++) {
        nrm_sqr += (x[i] - y[i]) * (x[i] - y[i]);
    }
    
    if (isnan(nrm_sqr)) {
      nrm_sqr = INFINITY;
    }
    
    return nrm_sqr;
}

/*
* Registers a user function to be tested by the driver program. Registers a
* string description of the function as well
*/
void add_function(comp_func f, string name, int flops)
{
    userFuncs.push_back(f);
    funcNames.emplace_back(name);
    funcFlops.push_back(flops);

    numFuncs++;
}

/*
* Checks the given function for validity. If valid, then computes and
* reports and returns the number of cycles required per iteration
*/
double perf_test(comp_func f, string desc, int flops)
{
    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    myInt64 start, end;

    // initialize input 
    std::vector<Tval> b(VERTICE, 0.0);
    init_vec(b);

    SparseMatrix A(VERTICE, EDGE);
    LLMatOrd llmat = LLMatOrd(A);
    cout << "Created random sparse matrix A.\n";

    flops = approxChol_count(llmat);

    // SparseMatrix lap_A;
    // laplacian(A, lap_A);

    // SolverParameter para;
    // para.maxits = 100;
    // para.maxtime = 1000;
    // para.verbose = false;
    // para.tolerance = 1e-10;
    // std::vector<Tval> sol(VERTICE, 0.0);

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            // f(A, lap_A, b, sol, para);   
            f(llmat);   
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);
        
    } while (multiplier > 2);

    list<double> cyclesList;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    for (size_t j = 0; j < REP; j++) {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            f(llmat); 
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;

        cyclesList.push_back(cycles);
    }

    cyclesList.sort();
    cycles = cyclesList.front();
    // return cycles; 
    return (1.0 * flops) / cycles;
}

int main(int argc, char **argv) {
    cout << "Starting program. ";
    double perf;
    int i;

    register_functions();

    if (numFuncs == 0){
        cout << endl;
        cout << "No functions registered - nothing for driver to do" << endl;
        cout << "Register functions by calling register_func(f, name)" << endl;
        cout << "in register_funcs()" << endl;

        return 0;
    }
    cout << numFuncs << " functions registered." << endl;

    //Check validity of functions. 
    // build(&x, n);
    // build(&y, n);

    // y_base = static_cast<double *>(malloc(n * sizeof(double)));
    // y_old  = static_cast<double *>(malloc(n * sizeof(double)));

    // memcpy(y_old,  y, n*sizeof(double));
    // // base(); 
    // // TODO: baseline function to be compared
    // memcpy(y_base, y, n*sizeof(double));

    // for (i = 0; i < numFuncs; i++) {
    //     memcpy(y, y_old, n*sizeof(double));
    //     comp_func f = userFuncs[i];
    //     // f(); // TODO
    //     double error = nrm_sqr_diff(y, y_base, n);

    //     if (error > TOLERANCE) {
    //         cout << error << endl;
    //         cout << "ERROR!!!!  the results for the " << i+1 << "th function are different to the previous" << std::endl;
    //     }
    // }

    for (i = 0; i < numFuncs; i++)
    {
        cout << "Starting performance test " << i << "...\n";
        perf = perf_test(userFuncs[i], funcNames[i], 12*EDGE);
        cout << endl << "Running: " << funcNames[i] << endl;
        cout << perf << " flops / cycles" << endl;
    }

    return 0;
}
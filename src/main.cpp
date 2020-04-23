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

using namespace std;

#define CYCLES_REQUIRED 1e7
#define REP 10
#define MAX_FUNCS 32
#define FLOPS (4.*n)
#define EPS (1e-3)

typedef void(*comp_func)(const SparseMatrix& a, const
                         std::vector<Tval>& b, std::vector<Tval>& sol, bool verbose);

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
   // add_function(&approxchol_lapGiven, "Base Solver", 1);
    // add_function(&approxchol_lapGiven1, "Fast Solver", 1);
}

/* Global vars, used to keep track of student functions */
vector<comp_func> userFuncs;
vector<string> funcNames;
vector<int> funcFlops;
int numFuncs = 0;

// TODO: generate a sparse SDDM matrix with dimension n 
//void init_matrices();
void rands(double * m, size_t row, size_t col) 
{
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (size_t i = 0; i < row*col; ++i)  
        m[i] = dist(gen);
}
void build(double **a, int n)
{
    // *a = static_cast<double *>(aligned_alloc(32, n * n * sizeof(double)));
    // rands(*a, n, n);
}

void destroy(double *m)
{
    free(m);
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

    // TODO: initialize input 
    // double *x,*y;
    // int n = NR*NR;

    // build(&x, n);
    // build(&y, n);

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            //f();  //TODO          
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
            // f(); // TODO
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;

        cyclesList.push_back(cycles);
    }

    // destroy(x);
    // destroy(y);
    cyclesList.sort();
    cycles = cyclesList.front();
    return  (1.0 * flops) / cycles;
}

int main(int argc, char **argv) {
    cout << "Starting program. ";
    double perf;
    int i;

    register_functions();
    //init_matrices();

    if (numFuncs == 0){
        cout << endl;
        cout << "No functions registered - nothing for driver to do" << endl;
        cout << "Register functions by calling register_func(f, name)" << endl;
        cout << "in register_funcs()" << endl;

        return 0;
    }
    cout << numFuncs << " functions registered." << endl;

    //Check validity of functions. 
    int n = NR*NR;
    double *x, *y, *y_old, *y_base;
    build(&x, n);
    build(&y, n);
    y_base = static_cast<double *>(malloc(n * sizeof(double)));
    y_old  = static_cast<double *>(malloc(n * sizeof(double)));

    memcpy(y_old,  y, n*sizeof(double));
    // base(); 
    // TODO: baseline function to be compared
    memcpy(y_base, y, n*sizeof(double));

    for (i = 0; i < numFuncs; i++) {
        memcpy(y, y_old, n*sizeof(double));
        comp_func f = userFuncs[i];
        // f(); // TODO
        double error = nrm_sqr_diff(y, y_base, n);

        if (error > EPS) {
            cout << error << endl;
            cout << "ERROR!!!!  the results for the " << i+1 << "th function are different to the previous" << std::endl;
        }
    }
    // destroy(x);
    // destroy(y);
    // destroy(y_base);


    for (i = 0; i < numFuncs; i++)
    {
        perf = perf_test(userFuncs[i], funcNames[i], 24);
        cout << endl << "Running: " << funcNames[i] << endl;
        cout << perf << " flops / cycles" << endl;
    }

    return 0;
}
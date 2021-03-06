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

#define CYCLES_REQUIRED 1e6
#define REP 5
#define MAX_FUNCS 32

// #define FLOPS (4.*n)

int VERTICE;
int EDGE;

typedef LDLinv(*comp_func)(LLMatOrd &a);
typedef LDLinv(*comp_func2)(LLMatOrd_vector2 &a);
typedef LDLinv(*comp_func3)(LLMatOrd_vector2_struct &a);

void add_function(comp_func f, std::string name, int flop);
void add_function(comp_func2 f, std::string name, int flop);
void add_function(comp_func3 f, std::string name, int flop);

/* prototype of the function you need to optimize */
double get_perf_score(comp_func f);
void register_functions();
double perf_test(comp_func f, string desc, int flops);
double perf_test(comp_func2 f, string desc, int flops);
double perf_test(comp_func3 f, string desc, int flops);
void build_x(double ** m, unsigned n1, unsigned n2);
/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions()
{
    // add_function(&approxchol_lapGiven, "Base iterative solver", 1);
    add_function(&approxChol, "Baseline", 1);
    add_function(&approxChol_opt, "Inline", 1);
    add_function(&approxChol_opt2, "Inline+simd", 1);
    add_function(&approxChol_vector2, "Vec", 1);
    add_function(&approxChol_vector2_merge, "VecMg", 1);
    // add_function(&approxChol_vector2_opt, "VecMgSIMD1", 1);
    // add_function(&approxChol_vector2_opt2, "VecMgSIMD2", 1);
    add_function(&approxChol_vector2_opt3, "VecMgSIMD3", 1);
    add_function(&approxChol_vector2_opt4, "VecMgSIMD3csum", 1);
    add_function(&approxChol_vector2_mergerand, "VecMgRand", 1);
    add_function(&approxChol_vector2_mergerand_simd, "VecMgRandSIMD", 1);
    add_function(&approxChol_vector2_merge_search, "VecMgSearch", 1);
    add_function(&approxChol_vector2_merge_search_opt2, "VecMgSearchSIMD", 1);
   // add_function(&approxChol_vector3, "approxChol 3", 1);
   add_function(&approxChol_vector2_struct_merge, "VecStructMg", 1);
   add_function(&approxChol_vector2_struct_merge_simd, "VecStructMgSIMD", 1);
   add_function(&approxChol_vector2_struct_merge_search, "VecStructMgSearch", 1);
   add_function(&approxChol_vector2_struct_merge_search_simd, "VecStructMgSearchSIMD", 1);
   add_function(&approxChol_vector2_struct_merge_search_nofree, "VecStructMgSearch_nofree", 1);

}

/* Global vars, used to keep track of student functions */
vector<comp_func> userFuncs;
vector<comp_func2> userFuncs2;
vector<comp_func3> userFuncs3;
vector<string> funcNames;
vector<string> funcNames2;
vector<string> funcNames3;
vector<int> funcFlops;
int numFuncs = 0;
vector<Tval> b;
SparseMatrix A;
int ops_count = 0;

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

void add_function(comp_func2 f, string name, int flops)
{
    userFuncs2.push_back(f);
    funcNames2.emplace_back(name);
    funcFlops.push_back(flops);

    numFuncs++;
}

void add_function(comp_func3 f, string name, int flops)
{
    userFuncs3.push_back(f);
    funcNames3.emplace_back(name);
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
    long num_runs = 3;
    double multiplier = 1;
    myInt64 start, end;
    LLMatOrd llmat = LLMatOrd(A);
    list<double> cyclesList;

    for (size_t j = 0; j < REP; j++) {
        num_runs = num_runs * multiplier;
        vector<LLMatOrd> llmats(num_runs, llmat);
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            f(llmats[i]); 
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        multiplier = (CYCLES_REQUIRED) / (cycles);
        cyclesList.push_back(cycles);
    }

    cyclesList.sort();
    cycles = cyclesList.front();
    cout << " " << cycles << " ";
    return (1.0 * ops_count) / cycles;
}

double perf_test(comp_func2 f, string desc, int flops)
{
    double cycles = 0.;
    long num_runs = 3;
    double multiplier = 1;
    myInt64 start, end;    

    list<double> cyclesList;
    LLMatOrd_vector2 llmat2(A);

    for (size_t j = 0; j < REP; j++) {
        num_runs = num_runs * multiplier;
        vector<LLMatOrd_vector2> llmats(num_runs, llmat2);
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            f(llmats[i]);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        multiplier = (CYCLES_REQUIRED) / (cycles);
        cyclesList.push_back(cycles);
    }

    cyclesList.sort();
    cycles = cyclesList.front();
    // return cycles; 
    cout << " " << cycles << " ";
    return (1.0 * ops_count) / cycles;
}

double perf_test(comp_func3 f, string desc, int flops)
{
    double cycles = 0.;
    long num_runs = 3;
    double multiplier = 1;
    myInt64 start, end;    

    list<double> cyclesList;
    LLMatOrd_vector2_struct llmat2(A);

    for (size_t j = 0; j < REP; j++) {
        num_runs = num_runs * multiplier;
        vector<LLMatOrd_vector2_struct> llmats(num_runs, llmat2);
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            f(llmats[i]);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;
        multiplier = (CYCLES_REQUIRED) / (cycles);
        cyclesList.push_back(cycles);
    }

    cyclesList.sort();
    cycles = cyclesList.front();
    // return cycles; 
    cout << " " << cycles << " ";
    return (1.0 * ops_count) / cycles;
}

int main(int argc, char **argv) {
    // cout << "Starting program. ";
    double perf;
    int i;
    VERTICE = atoi(argv[1]);
    EDGE = atoi(argv[2]);

    register_functions();

    if (numFuncs == 0){
        cout << endl;
        cout << "No functions registered - nothing for driver to do" << endl;
        cout << "Register functions by calling register_func(f, name)" << endl;
        cout << "in register_funcs()" << endl;

        return 0;
    }
    // cout << numFuncs << " functions registered." << endl;

    // b.reserve(VERTICE);
    // init_vec(b);

    SparseMatrix init_a(VERTICE, EDGE);
    A = init_a;

    LLMatOrd llmat = LLMatOrd(A);

    int flops_count, flcomp_count, intops_count, intcomp_count;
    approxChol_count(llmat, flops_count, flcomp_count, intops_count, intcomp_count);
    ops_count = flops_count + flcomp_count + intops_count + intcomp_count;
    // cout << "#fl ops: " << flops_count << std::endl;
    // cout << "#fl comparisons: " << flcomp_count << std::endl;
    // cout << "#int ops: " << intops_count << std::endl;
    // cout << "#int comparisons: " << intcomp_count << std::endl;
    // cout << "#total ops: " << ops_count << std::endl;

    for (i = 0; i < userFuncs.size(); i++)
    {
        cout << VERTICE << " " << EDGE << " " << flops_count << " " 
            << flcomp_count << " " << intops_count << " " << intcomp_count << " " << ops_count << " ";
        perf = perf_test(userFuncs[i], funcNames[i], 12*EDGE);
        cout << perf << " ";
        cout << funcNames[i] << "\n";
    }
    for (i = 0; i < userFuncs2.size(); i++)
    {
        cout << VERTICE << " " << EDGE << " " << flops_count << " " 
            << flcomp_count << " " << intops_count << " " << intcomp_count << " " << ops_count << " ";
        perf = perf_test(userFuncs2[i], funcNames2[i], 12*EDGE);
        cout << perf << " ";
        cout << funcNames2[i] << "\n";
    }
    for (i = 0; i < userFuncs3.size(); i++)
    {
        cout << VERTICE << " " << EDGE << " " << flops_count << " " 
            << flcomp_count << " " << intops_count << " " << intcomp_count << " " << ops_count << " ";
        perf = perf_test(userFuncs3[i], funcNames3[i], 12*EDGE);
        cout << perf << " ";
        cout << funcNames3[i] << "\n";
    }

    return 0;
}
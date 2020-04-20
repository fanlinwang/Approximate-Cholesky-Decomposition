#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include "approxCholTypes.h"
#define NR 10000

typedef void(*comp_func)(const SparseMatrix& a, const
                         std::vector<Tval>& b, std::vector<Tval>& sol, bool verbose);

void add_function(comp_func f, std::string name, int flop);

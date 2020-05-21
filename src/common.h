#ifndef __COMMON__
#define __COMMON__

#include <vector>
#include <set>
#include <algorithm>
#include <string>
#include <iostream>
#include <cassert>
#include <stdlib.h>
#include <boost/align/aligned_allocator.hpp>
typedef int Tind;
typedef double Tval;

typedef std::vector<Tind, boost::alignment::aligned_allocator<int, 4*sizeof(Tind)>> aligned_VI;
typedef std::vector<Tval, boost::alignment::
    aligned_allocator<Tval, 4*sizeof(Tval)>> aligned_VV;

// #define VERTICE 10000
// #define EDGE 100000


#endif

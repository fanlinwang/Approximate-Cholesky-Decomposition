#include "search.hpp"
#include <iostream>
#include <vector>

using namespace std;

int main()
{
    vector<double> x{0.1, 0.3, 0.6, 1.0};
    vector<double> r{.01, .1,.2,.3,.4,.5,.6,.7,.8,.9,1.0,1.1};
    levelorder_vector<double> x_level(x.begin(), x.end());
    for (int i = 0; i < r.size(); i++)
        cout << r[i] << " "<< x_level.lower_bound(r[i]) << endl;;
    return 0;
}
#include "search.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>

using namespace std;

std::ostream& operator << (std::ostream &out, __m256d d)
{
    out << d[0] << " " << d[1] << " "  << d[2] << " "  << d[3] << endl;
    return out;
}

std::ostream& operator << (std::ostream &out, __m256i d)
{
    out << d[0] << " " << d[1] << " "  << d[2] << " "  << d[3] << endl;
    return out;
}

int main()
{
    std::vector<double> x = {0.1, 0.3, 0.5, 0.8, 1.0, 1.0, 1.0, 1.0};
    std::vector<double> r = {-0.1, 0.2, 0.3, 0.75, 1.1};
    // std::vector<double> x = {1, 2};
    // std::vector<double> r = {1.4};
    for (int i = 0; i < r.size(); i++)
        std::cout << r[i] << ' ' << bitset_search(x, 8, r[i]) << std::endl;

    std::default_random_engine engine;
    std::uniform_real_distribution<double> u(0.0, 1.0);

    int x_len = 80000;
    int r_len = 80000;

    aligned_vector<double> x1(x_len);
    aligned_vector<double> r1(r_len);

    for (int i = 0; i < x_len; i++)
    {
        x1[i] = u(engine);
    }
    sort(x1.begin(), x1.end());

    for (int i = 0; i < r_len; i++)
    {
        r1[i] = u(engine);
    }

    int n = x1.size();
    unsigned nbits = 0;
    while (n >> nbits)
        ++nbits;
    int n_new = (1 << nbits) - 1;
    x1.resize(n_new);
    for (int i = n; i < n_new; i++)
        x1[i] = x1[n - 1];


    // // counting initializing vectors
    // auto start = std::chrono::steady_clock::now();

    // aligned_vector<double> x3(100000);
    // aligned_vector<double> r3(100000);
    // for (int k = 0; k < 100000; k++)
    // {

    //     for (int i = 0; i < r1.size(); i+=4)
    //     {
    //         __m256d rv = _mm256_load_pd(&r1[i]);
    //         __m256i res = bitset_search_simd(x1, n_new, rv);
    //     }
    // }
    // auto end = std::chrono::steady_clock::now();
    // long count1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();



    // vector<double> x(x_len);
    // vector<double> r(r_len);

    // for (int i = 0; i < x_len; i++)
    // {
    //     x[i] = u(engine);
    // }
    // sort(x.begin(), x.end());

    // for (int i = 0; i < r_len; i++)
    // {
    //     r[i] = u(engine);
    // }

    // n = x.size();
    // nbits = 0;
    // while (n >> nbits)
    //     ++nbits;
    // n_new = (1 << nbits) - 1;
    // x.resize(n_new);
    // for (int i = n; i < n_new; i++)
    //     x[i] = x[n - 1];

    // // counting initializing vectors
    // start = std::chrono::steady_clock::now();
    // vector<double> x2(100000);
    // vector<double> r2(100000);
    // for (int k = 0; k < 100000; k++)
    // {

    //     for (int i = 0; i < r.size(); i++)
    //     {
    //         int res = bitset_search(x, n_new, r[i]);
    //     }
    // }
    // end = std::chrono::steady_clock::now();
    // long count2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    // cout << "counting initializing vectors" << endl;
    // cout << "    simd: " << count1 << " ns" << endl;
    // cout << "not simd: " << count2 << " ns" << endl;
    // cout << endl;


    // not counting initializing vectors
    auto start = std::chrono::steady_clock::now();
    // for (int k = 0; k < 100; k++)
    // {
    //     for (int i = 0; i < r1.size(); i+=4)
    //     {
    //         __m256d rv = _mm256_load_pd(&r1[i]);
    //         __m256i res = bitset_search_simd(x1, n_new, rv);
    //     }
    // }
    auto end = std::chrono::steady_clock::now();
    auto count1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    // not counting initializing vectors
    start = std::chrono::steady_clock::now();
    for (int k = 0; k < 100; k++)
    {
        for (int i = 0; i < r1.size(); i++)
        {
            int res = bitset_search(x1, n_new, r1[i]);
        }
    }
    end = std::chrono::steady_clock::now();
    auto count2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    cout << "not counting initializing vectors" << endl;
    // cout << "    simd: " << count1 << " ns" << endl;
    cout << "not simd: " << count2 << " ns" << endl;

    return 0;
}
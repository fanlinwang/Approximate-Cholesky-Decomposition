#pragma once
//modification of http://bannalia.blogspot.com/2015/06/cache-friendly-binary-search.html
#include <algorithm>
#include <vector>
#include <iostream>
#include <immintrin.h>
#include <emmintrin.h>
#include <boost/align/aligned_allocator.hpp>

template <typename T>
using aligned_vector = std::vector<T, boost::alignment::aligned_allocator<T, 32>>;

template<typename T>
class levelorder_vector
{
  typedef std::vector<T> vector;
  
public:
  typedef typename vector::value_type             value_type;
  typedef typename vector::reference              reference;
  typedef typename vector::const_reference        const_reference;
  typedef typename vector::const_iterator         iterator;
  typedef typename vector::const_iterator         const_iterator;
  typedef typename vector::difference_type        difference_type;
  typedef typename vector::size_type              size_type;
  
  levelorder_vector(){}
  levelorder_vector(const levelorder_vector& x):impl(x.impl){}
  levelorder_vector& operator=(const levelorder_vector& x){impl=x.impl;return *this;}
  
  template<typename InputIterator>
  levelorder_vector(InputIterator first,InputIterator last)
  {
    vector aux(first,last);
    // std::sort(aux.begin(),aux.end());
    impl=aux;
    idx.resize(impl.size());
    insert(0,aux.size(),aux.begin());
  }
  
  const_iterator begin()const{return impl.begin();}
  const_iterator end()const{return impl.end();}
  const_iterator cbegin()const{return impl.cbegin();}
  const_iterator cend()const{return impl.cend();}
  friend bool    operator==(const levelorder_vector& x,const levelorder_vector& y)
                   {return x.impl==y.impl;}
  friend bool    operator!=(const levelorder_vector& x,const levelorder_vector& y)
                   {return x.impl!=y.impl;}
  void           swap(levelorder_vector& x){impl.swap(x.impl);}
  friend void    swap(levelorder_vector& x,levelorder_vector& y){x.swap(y);}
  size_type      size()const{return impl.size();}
  size_type      max_size()const{return impl.max_size();}
  bool           empty()const{return impl.empty();}
  
  int lower_bound(const T& x)const
  {
    size_type n=impl.size(),i=n,j=0;
    while(j<n){
    //   if(impl[j]<x){
    //     j=2*j+2;
    //   }
    //   else{
    //     i=j;
    //     j=2*j+1;
    //   }
        bool flag = (impl[j] < x);
        i = flag ? i : j;
        j = flag ? (j << 1) + 2 : (j << 1) + 1;
    }
    return idx[i];
  }
  
private:
  void insert(size_type i,size_type n,const_iterator first)
  {
    if(n){
      size_type h=root(n);
      impl[i]=*(first+h);
      idx[i]=first + h - begin();
      insert(2*i+1,h,first);
      insert(2*i+2,n-h-1,first+h+1);
    }
  }

  size_type root(size_type n)
  {
    if(n<=1)return 0;
    size_type i=2;
    while(i<=n)i*=2;
    return std::min(i/2-1,n-i/4);
  }
public: 
  vector impl;
  std::vector<int> idx;
};


// implementation of bitset/leadbit binary search 
extern "C" __attribute__((optimize("no-tree-vectorize")))
inline int bitset_search(std::vector<double>& x, const int n, double z)
{
  int i = n - 1;
  int k = n;
  int r;
  while(k >>= 1)
  {
    r = i ^ k;
    i = z <= x[r] ? r : i;
  }
  return i;
}

inline int bitset_search(aligned_vector<double>& x, const int n, double z)
{
  int i = n - 1;
  int k = n;
  int r;
  while(k >>= 1)
  {
    r = i ^ k;
    i = z <= x[r] ? r : i;
  }
  return i;
}

inline __m256i bitset_search_simd(aligned_vector<double>& x, const int n, __m256d z)
{
  __m256i i = _mm256_set1_epi64x(n - 1);
  int k = n;
  __m256i r;
  while(k >>= 1)
  {
    // r = i | k;
    r = _mm256_xor_si256(i, _mm256_set1_epi64x(k));
    //load from 4 different places
    // double xr0 = x[r[0]];
    // double xr1 = x[r[1]];
    // double xr2 = x[r[2]];
    // double xr3 = x[r[3]];
    // __m256d xr = _mm256_set_pd(xr3, xr2, xr1, xr0);
    __m256d xr = _mm256_i64gather_pd(&x[0], r, 8);
    __m256d mask = _mm256_cmp_pd(z, xr, _CMP_LE_OQ);
    i = _mm256_blendv_epi8(i, r, _mm256_castpd_si256(mask));
    // i = z >= x[r] ? r : i;
  }
  return i;
}
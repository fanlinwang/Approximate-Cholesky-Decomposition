//modification of http://bannalia.blogspot.com/2015/06/cache-friendly-binary-search.html
#include <algorithm>
#include <vector>
#include <iostream>

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
        i = flag ? j : i;
        j = flag ? (j << 1) + 2 : (j << 1) + 1;
    }
    return i;
  }
  
private:
  void insert(size_type i,size_type n,const_iterator first)
  {
    if(n){
      size_type h=root(n);
      impl[i]=*(first+h);
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
  
  vector impl;
};
#include "approxChol.hpp"
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <emmintrin.h>
#include <iostream>
#include <random>
#include <cmath>
#include <ctime>
#include <limits>

namespace std {

namespace sort_helper {

template <typename _Data, typename _Order>
struct value_reference_t;

template <typename _Data, typename _Order>
struct value_t {
    _Data data;
    _Order val;
    inline value_t(_Data _data, _Order _val) : data(_data), val(_val) {}
    inline value_t(const value_reference_t<_Data,_Order>& rhs);
};

template <typename _Data, typename _Order>
struct value_reference_t {
    _Data* pdata;
    _Order* pval;
    value_reference_t(_Data* _itData, _Order* _itVal) : pdata(_itData), pval(_itVal) {}
    inline value_reference_t& operator = (const value_reference_t& rhs) { *pdata = *rhs.pdata; *pval = *rhs.pval; return *this; }
    inline value_reference_t& operator = (const value_t<_Data,_Order>& rhs) { *pdata = rhs.data; *pval = rhs.val; return *this; }
    inline bool operator < (const value_reference_t& rhs) { return *pval < *rhs.pval; }
};

template <typename _Data, typename _Order>
struct value_iterator_t :
    iterator< random_access_iterator_tag, value_t<_Data,_Order>, ptrdiff_t, value_t<_Data,_Order>*, value_reference_t<_Data,_Order> >
{
    _Data* itData;
    _Order* itVal;
    value_iterator_t(_Data* _itData, _Order* _itVal) : itData(_itData), itVal(_itVal) {}
    inline ptrdiff_t operator - (const value_iterator_t& rhs) const { return itVal - rhs.itVal; }
    inline value_iterator_t operator + (ptrdiff_t off) const { return value_iterator_t(itData + off, itVal + off); }
    inline value_iterator_t operator - (ptrdiff_t off) const { return value_iterator_t(itData - off, itVal - off); }
    inline value_iterator_t& operator ++ () { ++itData; ++itVal; return *this; }
    inline value_iterator_t& operator -- () { --itData; --itVal; return *this; }
    inline value_iterator_t operator ++ (int) { return value_iterator_t(itData++, itVal++); }
    inline value_iterator_t operator -- (int) { return value_iterator_t(itData--, itVal--); }
    inline value_t<_Data,_Order> operator * () const { return value_t<_Data,_Order>(*itData, *itVal); }
    inline value_reference_t<_Data,_Order> operator * () { return value_reference_t<_Data,_Order>(itData, itVal); }
    inline bool operator  < (const value_iterator_t& rhs) const { return itVal  < rhs.itVal; }
    inline bool operator == (const value_iterator_t& rhs) const { return itVal == rhs.itVal; }
    inline bool operator != (const value_iterator_t& rhs) const { return itVal != rhs.itVal; }
};

template <typename _Data, typename _Order>
inline value_t<_Data,_Order>::value_t(const value_reference_t<_Data,_Order>& rhs)
    : data(*rhs.pdata), val(*rhs.pval) {}

template <typename _Data, typename _Order>
bool operator < (const value_t<_Data,_Order>& lhs, const value_reference_t<_Data,_Order>& rhs) {
    return lhs.val < *rhs.pval; }

template <typename _Data, typename _Order>
bool operator < (const value_reference_t<_Data,_Order>& lhs, const value_t<_Data,_Order>& rhs) {
    return *lhs.pval < rhs.val; }

template <typename _Data, typename _Order>
void swap(value_reference_t<_Data,_Order> lhs, value_reference_t<_Data,_Order> rhs) {
    std::swap(*lhs.pdata, *rhs.pdata);
    std::swap(*lhs.pval, *rhs.pval); }


} // namespace sort_helper

} // namespace std

void print_ll_col(LLMatOrd llmat, int i) {
    Tind ptr = llmat.cols[i];
    while (ptr != -1)
    {
        LLord ll = llmat.lles[ptr];
        std::cout << "col " << i << ", row " << ll.row << " : " << ll.val << std::endl;

        ptr = ll.next;
    }
}

int get_ll_col(LLMatOrd llmat, int i, std::vector<LLcol> &colspace) {
    
    Tind ptr = llmat.cols[i];
    int len = 0;

    while (ptr != -1) {
        LLcol item = {llmat.lles[ptr].row, ptr, llmat.lles[ptr].val};

        if (len >= colspace.size()) {
            colspace.push_back(item);
        } else {
            colspace[len] = item;
        }

        len += 1;
        ptr = llmat.lles[ptr].next;
    }

    return len;
}

bool cmp_row(const LLcol &a, const LLcol &b) {
    return a.row < b.row;
}

bool cmp_val(const LLcol &a, const LLcol &b) {
    return a.cval < b.cval;
}

Tind compressCol(std::vector<LLcol> &colspace, int len) {

    // sort colspace ? DY: just checked, yes! 
    std::sort(colspace.begin(), colspace.begin() + len, cmp_row);

    std::vector<LLcol> c = colspace;

    Tind ptr = -1;
    Tind currow = c[0].row;     // julia index start from 1? 
    Tval curval = c[0].cval;
    Tind curptr = c[0].ptr;

    // not consider inbounds
    for (int i = 1; i < len; i++) {

        if (c[i].row != currow) {

            ptr += 1;
            c[ptr].row = currow;
            c[ptr].cval = curval;
            c[ptr].ptr = curptr;

            currow = c[i].row;
            curval = c[i].cval;
            curptr = c[i].ptr;

        } else {
            
            curval += c[i].cval;

        }
    }

    // emit the last row
    ptr += 1;
    c[ptr].row = currow;
    c[ptr].cval = curval;
    c[ptr].ptr = curptr;

    std::sort(c.begin(), c.begin() + ptr + 1, cmp_val);
    
    colspace = c;
    return ptr+1; // if the returned value is used for vector length, then it should return ptr+1
}

LDLinv approxChol(LLMatOrd a) {
    auto n = a.n;

    // need to make custom one without col info later
    LDLinv ldli(a);
    Tind ldli_row_ptr = 0;

    std::vector<Tval> d(n, 0);

    std::vector<LLcol> colspace(n);
    std::vector<Tval> cumspace(n);

    // random engine and distribution
    std::default_random_engine engine;
    std::uniform_real_distribution<Tval> u(0.0, 1.0);


    for (long i = 0; i <= n-2; i++) {

        ldli.col[i] = i;
        ldli.colptr[i] = ldli_row_ptr;

        int len = get_ll_col(a, i, colspace);
        len = compressCol(colspace, len);
        // flop count: sort LLcol len*log(len)?

        Tval csum = 0;
        for (int ii = 0; ii < len; ii++) {
            csum += colspace[ii].cval;
            // flop count: 1 add
            cumspace[ii] = csum;
        }
        Tval wdeg = csum;

        Tval colScale = 1;

        Tind ks[len-1];
        Tval randnums[len-1];
        for (int joffset = 0; joffset <= len-2; joffset++)
        {
            randnums[joffset] = u(engine);
        }

        for (int joffset = 0; joffset <= len-2; joffset++) {
            Tval r = randnums[joffset];
            r = r * (csum - cumspace[joffset]) + cumspace[joffset];

            auto cumspace_last = cumspace.begin();
            std::advance(cumspace_last, len);
            int koff = std::distance(cumspace.begin(), std::lower_bound(cumspace.begin(), cumspace_last, r));

            ks[joffset] = colspace[koff].row;
        }

        Tval newEdgeVals[len-1];
        Tind js[len-1];
        Tind ptrs[len-1];
        for (int joffset = 0; joffset <= len-2; joffset++) {
            LLcol llcol = colspace[joffset];
            Tval w = llcol.cval * colScale;
            js[joffset] = llcol.row;
            ptrs[joffset] = llcol.ptr;

            Tval f = w/wdeg;
            // flop count: 1 mul

            newEdgeVals[joffset] = w*(1-f);

            colScale *= 1 - f;
            wdeg = wdeg - 2*w + w*f; // wdeg *= (1-f)^2
            // flop count: 3 mul 3 add

            ldli.rowval.push_back(js[joffset]);
            ldli.fval.push_back(f);
            ldli_row_ptr += 1;
        }

        for (int joffset = 0; joffset <= len-2; joffset++) {
            // create edge (j,k) with newEdgeVal
            // do it by reassigning ll
            Tind j = js[joffset];
            Tind k = ks[joffset];
            Tval newEdgeVal = newEdgeVals[joffset];
            Tind ptr = ptrs[joffset];
            if (j < k) {    // put it in col j
                Tind jhead = a.cols[j];
                a.lles[ptr].row = k;
                a.lles[ptr].next = jhead;
                a.lles[ptr].val = newEdgeVal;
                a.cols[j] = ptr;
            } else {        // put it in col k
                Tind khead = a.cols[k];
                a.lles[ptr].row = j;
                a.lles[ptr].next = khead;
                a.lles[ptr].val = newEdgeVal;
                a.cols[k] = ptr;
            }
        }

        LLcol llcol = colspace[len-1];
        Tval w = llcol.cval * colScale;
        // flop count: 1 mul
        Tind j = llcol.row;

        ldli.rowval.push_back(j);
        ldli.fval.push_back(1);
        ldli_row_ptr += 1;

        d[i] = w;
    }

    ldli.colptr[n-1] = ldli_row_ptr;
    ldli.d = d;

    return ldli;
}

LDLinv approxChol_vector2(LLMatOrd_vector2 a) {
    auto n = a.n;

    // need to make custom one without col info later
    LDLinv ldli(a);
    Tind ldli_row_ptr = 0;

    std::vector<Tval> d(n, 0);

    // std::vector<LLcol> colspace(n);


    std::vector<Tval> cumspace(n);

    // random engine and distribution
    std::default_random_engine engine;
    std::uniform_real_distribution<Tval> u(0.0, 1.0);


    for (long i = 0; i <= n-2; i++) {

        ldli.col[i] = i;
        ldli.colptr[i] = ldli_row_ptr;

        int len = a.row[i].size();
        // TODO: is sort by val necessary
        typedef std::sort_helper::value_iterator_t<Tind,Tval> IndexIt;
        std::sort(IndexIt(&a.row[i][0], &a.val[i][0]), IndexIt(&a.row[i][0] + len, &a.val[i][0] + len) );
        // flop count: sort LLcol len*log(len)?

        Tval csum = 0;
        for (int ii = 0; ii < len; ii++) {
            csum += a.val[i][ii];
            // flop count: 1 add
            cumspace[ii] = csum;
        }
        Tval wdeg = csum;

        Tval colScale = 1;

        Tind ks[len-1];
        Tval randnums[len-1];
        for (int joffset = 0; joffset <= len-2; joffset++)
        {
            randnums[joffset] = u(engine);
        }

        for (int joffset = 0; joffset <= len-2; joffset++) {
            Tval r = randnums[joffset];
            r = r * (csum - cumspace[joffset]) + cumspace[joffset];

            auto cumspace_last = cumspace.begin();
            std::advance(cumspace_last, len);
            int koff = std::distance(cumspace.begin(), std::lower_bound(cumspace.begin(), cumspace_last, r));

            ks[joffset] = a.row[i][koff];
        }

        Tval newEdgeVals[len-1];
        Tind js[len-1];
        Tind ptrs[len-1];
        for (int joffset = 0; joffset <= len-2; joffset++) {
            // LLcol llcol = colspace[joffset];
            Tval w = a.val[i][joffset] * colScale;
            js[joffset] = a.row[i][joffset];
            // ptrs[joffset] = llcol.ptr;

            Tval f = w/wdeg;
            // flop count: 1 mul

            newEdgeVals[joffset] = w*(1-f);

            colScale *= 1 - f;
            wdeg = wdeg - 2*w + w*f; // wdeg *= (1-f)^2
            // flop count: 3 mul 3 add

            ldli.rowval.push_back(js[joffset]);
            ldli.fval.push_back(f);
            ldli_row_ptr += 1;
        }

        for (int joffset = 0; joffset <= len-2; joffset++) {
            // create edge (j,k) with newEdgeVal
            // do it by reassigning ll
            Tind j = js[joffset];
            Tind k = ks[joffset];
            Tval newEdgeVal = newEdgeVals[joffset];
            Tind ptr = ptrs[joffset];
            if (j < k) {    // put it in col j
                // Tind jhead = a.cols[j];
                // a.lles[ptr].row = k;
                // a.lles[ptr].next = jhead;
                // a.lles[ptr].val = newEdgeVal;
                // a.cols[j] = ptr;
                // std::cout << "start create j < k" << j << " " << k << std::endl;
                auto iter = std::lower_bound(a.row[j].begin(), a.row[j].end(), k);
                int idx = iter - a.row[j].begin();
                // std::cout << "iter" << int(iter - a.row[j].begin()) << iter - a.row[j].begin() <<  std::endl;
                if (iter != a.row[j].end() && *iter == k)
                {
                    // std::cout << "exist" << std::endl;
                    a.val[j][idx] += newEdgeVal;
                }
                else
                {
                    // std::cout << "new" << std::endl;
                    a.row[j].insert(iter, k);
                    // std::cout << "new row finish" << std::endl;
                    a.val[j].insert(a.val[j].begin() + idx, newEdgeVal);
                    // std::cout << "new val finish" << std::endl;
                } 
                // std::cout << "created j < k" << std::endl;
            } else {        // put it in col k
                auto iter = std::lower_bound(a.row[k].begin(), a.row[k].end(), j);
                int idx = iter - a.row[k].begin();
                if (iter != a.row[k].end() && *iter == j)
                    a.val[k][idx] += newEdgeVal;
                else
                {
                    a.row[k].insert(iter, j);
                    a.val[k].insert(a.val[k].begin() + idx, newEdgeVal);
                } 
            }
        }

        // LLcol llcol = colspace[len-1];
        Tval w = a.val[i][len - 1] * colScale;
        // flop count: 1 mul
        Tind j = a.row[i][len - 1];

        ldli.rowval.push_back(j);
        ldli.fval.push_back(1);
        ldli_row_ptr += 1;

        d[i] = w;

        //free column
        a.row[i].clear();
        a.row[i].shrink_to_fit();
        a.val[i].clear();
        a.val[i].shrink_to_fit();
    }

    ldli.colptr[n-1] = ldli_row_ptr;
    ldli.d = d;

    return ldli;
}

LDLinv approxChol_vector2_opt(LLMatOrd_vector2 a) {
    auto n = a.n;

    LDLinv ldli(a);
    Tind ldli_row_ptr = 0;

    std::vector<Tval> d(n, 0);


    std::vector<Tval> cumspace(n);

    std::default_random_engine engine;
    std::uniform_real_distribution<Tval> u(0.0, 1.0);


    for (long i = 0; i <= n-2; i++) {

        ldli.col[i] = i;
        ldli.colptr[i] = ldli_row_ptr;

        int len = a.row[i].size();
        typedef std::sort_helper::value_iterator_t<Tind,Tval> IndexIt;
        std::sort(IndexIt(&a.row[i][0], &a.val[i][0]), IndexIt(&a.row[i][0] + len, &a.val[i][0] + len) );
        // flop count: sort LLcol len*log(len)?

        Tval csum = 0;
        for (int ii = 0; ii < len; ii++) {
            csum += a.val[i][ii];
            // flop count: 1 add
            cumspace[ii] = csum;
        }
        Tval wdeg = csum;

        Tval colScale = 1;

        int newlen = ceil(double(len-1)/4)*4;
        Tind js[newlen];
        Tind ks[newlen];
        
        /*Tind *js, *ks;
        int newlen = ceil(double(len-1)/4)*4;
        js = static_cast<Tind *>(aligned_alloc(4*sizeof(Tind), newlen * sizeof(Tind)));
        ks = static_cast<Tind *>(aligned_alloc(4*sizeof(Tind), newlen * sizeof(Tind)));*/
        
        Tval randnums[len-1];
        for (int joffset = 0; joffset <= len-2; joffset++)
        {
            randnums[joffset] = u(engine);
        }

        for (int joffset = 0; joffset <= len-2; joffset++) {
            Tval r = randnums[joffset];
            r = r * (csum - cumspace[joffset]) + cumspace[joffset];

            auto cumspace_last = cumspace.begin();
            std::advance(cumspace_last, len);
            int koff = std::distance(cumspace.begin(), std::lower_bound(cumspace.begin(), cumspace_last, r));

            ks[joffset] = a.row[i][koff];
        }

        Tval newEdgeVals[len-1];
        Tind ptrs[len-1];
        for (int joffset = 0; joffset <= len-2; joffset++) {
            // LLcol llcol = colspace[joffset];
            Tval w = a.val[i][joffset] * colScale;
            js[joffset] = a.row[i][joffset];
            // ptrs[joffset] = llcol.ptr;

            Tval f = w/wdeg;
            // flop count: 1 mul

            newEdgeVals[joffset] = w*(1-f);

            colScale *= 1 - f;
            wdeg = wdeg - 2*w + w*f; // wdeg *= (1-f)^2
            // flop count: 3 mul 3 add

            ldli.rowval.push_back(js[joffset]);
            ldli.fval.push_back(f);
            ldli_row_ptr += 1;
        }
        
        __m128i allone = _mm_set1_epi32(0xffffffff);
        for (int joffset = 0; joffset < newlen; joffset+=4) {
            __m128i j_4 = _mm_loadu_si128((__m128i *)(js+joffset));
            __m128i k_4 = _mm_loadu_si128((__m128i *)(ks+joffset));

            __m128i mask1 = _mm_cmplt_epi32(k_4, j_4);
            __m128i res1 = _mm_blendv_epi8(j_4, k_4, mask1);

            //__m128i mask2 = _mm_cmplt_epi32(j_4, k_4);
            __m128i mask2 = _mm_xor_si128(mask1, allone);
            __m128i res2 = _mm_blendv_epi8(j_4, k_4, mask2);

            //__m128i res1 = _mm_min_epi32(j_4, k_4);
            //__m128i res2 = _mm_max_epi32(j_4, k_4);

            _mm_storeu_si128((__m128i *)(js+joffset), res1);
            _mm_storeu_si128((__m128i *)(ks+joffset), res2);
        }
        
        for (int joffset = 0; joffset <= len-2; joffset++) {
            // create edge (j,k) with newEdgeVal
            // do it by reassigning ll
            Tind j = js[joffset];
            Tind k = ks[joffset];
            Tval newEdgeVal = newEdgeVals[joffset];
            Tind ptr = ptrs[joffset];

            auto iter = std::lower_bound(a.row[j].begin(), a.row[j].end(), k);
            int idx = iter - a.row[j].begin();
            if (iter != a.row[j].end() && *iter == k)
            {
                a.val[j][idx] += newEdgeVal;
            }
            else
            {
                a.row[j].insert(iter, k);
                a.val[j].insert(a.val[j].begin() + idx, newEdgeVal);
            }
        }
        
        //free(js);
        //free(ks);
        // LLcol llcol = colspace[len-1];
        Tval w = a.val[i][len - 1] * colScale;
        // flop count: 1 mul
        Tind j = a.row[i][len - 1];

        ldli.rowval.push_back(j);
        ldli.fval.push_back(1);
        ldli_row_ptr += 1;

        d[i] = w;

        //free column
        a.row[i].clear();
        a.row[i].shrink_to_fit();
        a.val[i].clear();
        a.val[i].shrink_to_fit();
    }

    ldli.colptr[n-1] = ldli_row_ptr;
    ldli.d = d;

    return ldli;
}

LDLinv approxChol_opt(LLMatOrd a) {
    auto n = a.n;

    // need to make custom one without col info later
    LDLinv ldli(a);
    Tind ldli_row_ptr = 0;

    std::vector<Tval> d(n, 0);

    std::vector<LLcol> colspace(n);
    std::vector<Tval> cumspace(n);

    // random engine and distribution
    std::default_random_engine engine;
    std::uniform_real_distribution<Tval> u(0.0, 1.0);


    for (long i = 0; i <= n-2; i++) {

        ldli.col[i] = i;
        ldli.colptr[i] = ldli_row_ptr;

        //int len = get_ll_col(a, i, colspace);
        Tind ptr = a.cols[i];
        int len = 0;

        while (ptr != -1) {
            LLcol item = {a.lles[ptr].row, ptr, a.lles[ptr].val};

            if (len >= colspace.size()) {
                colspace.push_back(item);
            } else {
                colspace[len] = item;
            }

            len += 1;
            ptr = a.lles[ptr].next;
        }

        //len = compressCol(colspace, len);
        std::sort(colspace.begin(), colspace.begin() + len, cmp_row);
        std::vector<LLcol> c = colspace;
        ptr = -1;
        Tind currow = c[0].row;     // julia index start from 1? 
        Tval curval = c[0].cval;
        Tind curptr = c[0].ptr;

        for (int i = 1; i < len; i++) {

            if (c[i].row != currow) {

                ptr += 1;
                c[ptr].row = currow;
                c[ptr].cval = curval;
                c[ptr].ptr = curptr;

                currow = c[i].row;
                curval = c[i].cval;
                curptr = c[i].ptr;

            } else {
                
                curval += c[i].cval;

            }
        }

        ptr += 1;
        c[ptr].row = currow;
        c[ptr].cval = curval;
        c[ptr].ptr = curptr;

        std::sort(c.begin(), c.begin() + ptr + 1, cmp_val);
        
        colspace = c;
        len = ptr+1;



        Tval csum = 0;
        for (int ii = 0; ii < len; ii++) {
            csum += colspace[ii].cval;
            // flop count: 1 add
            cumspace[ii] = csum;
        }
        Tval wdeg = csum;

        Tval colScale = 1;

        Tind ks[len-1];
        Tval randnums[len-1];
        int joffset;
        for (joffset = 0; joffset <= len-2; joffset++)
        {
            randnums[joffset] = u(engine);
        }

        for (joffset = 0; joffset <= len-2; joffset++) {
            Tval r = randnums[joffset];
            r = r * (csum - cumspace[joffset]) + cumspace[joffset];

            auto cumspace_last = cumspace.begin();
            std::advance(cumspace_last, len);
            int koff = std::distance(cumspace.begin(), std::lower_bound(cumspace.begin(), cumspace_last, r));

            ks[joffset] = colspace[koff].row;
        }

        Tval newEdgeVals[len-1];
        Tind js[len-1];
        Tind ptrs[len-1];
        for (joffset = 0; joffset <= len-3; joffset+=2) {
            LLcol llcol1 = colspace[joffset];
            LLcol llcol2 = colspace[joffset+1];

            Tval w1 = llcol1.cval * colScale;
            Tval f1 = 1-w1/wdeg;
            wdeg *= f1*f1;

            Tval w2 = llcol2.cval * colScale * f1;
            Tval f2 = 1 - w2/wdeg;
            wdeg *= f2*f2;
            colScale *= f1*f2;

            js[joffset] = llcol1.row;
            ptrs[joffset] = llcol1.ptr;
            js[joffset+1] = llcol2.row;
            ptrs[joffset+1] = llcol2.ptr;


            newEdgeVals[joffset] = w1*f1;
            newEdgeVals[joffset+1] = w2*f2;


            ldli.rowval.push_back(js[joffset]);
            ldli.fval.push_back(1-f1);
            ldli.rowval.push_back(js[joffset+1]);
            ldli.fval.push_back(1-f2);
        }
        for (; joffset <= len-2; joffset++)
        {
            LLcol llcol = colspace[joffset];
            Tval w = llcol.cval * colScale;
            js[joffset] = llcol.row;
            ptrs[joffset] = llcol.ptr;

            Tval f = w/wdeg;
            // flop count: 1 mul

            newEdgeVals[joffset] = w*(1-f);

            colScale *= 1 - f;
            wdeg = wdeg - 2*w + w*f; // wdeg *= (1-f)^2
            // flop count: 3 mul 3 add

            ldli.rowval.push_back(js[joffset]);
            ldli.fval.push_back(f);
        }

        for (joffset = 0; joffset <= len-2; joffset++) {
            // create edge (j,k) with newEdgeVal
            // do it by reassigning ll
            Tind j = js[joffset];
            Tind k = ks[joffset];
            Tval newEdgeVal = newEdgeVals[joffset];
            Tind ptr = ptrs[joffset];
            if (j < k) {    // put it in col j
                Tind jhead = a.cols[j];
                a.lles[ptr].row = k;
                a.lles[ptr].next = jhead;
                a.lles[ptr].val = newEdgeVal;
                a.cols[j] = ptr;
            } else {        // put it in col k
                Tind khead = a.cols[k];
                a.lles[ptr].row = j;
                a.lles[ptr].next = khead;
                a.lles[ptr].val = newEdgeVal;
                a.cols[k] = ptr;
            }
        }

        LLcol llcol = colspace[len-1];
        Tval w = llcol.cval * colScale;
        // flop count: 1 mul
        Tind j = llcol.row;

        ldli.rowval.push_back(j);
        ldli.fval.push_back(1);
        ldli_row_ptr += len;

        d[i] = w;
    }

    ldli.colptr[n-1] = ldli_row_ptr;
    ldli.d = d;

    return ldli;
}

int get_ll_col_count(LLMatOrd llmat, int i, std::vector<LLcol> &colspace, int &intops, int &intcomp) {
    
    Tind ptr = llmat.cols[i];
    int len = 0;

    while (ptr != -1) {
        LLcol item = {llmat.lles[ptr].row, ptr, llmat.lles[ptr].val};

        if (len >= colspace.size()) {
            colspace.push_back(item);
        } else {
            colspace[len] = item;
        }

        len += 1;
        ptr = llmat.lles[ptr].next;
        intops++;
        intcomp++;
    }

    return len;
}

int getPartition_row(std::vector<LLcol> &nums, int low, int height, int &intops, int &intcomp)
{
    LLcol keyVal = nums[low];
    while (low<height)
    {
        intcomp++;
        while (low < height && nums[height].row >= keyVal.row) {
            height--;
            intcomp += 2;
            intops++;
        }
        intcomp += 2;
        nums[low] = nums[height];
        while (low < height && nums[low].row <= keyVal.row) {
            low++;
            intcomp += 2;
            intops++;
        }
        intcomp += 2;
        nums[height] = nums[low];
    }
    intcomp++;
    nums[low] = keyVal;
    return low;
}

void sort_row(std::vector<LLcol> &c, int low, int high, int &intops, int &intcomp) {

    if (low < high)
    {
        int mid = getPartition_row(c, low, high, intops, intcomp);
        sort_row(c, low, mid-1, intops, intcomp);
        sort_row(c, mid + 1, high, intops, intcomp);
    }
    intcomp++;
}

int getPartition_val(std::vector<LLcol> &nums, int low, int height, int &intops, int &intcomp, int &flcomp, int &flops)
{
    LLcol keyVal = nums[low];
    while (low<height)
    {
        intcomp++;
        while (low < height && nums[height].cval >= keyVal.cval) {
            height--;
            intcomp++;
            intops++;
            flcomp++;
        }
        intcomp++;
        flcomp++;
        nums[low] = nums[height];
        while (low < height && nums[low].cval <= keyVal.cval) {
            low++;
            intcomp++;
            intops++;
            flcomp++;
        }
        intcomp++;
        flcomp++;
        nums[height] = nums[low];
    }
    intcomp++;
    nums[low] = keyVal;
    return low;
}

void sort_val(std::vector<LLcol> &c, int low, int high, int &intops, int &intcomp, int &flcomp, int &flops) {

    if (low < high)
    {
        int mid = getPartition_val(c, low, high, intops, intcomp, flcomp, flops);
        sort_val(c, low, mid-1, intops, intcomp, flcomp, flops);
        sort_val(c, mid + 1, high, intops, intcomp, flcomp, flops);
    }
    intcomp++;
}

Tind compressCol_count(std::vector<LLcol> &colspace, int len, int &intops, int &intcomp, int &flcomp, int &flops) {

    // sort colspace ? DY: just checked, yes! 
    sort_row(colspace, 0, len-1, intops, intcomp);

    std::vector<LLcol> c = colspace;

    Tind ptr = -1;
    Tind currow = c[0].row;     // julia index start from 1? 
    Tval curval = c[0].cval;
    Tind curptr = c[0].ptr;
    intops++;

    // not consider inbounds
    for (int i = 1; i < len; i++) {

        if (c[i].row != currow) {

            ptr += 1;
            c[ptr].row = currow;
            c[ptr].cval = curval;
            c[ptr].ptr = curptr;

            currow = c[i].row;
            curval = c[i].cval;
            curptr = c[i].ptr;
            intops++;

        } else {
            
            curval += c[i].cval;
            flops++;

        }
        intcomp++;
    }

    // emit the last row
    ptr += 1;
    c[ptr].row = currow;
    c[ptr].cval = curval;
    c[ptr].ptr = curptr;
    intops++;

    sort_val(colspace, 0, len-1, intops, intcomp, flcomp, flops);
    
    colspace = c;
    return ptr+1; // if the returned value is used for vector length, then it should return ptr+1
}

void approxChol_count(LLMatOrd a, int &flops_count, int &flcomp_count, int &intops_count, int &intcomp_count) {
    auto n = a.n;
    flops_count = 0;
    flcomp_count = 0;
    intops_count = 0;
    intcomp_count = 0;

    // need to make custom one without col info later
    LDLinv ldli(a);
    Tind ldli_row_ptr = 0;

    std::vector<Tval> d(n, 0);

    std::vector<LLcol> colspace(n);
    std::vector<Tval> cumspace(n);

    // random engine and distribution
    std::default_random_engine engine;
    std::uniform_real_distribution<Tval> u(0.0, 1.0);

    for (long i = 0; i <= n-2; i++) {

        ldli.col[i] = i;
        ldli.colptr[i] = ldli_row_ptr;

        int len = get_ll_col_count(a, i, colspace, intops_count, intcomp_count);
        len = compressCol_count(colspace, len, intops_count, intcomp_count, flcomp_count, flops_count);
        // flop count: sort LLcol len*log(len)?

        Tval csum = 0;
        for (int ii = 0; ii < len; ii++) {
            csum += colspace[ii].cval;
            flops_count++;
            cumspace[ii] = csum;
        }
        Tval wdeg = csum;

        Tval colScale = 1;

        for (int joffset = 0; joffset <= len-2; joffset++) {

            LLcol llcol = colspace[joffset];
            Tval w = llcol.cval * colScale;
            Tind j = llcol.row;

            Tval f = w/wdeg;
            flops_count++;
            // flop count: 1 mul

            Tval r = u(engine);
            //Tval r = randnums[randptr++];
            r = r * (csum - cumspace[joffset]) + cumspace[joffset];
            flops_count += 3;
            // flop count: 1 mul 2 add

            auto cumspace_last = cumspace.begin();
            std::advance(cumspace_last, len);
            int koff = std::distance(cumspace.begin(), std::lower_bound(cumspace.begin(), cumspace_last, r));
            // flop count: len?

            Tind k = colspace[koff].row;

            Tval newEdgeVal = w*(1-f);
            flops_count += 2;
            // flop count: 1 mul 1 add

            // create edge (j,k) with newEdgeVal
            // do it by reassigning ll
            if (j < k) {    // put it in col j
                Tind jhead = a.cols[j];
                a.lles[llcol.ptr].row = k;
                a.lles[llcol.ptr].next = jhead;
                a.lles[llcol.ptr].val = newEdgeVal;
                a.cols[j] = llcol.ptr;
            } else {        // put it in col k
                Tind khead = a.cols[k];
                a.lles[llcol.ptr].row = j;
                a.lles[llcol.ptr].next = khead;
                a.lles[llcol.ptr].val = newEdgeVal;
                a.cols[k] = llcol.ptr;
            }
            intcomp_count++;

            colScale *= 1 - f;
            wdeg = wdeg - 2*w + w*f;
            flops_count += 6;
            // flop count: 3 mul 3 add

            ldli.rowval.push_back(j);
            ldli.fval.push_back(f);
            ldli_row_ptr += 1;
            intops_count++;
        }

        LLcol llcol = colspace[len-1];
        Tval w = llcol.cval * colScale;
        flops_count += 1;
        // flop count: 1 mul
        Tind j = llcol.row;

        ldli.rowval.push_back(j);
        ldli.fval.push_back(1);
        ldli_row_ptr += 1;

        d[i] = w;
    }

    ldli.colptr[n-1] = ldli_row_ptr;
    ldli.d = d;
}

void forward(const LDLinv& ldli, std::vector<Tval>& y) {
    int len = ldli.col.size(); 
    for (int ii = 0; ii < len; ii++){ 
        Tind i = ldli.col[ii]; 

        Tind j0 = ldli.colptr[ii]; 
        Tind j1 = ldli.colptr[ii+1]-1;

        Tval yi = y[i];

        for (Tind jj = j0; jj <= j1-1; jj++){
            Tval j = ldli.rowval[jj];
            y[j] += ldli.fval[jj] * yi;
            yi *= ((Tval)(1)-ldli.fval[jj]);
        }
        Tval j = ldli.rowval[j1];
        y[j] += yi;
        y[i] = yi;
    }
} 

void backward(const LDLinv& ldli, std::vector<Tval>& y) {
    int len = ldli.col.size(); 
    for (int ii = len-1; ii >= 0; ii--) {
        Tind i = ldli.col[ii];

        Tind j0 = ldli.colptr[ii];
        Tind j1 = ldli.colptr[ii+1]-1;

        Tval j = ldli.rowval[j1];
        Tval yi = y[i];
        yi = yi + y[j];

        for (Tind jj = j1-1; jj >= j0; jj--) {
            j = ldli.rowval[jj];
            yi = ((Tind)1-ldli.fval[jj])*yi + ldli.fval[jj]*y[j];
        }
        y[i] = yi;
    }
} 

std::vector<Tval> LDLsolver(const LDLinv& ldli, const std::vector<Tval>& b){
    std::vector<Tval> y(b.begin(), b.end());

    forward(ldli, y);

    for (int i = 0; i < ldli.d.size(); i++){
        if (ldli.d[i] != 0)
            y[i] /= ldli.d[i];
    }

    backward(ldli, y);

    Tval mu = mean(y);
    for (auto& elem: y){
        elem = elem - mu;
    }

    return y;
}

void approxchol_lapGiven(SparseMatrix& A, SparseMatrix& lap_A,
                         const std::vector<Tval>& b,          
                         std::vector<Tval>& sol, SolverParameter paras) 
{
    LLMatOrd llmat = LLMatOrd(A);
    LDLinv ldli = approxChol(llmat);
    sol = LDLsolver(ldli, b);

    Tval error = norm(lap_A*sol-b)/norm(b); 

    if (paras.verbose){
        // std::cout << "Ratio of operator edges to original edges: " << 2 * ldli.fval.size() / nnz(a) << "\n";
        // std::cout << "ratio of max to min diagonal of laplacian : " << maximum(diag(la))/minimum(diag(la))) << "\n";
        std::cout << "Error after first iteration: " << error << "\n";
    }

    // iteratively solve the equation
    sol = pcg(lap_A, b, LDLsolver, sol, ldli, paras);
}


void axpy2(const Tval al, const std::vector<Tval>& p, std::vector<Tval>& x){
    int n = x.size();
    for (int i = 0; i < n; ++i)
        x[i] = x[i] + al*p[i];
}

void bzbeta(const Tval beta, const std::vector<Tval>& z, std::vector<Tval>& p){
    int n = p.size();
    for (int i = 0; i < n; ++i)
        p[i] = z[i] + beta*p[i];
}

Tval mean(const std::vector<Tval>& y){
    Tval sum = 0.0;
    for (auto& elem: y){
        sum += elem;
    }
    return sum/y.size();
}


Tval norm(const std::vector<Tval>& y){
    Tval sum = 0.0;
    for (auto& elem: y){
        sum += elem*elem;
    }
    return std::sqrt(sum);
}

Tval dot(const std::vector<Tval>& x, const std::vector<Tval>& y){
    Tval sum = 0.0;
    for (int i = 0; i < x.size(); ++i){
        sum += x[i]*y[i];
    }
    return sum;
}


std::vector<Tval> pcg(const SparseMatrix& la, const std::vector<Tval>& b,       
                      solver ldlsolver, const std::vector<Tval>& sol, 
                      const LDLinv& ldli, SolverParameter& paras)
{ 
    Tind n = b.size();
    Tval nb = norm(b);

    // If input vector is zero, quit
    if (nb == 0)
      return sol;

    std::vector<Tval> x(n, 0.0);
    std::vector<Tval> bestx(n, 0.0);
    Tval bestnr = 1.0;

    std::vector<Tval> r(b.begin(), b.end());
    std::vector<Tval> z(sol.begin(), sol.end());
    std::vector<Tval> p(sol.begin(), sol.end());
    std::vector<Tval> q(n, 0.0);

    Tval rho = dot(r, sol);
    Tval best_rho = rho;
    Tind stag_count = 0;
    time_t t1 = time(NULL);

    Tind itcnt = 0;
    Tval oldrho = rho; 
    while (itcnt < paras.maxits){
        itcnt = itcnt+1;

        q = la*p; 
        Tval pq = dot(p,q);
        // if (paras.verbose) 
        //     std::cout << "Current iteration: " << itcnt << ", pq = " 
        //     << pq << ", rho = " << rho << "\n";
        
        if (pq < EPS || pq > INF){ 
          if (paras.verbose)
            std::cout << "PCG Stopped due to small or large pq.\n";
          break;
        }

        Tval al = rho/pq;
        // if (paras.verbose) 
        //     std::cout << "al: " << al << ", nb = " 
        //     << nb << "\n";

        // the following line could cause slowdown
        if (al*norm(p) < EPS*norm(x)){
            if (paras.verbose) 
                std::cout << "PCG: Stopped due to stagnation.\n";
            break;
        }
        
        axpy2(al,p,x);
        axpy2(-al,q,r);

        Tval nr = norm(r)/nb;
        if (nr < bestnr){
            bestnr = nr;
            for (int i = 0; i < n; ++i)
                bestx[i] = x[i];
        }
                    
        if (nr < paras.tolerance) 
            break;
        z = ldlsolver(ldli, r);

        oldrho = rho;
        rho = dot(z, r); 

        if (paras.stag_test != 0) {
            if (rho < best_rho*(1-1/paras.stag_test)){
                best_rho = rho;
                stag_count = 0;
            } else {
                if (paras.stag_test > 0) {
                    if (best_rho > (1-1/paras.stag_test)*rho) {
                        stag_count += 1;
                        if (stag_count > paras.stag_test) {
                            std::cout << "PCG Stopped by stagnation test " << paras.stag_test << "\n";
                            break;
                        }
                    }
                }
            }
        }

        if (rho < EPS || rho > INF){
          if (paras.verbose)
            std::cout << "PCG Stopped due to small or large rho.\n"; 
          break;
        }


        Tval beta = rho/oldrho; 
        if (beta < EPS || beta > INF) {
          if (paras.verbose)
            std::cout << "PCG Stopped due to small or large beta.\n"; 
          break;
        }

        bzbeta(beta,z,p);

        if ((time(NULL) - t1) > paras.maxtime){
            if (paras.verbose) 
                std::cout << "PCG New stopped at maxtime.\n"; 
            break;
        }

    }

    if (paras.verbose)
        std::cout << "PCG stopped after: " << std::round((time(NULL) - t1)) << " seconds and " << itcnt << " iterations with relative error " << (norm(la*bestx-b)/norm(b)) << ".\n";

    if (paras.pcgIts.size() > 0)
        paras.pcgIts[0] = itcnt; 

    return bestx;
}

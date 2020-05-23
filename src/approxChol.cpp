#include "approxChol.hpp"
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <emmintrin.h>
#include <iostream>
#include <random>
#include "pcg_random.hpp"
#include <cmath>
#include <ctime>
#include <limits>

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

bool cmp_row_elem(const Elem &a, const Elem &b) {
    return a.row < b.row;
}

bool cmp_val_elem(const Elem &a, const Elem &b) {
    return a.val < b.val;
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

LDLinv approxChol_vector2_merge(LLMatOrd_vector2 a) {
    auto n = a.n;

    // need to make custom one without col info later
    LDLinv ldli(a);
    Tind ldli_row_ptr = 0;

    std::vector<Tval> d(n, 0);

    // std::vector<LLcol> colspace(n);
    std::vector<Tind> row(n);
    std::vector<Tval> val(n);


    std::vector<Tval> cumspace(n);

    // random engine and distribution
    std::default_random_engine engine;
    std::uniform_real_distribution<Tval> u(0.0, 1.0);


    for (long i = 0; i <= n-2; i++) {

        ldli.col[i] = i;
        ldli.colptr[i] = ldli_row_ptr;

        int len = a.row[i].size();
        // TODO: is sort by val necessary
        if (len > 0)
        {
            typedef std::sort_helper::value_iterator_t<Tval,Tind> IndexByRow;
            std::sort(IndexByRow(&a.val[i][0], &a.row[i][0]), IndexByRow(&a.val[i][0] + len, &a.row[i][0] + len) );
            row.clear();
            val.clear();

            int idx = 1;
            Tval last_val = a.val[i][0];
            Tind last_row = a.row[i][0];
            // std::cout << "merge" << std::endl;
            while (idx < len)
            {
                while (idx < len && a.row[i][idx] == last_row)
                    last_val += a.val[i][idx++];
                row.push_back(last_row);
                val.push_back(last_val);
                if (idx < len)
                {
                    last_val = a.val[i][idx];
                    last_row = a.row[i][idx];
                }
            }
            // std::cout << "merge done" << std::endl;
            len = row.size();
            typedef std::sort_helper::value_iterator_t<Tind,Tval> IndexByVal;
            std::sort(IndexByVal(&row[0], &val[0]), IndexByVal(&row[0] + len, &val[0] + len) );
        }
        // flop count: sort LLcol len*log(len)?

        Tval csum = 0;
        for (int ii = 0; ii < len; ii++) {
            csum += val[ii];
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

            ks[joffset] = row[koff];
        }

        Tval newEdgeVals[len-1];
        Tind js[len-1];
        Tind ptrs[len-1];
        for (int joffset = 0; joffset <= len-2; joffset++) {
            // LLcol llcol = colspace[joffset];
            Tval w = val[joffset] * colScale;
            js[joffset] = row[joffset];
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
                a.row[j].push_back(k);
                a.val[j].push_back(newEdgeVal);
            } else {        // put it in col k
                a.row[k].push_back(j);
                a.val[k].push_back(newEdgeVal);
            }
        }

        // LLcol llcol = colspace[len-1];
        Tval w = val[len - 1] * colScale;
        // flop count: 1 mul
        Tind j = row[len - 1];

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

LDLinv approxChol_vector2_struct_merge(LLMatOrd_vector2_struct a) {
    auto n = a.n;

    // need to make custom one without col info later
    LDLinv ldli(a);
    Tind ldli_row_ptr = 0;

    std::vector<Tval> d(n, 0);

    // std::vector<LLcol> colspace(n);
    std::vector<Elem> elems(n);


    std::vector<Tval> cumspace(n);

    // random engine and distribution
    std::default_random_engine engine;
    std::uniform_real_distribution<Tval> u(0.0, 1.0);


    for (long i = 0; i <= n-2; i++) {

        ldli.col[i] = i;
        ldli.colptr[i] = ldli_row_ptr;

        int len = a.elems[i].size();
        // TODO: is sort by val necessary
        if (len > 0)
        {
            std::sort(a.elems[i].begin(), a.elems[i].end(), cmp_row_elem);
            elems.clear();

            int idx = 1;
            Tval last_val = a.elems[i][0].val;
            Tind last_row = a.elems[i][0].row;
            // std::cout << "merge" << std::endl;
            while (idx < len)
            {
                while (idx < len && a.elems[i][idx].row == last_row)
                    last_val += a.elems[i][idx++].val;
                elems.push_back(Elem(last_row, last_val));
                if (idx < len)
                {
                    last_val = a.elems[i][idx].val;
                    last_row = a.elems[i][idx].row;
                }
            }
            // std::cout << "merge done" << std::endl;
            len = elems.size();
            std::sort(elems.begin(), elems.end(), cmp_val_elem);
        }
        // flop count: sort LLcol len*log(len)?

        Tval csum = 0;
        for (int ii = 0; ii < len; ii++) {
            csum += elems[ii].val;
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

            ks[joffset] = elems[koff].row;
        }

        Tval newEdgeVals[len-1];
        Tind js[len-1];
        Tind ptrs[len-1];
        for (int joffset = 0; joffset <= len-2; joffset++) {
            // LLcol llcol = colspace[joffset];
            Tval w = elems[joffset].val * colScale;
            js[joffset] = elems[joffset].row;
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
                a.elems[j].push_back(Elem(k, newEdgeVal));
            } else {        // put it in col k
                a.elems[k].push_back(Elem(j, newEdgeVal));
            }
        }

        // LLcol llcol = colspace[len-1];
        Tval w = elems[len - 1].val * colScale;
        // flop count: 1 mul
        Tind j = elems[len - 1].row;

        ldli.rowval.push_back(j);
        ldli.fval.push_back(1);
        ldli_row_ptr += 1;

        d[i] = w;

        //free column
        a.elems[i].clear();
        a.elems[i].shrink_to_fit();
    }

    ldli.colptr[n-1] = ldli_row_ptr;
    ldli.d = d;

    return ldli;
}

LDLinv approxChol_vector2_mergerand(LLMatOrd_vector2 a) {
    auto n = a.n;

    // need to make custom one without col info later
    LDLinv ldli(a);
    Tind ldli_row_ptr = 0;

    std::vector<Tval> d(n, 0);

    // std::vector<LLcol> colspace(n);
    std::vector<Tind> row(n);
    std::vector<Tval> val(n);


    std::vector<Tval> cumspace(n);

    // random engine and distribution
    pcg_extras::seed_seq_from<std::random_device> seed_source;
    pcg32_fast engine(seed_source);
    std::uniform_real_distribution<Tval> u(0.0, 1.0);
    std::uniform_real_distribution<Tval> uidx(0, n);
    Tval randnums[n];
    Tind randidx[n];
    for (int joffset = 0; joffset < n; joffset++)
    {
        randnums[joffset] = u(engine);
        randidx[joffset] = u(engine);
    }
    for (long i = 0; i <= n-2; i++) {

        ldli.col[i] = i;
        ldli.colptr[i] = ldli_row_ptr;

        int len = a.row[i].size();
        // TODO: is sort by val necessary
        if (len > 0)
        {
            typedef std::sort_helper::value_iterator_t<Tval,Tind> IndexByRow;
            std::sort(IndexByRow(&a.val[i][0], &a.row[i][0]), IndexByRow(&a.val[i][0] + len, &a.row[i][0] + len) );
            row.clear();
            val.clear();

            int idx = 1;
            Tval last_val = a.val[i][0];
            Tind last_row = a.row[i][0];
            // std::cout << "merge" << std::endl;
            while (idx < len)
            {
                while (idx < len && a.row[i][idx] == last_row)
                    last_val += a.val[i][idx++];
                row.push_back(last_row);
                val.push_back(last_val);
                if (idx < len)
                {
                    last_val = a.val[i][idx];
                    last_row = a.row[i][idx];
                }
            }
            // std::cout << "merge done" << std::endl;
            len = row.size();
            typedef std::sort_helper::value_iterator_t<Tind,Tval> IndexByVal;
            std::sort(IndexByVal(&row[0], &val[0]), IndexByVal(&row[0] + len, &val[0] + len) );
        }
        // flop count: sort LLcol len*log(len)?

        Tval csum = 0;
        for (int ii = 0; ii < len; ii++) {
            csum += val[ii];
            // flop count: 1 add
            cumspace[ii] = csum;
        }
        Tval wdeg = csum;

        Tval colScale = 1;

        Tind ks[len-1];

        for (int joffset = 0; joffset <= len-2; joffset++) {
            Tval r = randnums[randidx[(i+joffset)%n]];
            r = r * (csum - cumspace[joffset]) + cumspace[joffset];

            auto cumspace_last = cumspace.begin();
            std::advance(cumspace_last, len);
            int koff = std::distance(cumspace.begin(), std::lower_bound(cumspace.begin(), cumspace_last, r));

            ks[joffset] = row[koff];
        }

        Tval newEdgeVals[len-1];
        Tind js[len-1];
        Tind ptrs[len-1];
        for (int joffset = 0; joffset <= len-2; joffset++) {
            // LLcol llcol = colspace[joffset];
            Tval w = val[joffset] * colScale;
            js[joffset] = row[joffset];
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
                a.row[j].push_back(k);
                a.val[j].push_back(newEdgeVal);
            } else {        // put it in col k
                a.row[k].push_back(j);
                a.val[k].push_back(newEdgeVal);
            }
        }

        // LLcol llcol = colspace[len-1];
        Tval w = val[len - 1] * colScale;
        // flop count: 1 mul
        Tind j = row[len - 1];

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

    std::vector<Tind> row(n);
    std::vector<Tval> val(n);

    std::vector<Tval> cumspace(n);

    std::default_random_engine engine;
    std::uniform_real_distribution<Tval> u(0.0, 1.0);


    for (long i = 0; i <= n-2; i++) {

        ldli.col[i] = i;
        ldli.colptr[i] = ldli_row_ptr;

        int len = a.row[i].size();
        
        if (len > 0)
        {
            typedef std::sort_helper::value_iterator_t<Tval,Tind> IndexByRow;
            std::sort(IndexByRow(&a.val[i][0], &a.row[i][0]), IndexByRow(&a.val[i][0] + len, &a.row[i][0] + len) );
            row.clear();
            val.clear();

            int idx = 1;
            Tval last_val = a.val[i][0];
            Tind last_row = a.row[i][0];
            // std::cout << "merge" << std::endl;
            while (idx < len)
            {
                while (idx < len && a.row[i][idx] == last_row)
                    last_val += a.val[i][idx++];
                row.push_back(last_row);
                val.push_back(last_val);
                if (idx < len)
                {
                    last_val = a.val[i][idx];
                    last_row = a.row[i][idx];
                }
            }
            // std::cout << "merge done" << std::endl;
            len = row.size();
            typedef std::sort_helper::value_iterator_t<Tind,Tval> IndexByVal;
            std::sort(IndexByVal(&row[0], &val[0]), IndexByVal(&row[0] + len, &val[0] + len) );
        }
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
        
        Tval randnums[newlen];
        for (int joffset = 0; joffset <= len-2; joffset++)
        {
            randnums[joffset] = u(engine);
        }

        auto cumspace_last = cumspace.begin();
        std::advance(cumspace_last, len);

        for (int joffset = 0; joffset <= len-2; joffset++) {
            Tval r = randnums[joffset];
            r = r * (csum - cumspace[joffset]) + cumspace[joffset];

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

// with jkswap and sampling simd
LDLinv approxChol_vector2_opt2(LLMatOrd_vector2 a) {
    auto n = a.n;

    LDLinv ldli(a);
    Tind ldli_row_ptr = 0;

    std::vector<Tval> d(n, 0);

    std::vector<Tind> row(n);
    std::vector<Tval> val(n);

    std::vector<Tval> cumspace(n);

    std::default_random_engine engine;
    std::uniform_real_distribution<Tval> u(0.0, 1.0);


    for (long i = 0; i <= n-2; i++) {

        ldli.col[i] = i;
        ldli.colptr[i] = ldli_row_ptr;

        int len = a.row[i].size();
        if (len > 0)
        {
            typedef std::sort_helper::value_iterator_t<Tval,Tind> IndexByRow;
            std::sort(IndexByRow(&a.val[i][0], &a.row[i][0]), IndexByRow(&a.val[i][0] + len, &a.row[i][0] + len) );
            row.clear();
            val.clear();

            int idx = 1;
            Tval last_val = a.val[i][0];
            Tind last_row = a.row[i][0];
            // std::cout << "merge" << std::endl;
            while (idx < len)
            {
                while (idx < len && a.row[i][idx] == last_row)
                    last_val += a.val[i][idx++];
                row.push_back(last_row);
                val.push_back(last_val);
                if (idx < len)
                {
                    last_val = a.val[i][idx];
                    last_row = a.row[i][idx];
                }
            }
            // std::cout << "merge done" << std::endl;
            len = row.size();
            typedef std::sort_helper::value_iterator_t<Tind,Tval> IndexByVal;
            std::sort(IndexByVal(&row[0], &val[0]), IndexByVal(&row[0] + len, &val[0] + len) );
        }
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
        
        Tval randnums[newlen];
        for (int joffset = 0; joffset <= len-2; joffset++)
        {
            randnums[joffset] = u(engine);
        }

        __m256d csum_4 = _mm256_set1_pd(csum);
        for (int joffset = 0; joffset < newlen; joffset+=4) {
            __m256d r_4 = _mm256_loadu_pd(randnums+joffset);
            __m256d cum_4 = _mm256_loadu_pd(&cumspace[joffset]);

            __m256d res1_4 = _mm256_sub_pd(csum_4, cum_4);
            __m256d res2_4 = _mm256_fmadd_pd(r_4, res1_4, cum_4);

            _mm256_storeu_pd(randnums+joffset, res2_4);
        }

        auto cumspace_last = cumspace.begin();
        std::advance(cumspace_last, len);

        for (int joffset = 0; joffset <= len-2; joffset++) {
            Tval r = randnums[joffset];
            //r = r * (csum - cumspace[joffset]) + cumspace[joffset];

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

            __m128i mask2 = _mm_xor_si128(mask1, allone);
            __m128i res2 = _mm_blendv_epi8(j_4, k_4, mask2);

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


// inline
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
        for (joffset = 0; joffset <= len-2; joffset++) {
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

// inline + simd
LDLinv approxChol_opt2(LLMatOrd a) {
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

        int newlen = ceil(double(len-1)/4)*4;
        Tind js[newlen];
        Tind ks[newlen];
        
        Tval randnums[newlen];
        for (int joffset = 0; joffset <= len-2; joffset++)
        {
            randnums[joffset] = u(engine);
        }

        __m256d csum_4 = _mm256_set1_pd(csum);
        for (int joffset = 0; joffset < newlen; joffset+=4) {
            __m256d r_4 = _mm256_loadu_pd(randnums+joffset);
            __m256d cum_4 = _mm256_loadu_pd(&cumspace[joffset]);

            __m256d res1_4 = _mm256_sub_pd(csum_4, cum_4);
            __m256d res2_4 = _mm256_fmadd_pd(r_4, res1_4, cum_4);

            _mm256_storeu_pd(randnums+joffset, res2_4);
        }

        auto cumspace_last = cumspace.begin();
        std::advance(cumspace_last, len);

        for (int joffset = 0; joffset <= len-2; joffset++) {
            Tval r = randnums[joffset];
            //r = r * (csum - cumspace[joffset]) + cumspace[joffset];

            int koff = std::distance(cumspace.begin(), std::lower_bound(cumspace.begin(), cumspace_last, r));

            ks[joffset] = colspace[koff].row;
        }

        Tval newEdgeVals[len-1];
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

        __m128i allone = _mm_set1_epi32(0xffffffff);
        for (int joffset = 0; joffset < newlen; joffset+=4) {
            __m128i j_4 = _mm_loadu_si128((__m128i *)(js+joffset));
            __m128i k_4 = _mm_loadu_si128((__m128i *)(ks+joffset));

            __m128i mask1 = _mm_cmplt_epi32(k_4, j_4);
            __m128i res1 = _mm_blendv_epi8(j_4, k_4, mask1);

            __m128i mask2 = _mm_xor_si128(mask1, allone);
            __m128i res2 = _mm_blendv_epi8(j_4, k_4, mask2);

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
            
            Tind jhead = a.cols[j];
            a.lles[ptr].row = k;
            a.lles[ptr].next = jhead;
            a.lles[ptr].val = newEdgeVal;
            a.cols[j] = ptr;
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

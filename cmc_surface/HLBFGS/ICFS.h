#ifndef ICFS_H
#define ICFS_H

#include <cmath>
#include <numeric>
#include <algorithm>

int dicfs_(int *n, int *nnz, double *a,
           double *adiag, int *acol_ptr__, int *arow_ind__,
           double *l, double *ldiag, int *lcol_ptr__, int *
           lrow_ind__, int *p, double *alpha, int *iwa, double *
           wa1, double *wa2);

int dicf_(int *n, int *nnz, double *a,
          double *diag, int *col_ptr__, int *row_ind__, int *p,
          int *info, int *indr, int *indf, int *list,
          double *w);

int dsel2_(int *n, double *x, int *keys, int *k);

int insort_(int *n, int *keys);

int ihsort_(int *n, int *keys);

int dstrsol_(int *n, double *l, double *ldiag,
             int *jptr, int *indr, double *r__, char *task);

#endif

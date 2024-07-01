//if you have fast blas like MKL, ACML, ATLAS, GOTOBLAS
//define FAST_BLAS, and link them

#ifndef FAST_BLAS

#include "HLBFGS_BLAS.h"
#include <cmath>

double DDOT(const  int *n, const double *x, const int *incx,
            const double *y, const int *incy)
{
    int i__1;
    double ret_val;
    static int i__, m, ix, iy, mp1;
    static double dtemp;

    --y;
    --x;

    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0)
    {
        return ret_val;
    }
    if (*incx == 1 && *incy == 1)
    {
        goto L20;
    }


    ix = 1;
    iy = 1;
    if (*incx < 0)
    {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0)
    {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        dtemp += x[ix] * y[iy];
        ix += *incx;
        iy += *incy;
    }
    ret_val = dtemp;
    return ret_val;

L20:
    m = *n % 5;
    if (m == 0)
    {
        goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        dtemp += x[i__] * y[i__];
    }
    if (*n < 5)
    {
        goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5)
    {
        dtemp = dtemp + x[i__] * y[i__] + x[i__ + 1] * y[i__ + 1] + x[
            i__ + 2] * y[i__ + 2] + x[i__ + 3] * y[i__ + 3] + x[i__ +
                4] * y[i__ + 4];
    }
L60:
    ret_val = dtemp;
    return ret_val;
}

void DAXPY(const int *n, const double *alpha, const double *x, const int *incx,
           double *y, const int *incy)
{
    int i__1;

    static int i__, m, ix, iy, mp1;

    --y;
    --x;

    if (*n <= 0)
    {
        return;
    }
    if (*alpha == 0.)
    {
        return;
    }
    if (*incx == 1 && *incy == 1)
    {
        goto L20;
    }

    ix = 1;
    iy = 1;
    if (*incx < 0)
    {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0)
    {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        y[iy] += *alpha * x[ix];
        ix += *incx;
        iy += *incy;
    }
    return;

L20:
    m = *n % 4;
    if (m == 0)
    {
        goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        y[i__] += *alpha * x[i__];
    }
    if (*n < 4)
    {
        return;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4)
    {
        y[i__] += *alpha * x[i__];
        y[i__ + 1] += *alpha * x[i__ + 1];
        y[i__ + 2] += *alpha * x[i__ + 2];
        y[i__ + 3] += *alpha * x[i__ + 3];
    }
    return;
}

double DNRM2(const int *n, const double *x, const int *incx)
{
    int i__1, i__2;
    double ret_val, d__1;

    static int ix;
    static double ssq, norm, scale, absxi;

    --x;

    if (*n < 1 || *incx < 1)
    {
        norm = 0.;
    }
    else if (*n == 1)
    {
        norm = fabs(x[1]);
    }
    else
    {
        scale = 0.;
        ssq = 1.;

        i__1 = (*n - 1) * *incx + 1;
        i__2 = *incx;
        for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2)
        {
            if (x[ix] != 0.)
            {
                absxi = (d__1 = x[ix], fabs(d__1));
                if (scale < absxi)
                {
                    d__1 = scale / absxi;
                    ssq = ssq * (d__1 * d__1) + 1.;
                    scale = absxi;
                }
                else
                {
                    d__1 = absxi / scale;
                    ssq += d__1 * d__1;
                }
            }
        }
        norm = scale * std::sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;
}

void DSCAL(const int *n, const double *a, double *x, const int *incx)
{
    int i__1, i__2;

    static int i__, m, mp1, nincx;

    --x;

    if (*n <= 0 || *incx <= 0)
    {
        return;
    }
    if (*incx == 1)
    {
        goto L20;
    }


    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
    {
        x[i__] = *a * x[i__];
    }
    return;

L20:
    m = *n % 5;
    if (m == 0)
    {
        goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__)
    {
        x[i__] = *a * x[i__];
    }
    if (*n < 5)
    {
        return;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 5)
    {
        x[i__] = *a * x[i__];
        x[i__ + 1] = *a * x[i__ + 1];
        x[i__ + 2] = *a * x[i__ + 2];
        x[i__ + 3] = *a * x[i__ + 3];
        x[i__ + 4] = *a * x[i__ + 4];
    }
    return;
}

#endif


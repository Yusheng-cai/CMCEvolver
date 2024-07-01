#include "LineSearch.h"

#define SAFE_SEARCH

/* Subroutine */ int MCSRCH(int *n, double *x, double *f,
                            double *g, double *s, double *stp, double *ftol,
                            double *gtol, double *xtol, double *stpmin, double *
                            stpmax, int *maxfev, int *info, int *nfev, double *wa,
                            int *keep, double *rkeep, double *cg_dginit)
{
    /* System generated locals */
    double d__1;
    static int c__1 = 1;
    /* Local variables */
    static double dg, fm, fx, fy, dgm, dgx, dgy, fxm, fym, stx, sty;
    static double dgxm, dgym;

    static int infoc;
    static double finit;
    static double width;
    static double stmin, stmax;
    static bool stage1;
    static double width1, ftest1;
    static bool brackt;
    static double dginit, dgtest;

    /* Parameter adjustments */
    --rkeep;
    --keep;
    --wa;
    --s;
    --g;
    --x;

    /* Function Body */
    dg = rkeep[14];
    dginit = rkeep[15];
    dgm = rkeep[16];
    dgtest = rkeep[17];
    dgx = rkeep[18];
    dgxm = rkeep[19];
    dgy = rkeep[20];
    dgym = rkeep[21];
    finit = rkeep[22];
    fm = rkeep[23];
    ftest1 = rkeep[24];
    fx = rkeep[25];
    fxm = rkeep[26];
    fy = rkeep[27];
    fym = rkeep[28];
    stmax = rkeep[29];
    stmin = rkeep[30];
    stx = rkeep[31];
    sty = rkeep[32];
    width = rkeep[33];
    width1 = rkeep[34];
    infoc = keep[15];
    brackt = keep[16] != 0;
    stage1 = keep[17] != 0;
    if (*info == -1)
    {
        goto L45;
    }
    infoc = 1;
    if (*n <= 0 || *stp <= 0. || *ftol < 0. || *gtol < 0. || *xtol < 0. || *
        stpmin < 0. || *stpmax < *stpmin || *maxfev <= 0)
    {
        goto L1010;
    }
    dginit = cg_dginit? *cg_dginit : DDOT(n, &g[1], &c__1, &s[1], &c__1);
    if (dginit >= 0.)
    {
        goto L1010;
    }
    brackt = false;
    stage1 = true;
    *nfev = 0;
    finit = *f;
    dgtest = *ftol * dginit;
    width = *stpmax - *stpmin;
    width1 = width / .5;
    memcpy(&wa[1], &x[1], sizeof(double)*(*n));
    stx = 0.;
    fx = finit;
    dgx = dginit;
    sty = 0.;
    fy = finit;
    dgy = dginit;
L30:
    if (brackt)
    {
        stmin = std::min<double>(stx,sty);
        stmax = std::max<double>(stx,sty);
    }
    else
    {
        stmin = stx;
        stmax = *stp + (*stp - stx) * 4.;
    }
    *stp = std::max<double>(*stp,*stpmin);
    *stp = std::min<double>(*stp,*stpmax);
    if (brackt && (*stp <= stmin || *stp >= stmax) || *nfev >= *maxfev - 1 ||
        infoc == 0 || brackt && stmax - stmin <= *xtol * stmax)
    {
        *stp = stx;
    }
    memcpy(&x[1], &wa[1], sizeof(double)*(*n));
    DAXPY(n, stp, &s[1], &c__1, &x[1], &c__1);
    *info = -1;
    goto L1010;
L45:
    *info = 0;
    ++(*nfev);
    dg = DDOT(n, &g[1], &c__1, &s[1], &c__1);
    ftest1 = finit + *stp * dgtest;
    if (brackt && (*stp <= stmin || *stp >= stmax) || infoc == 0)
    {
        *info = 6;
    }
    if (*stp == *stpmax && *f <= ftest1 && dg <= dgtest)
    {
        *info = 5;
    }
    if (*stp == *stpmin && (*f > ftest1 || dg >= dgtest))
    {
        *info = 4;
    }
    if (*nfev >= *maxfev)
    {
        *info = 3;
    }
    if (brackt && stmax - stmin <= *xtol * stmax)
    {
        *info = 2;
    }
    //  if (*f <= ftest1 &&  (dg >= *gtol * dginit ||
    //      dg <= (2.0- (*gtol) ) * dginit) )  //TNPACK's C2 version, it seems not good.

    if (cg_dginit)
    {
        if (*f <= ftest1 &&  (dg >= *gtol * dginit ||
            dg <= (2.0- (*gtol) ) * dginit) )
        {
            *info = 1;
        }
    }
    else
    {
        if (*f <= ftest1 && std::fabs(dg) <= *gtol * (-dginit) )
        {
            *info = 1;
        }
    }

    if (*info != 0)
    {
        goto L1010;
    }
    if (stage1 && *f <= ftest1 && dg >= std::min<double>(*ftol,*gtol) * dginit)
    {
        stage1 = false;
    }
    if (stage1 && *f <= fx && *f > ftest1)
    {
        fm = *f - *stp * dgtest;
        fxm = fx - stx * dgtest;
        fym = fy - sty * dgtest;
        dgm = dg - dgtest;
        dgxm = dgx - dgtest;
        dgym = dgy - dgtest;
        MCSTEP(&stx, &fxm, &dgxm, &sty, &fym, &dgym, stp, &fm, &dgm, &brackt,
            &stmin, &stmax, &infoc);
        fx = fxm + stx * dgtest;
        fy = fym + sty * dgtest;
        dgx = dgxm + dgtest;
        dgy = dgym + dgtest;
    }
    else
    {
        MCSTEP(&stx, &fx, &dgx, &sty, &fy, &dgy, stp, f, &dg, &brackt, &
            stmin, &stmax, &infoc);
    }
    if (brackt)
    {
        if ((d__1 = sty - stx, std::fabs(d__1)) >= width1 * .66)
        {
            *stp = stx + (sty - stx) * .5;
        }
        width1 = width;
        width = (d__1 = sty - stx, std::fabs(d__1));
    }
    goto L30;
L1010:
    rkeep[14] = dg;
    rkeep[15] = dginit;
    rkeep[16] = dgm;
    rkeep[17] = dgtest;
    rkeep[18] = dgx;
    rkeep[19] = dgxm;
    rkeep[20] = dgy;
    rkeep[21] = dgym;
    rkeep[22] = finit;
    rkeep[23] = fm;
    rkeep[24] = ftest1;
    rkeep[25] = fx;
    rkeep[26] = fxm;
    rkeep[27] = fy;
    rkeep[28] = fym;
    rkeep[29] = stmax;
    rkeep[30] = stmin;
    rkeep[31] = stx;
    rkeep[32] = sty;
    rkeep[33] = width;
    rkeep[34] = width1;
    keep[15] = infoc;
    keep[16] = 0;
    if (brackt)
    {
        keep[16] = 1;
    }
    keep[17] = 0;
    if (stage1)
    {
        keep[17] = 1;
    }
    return 0;
} /* MCSRCH */

/* Subroutine */
int MCSTEP(double *stx, double *fx, double *dx,
           double *sty, double *fy, double *dy, double *stp,
           double *fp, double *dp, bool *brackt, double *stpmin,
           double *stpmax, int *info)
{
    /* System generated locals */
    double d__1, d__2, d__3;

    /* Local variables */
    static double p, q, r__, s, gama, sgnd, stpc, stpf, stpq, theta;
    static bool bound;

#ifdef SAFE_SEARCH
    /*     copy from TNPACK */
    static double xsafe;
    xsafe = .001;
#endif
    *info = 0;
    if (*brackt && (*stp <= std::min<double>(*stx,*sty) || *stp >= std::max<double>(*stx,*sty)) || *dx *
        (*stp - *stx) >= 0. || *stpmax < *stpmin)
    {
        return 0;
    }
    sgnd = *dp * (*dx / std::fabs(*dx));
    if (*fp > *fx)
    {
        *info = 1;
        bound = true;
        theta = (*fx - *fp) * 3 / (*stp - *stx) + *dx + *dp;
        /* Computing MAX */
        d__1 = std::fabs(theta), d__2 = std::fabs(*dx), d__1 = std::max<double>(d__1,d__2), d__2 = std::fabs(
            *dp);
        s = std::max<double>(d__1,d__2);
        /* Computing 2nd power */
        d__1 = theta / s;
        gama = s * std::sqrt(d__1 * d__1 - *dx / s * (*dp / s));
        if (*stp < *stx)
        {
            gama = -gama;
        }
        p = gama - *dx + theta;
        q = gama - *dx + gama + *dp;
        r__ = p / q;
        stpc = *stx + r__ * (*stp - *stx);
        stpq = *stx + *dx / ((*fx - *fp) / (*stp - *stx) + *dx) / 2 * (*stp -
            *stx);
        if ((d__1 = stpc - *stx, std::fabs(d__1)) < (d__2 = stpq - *stx, std::fabs(d__2)))
        {
            stpf = stpc;
        }
        else
        {
            stpf = stpc + (stpq - stpc) / 2;
        }
#ifdef SAFE_SEARCH
        /*       copy from TNPACK */
        if (*stp > *stx)
        {
            /* Computing MAX */
            d__1 = *stx + xsafe * (*stp - *stx);
            stpf = std::max<double>(d__1,stpf);
        }
        else
        {
            /* Computing MIN */
            d__1 = *stx + xsafe * (*stp - *stx);
            stpf = std::min<double>(d__1,stpf);
        }
#endif

        *brackt = true;
    }
    else if (sgnd < 0.)
    {
        *info = 2;
        bound = false;
        theta = (*fx - *fp) * 3 / (*stp - *stx) + *dx + *dp;
        /* Computing MAX */
        d__1 = std::fabs(theta), d__2 = std::fabs(*dx), d__1 = std::max<double>(d__1,d__2), d__2 = std::fabs(
            *dp);
        s = std::max<double>(d__1,d__2);
        /* Computing 2nd power */
        d__1 = theta / s;
        gama = s * std::sqrt(d__1 * d__1 - *dx / s * (*dp / s));
        if (*stp > *stx)
        {
            gama = -gama;
        }
        p = gama - *dp + theta;
        q = gama - *dp + gama + *dx;
        r__ = p / q;
        stpc = *stp + r__ * (*stx - *stp);
        stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
        if ((d__1 = stpc - *stp, std::fabs(d__1)) > (d__2 = stpq - *stp, std::fabs(d__2)))
        {
            stpf = stpc;
        }
        else
        {
            stpf = stpq;
        }
        *brackt = true;
    }
    else if (std::fabs(*dp) < std::fabs(*dx))
    {
        *info = 3;
        bound = true;
        theta = (*fx - *fp) * 3 / (*stp - *stx) + *dx + *dp;
        /* Computing MAX */
        d__1 = std::fabs(theta), d__2 = std::fabs(*dx), d__1 = std::max<double>(d__1,d__2), d__2 = std::fabs(
            *dp);
        s = std::max<double>(d__1,d__2);
        /* Computing MAX */
        /* Computing 2nd power */
        d__3 = theta / s;
        d__1 = 0., d__2 = d__3 * d__3 - *dx / s * (*dp / s);
        gama = s * std::sqrt((std::max<double>(d__1,d__2)));
        if (*stp > *stx)
        {
            gama = -gama;
        }
        p = gama - *dp + theta;
        q = gama + (*dx - *dp) + gama;
        r__ = p / q;
        if (r__ < 0. && gama != 0.)
        {
            stpc = *stp + r__ * (*stx - *stp);
        }
        else if (*stp > *stx)
        {
            stpc = *stpmax;
        }
        else
        {
            stpc = *stpmin;
        }
        stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
        if (*brackt)
        {
            if ((d__1 = *stp - stpc, std::fabs(d__1)) < (d__2 = *stp - stpq, std::fabs(
                d__2)))
            {
                stpf = stpc;
            }
            else
            {
                stpf = stpq;
            }
        }
        else
        {
            if ((d__1 = *stp - stpc, std::fabs(d__1)) > (d__2 = *stp - stpq, std::fabs(
                d__2)))
            {
                stpf = stpc;
            }
            else
            {
                stpf = stpq;
            }
        }
    }
    else
    {
        *info = 4;
        bound = false;
        if (*brackt)
        {
            theta = (*fp - *fy) * 3 / (*sty - *stp) + *dy + *dp;
            /* Computing MAX */
            d__1 = std::fabs(theta), d__2 = std::fabs(*dy), d__1 = std::max<double>(d__1,d__2), d__2 =
                std::fabs(*dp);
            s = std::max<double>(d__1,d__2);
            /* Computing 2nd power */
            d__1 = theta / s;
            gama = s * std::sqrt(d__1 * d__1 - *dy / s * (*dp / s));
            if (*stp > *sty)
            {
                gama = -gama;
            }
            p = gama - *dp + theta;
            q = gama - *dp + gama + *dy;
            r__ = p / q;
            stpc = *stp + r__ * (*sty - *stp);
            stpf = stpc;
        }
        else if (*stp > *stx)
        {
            stpf = *stpmax;
        }
        else
        {
            stpf = *stpmin;
        }
    }
#ifdef SAFE_SEARCH
    /*     copy from tnpack */
    sgnd = *dp * (*stx - *stp);
#endif
    if (*fp > *fx)
    {
        *sty = *stp;
        *fy = *fp;
        *dy = *dp;
    }
    else
    {
        if (sgnd < 0.)
        {
            *sty = *stx;
            *fy = *fx;
            *dy = *dx;
        }
        *stx = *stp;
        *fx = *fp;
        *dx = *dp;
    }
    stpf = std::min<double>(*stpmax,stpf);
    stpf = std::max<double>(*stpmin,stpf);
    *stp = stpf;
    if (*brackt && bound)
    {
        if (*sty > *stx)
        {
            /* Computing MIN */
            d__1 = *stx + (*sty - *stx) * .66;
            *stp = std::min<double>(d__1,*stp);
        }
        else
        {
            /* Computing MAX */
            d__1 = *stx + (*sty - *stx) * .66;
            *stp = std::max<double>(d__1,*stp);
        }
    }
    return 0;
} /* MCSTEP */


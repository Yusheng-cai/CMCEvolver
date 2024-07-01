#ifndef _HLBFGS_BLAS_H_
#define _HLBFGS_BLAS_H_

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

    double DDOT(const  int *n, const double *x, const int *incx,
        const double *y, const int *incy);

    void DAXPY(const int *n, const double *alpha, const double *x, const int *incx,
        double *y, const int *incy);

    double DNRM2(const int *n, const double *x, const int *incx);

    void  DSCAL(const int *n, const double *a, double *x, const int *incx);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _HLBFGS_BLAS_H_ */

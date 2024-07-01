#ifndef HLBFGS_H
#define HLBFGS_H

//////////////////////////////////////////////////////////////////////////

//! ICFS_INFO stores ICFS's working arrays
class ICFS_INFO
{
public:
    ICFS_INFO()
    {
        lcol_ptr = 0;
        lrow_ind = 0;
        ldiag = 0;
        l = 0;
        iwa = 0;
        wa1 = 0;
        wa2 = 0;
        p = 15;
        r = 0;
        CGP = 0;
        CGQ = 0;
        CGR = 0;
        CGZ = 0;
    }
    ~ICFS_INFO()
    {
        if (lcol_ptr)
        {
            delete[] lcol_ptr;
        }
        if (lrow_ind)
        {
            delete[] lrow_ind;
        }
        if (ldiag)
        {
            delete[] ldiag;
        }
        if (l)
        {
            delete[] l;
        }

        if (wa1)
        {
            delete[] wa1;
        }
        if (wa2)
        {
            delete[] wa2;
        }
        if (iwa)
        {
            delete[] iwa;
        }
        if (r)
        {
            delete[] r;
        }
        if (CGP)
            delete[] CGP;
        if (CGQ)
            delete[] CGQ;
        if (CGR)
            delete[] CGR;
        if (CGZ)
            delete[] CGZ;
    }
    void allocate_mem(int N)
    {
        if ( N > 0)
        {
            lcol_ptr = new int[N+1];
            ldiag = new double[N];
            iwa = new int[3*N];
            wa1 = new double[N];
            wa2 = new double[N];
            r = new double[N];
            p = 15;
            CGP = new double[N];
            CGR = new double[N];
            CGQ = new double[N];
            CGZ = new double[N];
        }
    }
public:
    int *lcol_ptr;
    int *lrow_ind;
    double *ldiag;
    double *l;
    int *iwa;
    double *wa1;
    double *wa2;
    int p;
    double *r;
    double *CGP;
    double *CGQ;
    double *CGR;
    double *CGZ;
    double icfs_alpha;

};
//////////////////////////////////////////////////////////////////////////
//! Stores the pointers of hessian matrix
class HESSIAN_MATRIX
{
public:
    HESSIAN_MATRIX(int N)
    {
        n = N;
        nnz = 0;
        values = 0;
        rowind = 0;
        colptr = 0;
        diag=0;
    }
    ~HESSIAN_MATRIX()
    {

    }
public:

    int n;
    int nnz;
    double *values;
    int *rowind;
    int *colptr;

    double *diag;
    ICFS_INFO l_info;

};
//////////////////////////////////////////////////////////////////////////
//! HLBFGS initialization
//Dimension of arrays: 20, 20
void INIT_HLBFGS(double PARAMETERS[], int INFO[]);

void HLBFGS_MESSAGE(bool print, int id, const double PARAMETERS[]);
//////////////////////////////////////////////////////////////////////////
void HLBFGS_UPDATE_First_Step(int N, int M, double *q, double *s, double *y,
                              double *rho, double *alpha, int bound, int cur_pos, int iter);

void HLBFGS_UPDATE_Hessian(int N, int M, double *q, double *s, double *y, int cur_pos, double *diag, int INFO[]);

void HLBFGS_UPDATE_Second_Step(int N, int M, double *q, double *s, double *y,
                               double *rho, double *alpha, int bound, int cur_pos, int iter);

void CONJUGATE_GRADIENT_UPDATE(int N, double *q, double *prev_q_update, double *prev_q_first_stage, int INFO[]);
//////////////////////////////////////////////////////////////////////////
void HLBFGS_BUILD_HESSIAN_INFO(HESSIAN_MATRIX& m_hessian, int INFO[]);
//////////////////////////////////////////////////////////////////////////
//! HLBFGS functions
void HLBFGS(int N, int M, double *x,
            void EVALFUNC(int,double*,double*,double*,double*),
            void EVALFUNC_H(int,double*,double*,double*,double*,HESSIAN_MATRIX&),
            void USER_DEFINED_HLBFGS_UPDATE_H(int,int,double*,double*,double*,int,double*, int[]),
            void NEWITERATION(int,int,double*,double*,double*,double*),
            double PARAMETERS[], int INFO[] );
//////////////////////////////////////////////////////////////////////////

#endif

//version 1.1
//Author: LIU Yang

#include <iostream>

#include "HLBFGS.h"

#include "LineSearch.h"
#include "ICFS.h"


//////////////////////////////////////////////////////////////////////////
void INIT_HLBFGS(double PARAMETERS[], int INFO[])
{
    PARAMETERS[0] = 1.0e-4; //ftol
    PARAMETERS[1] = 1.0e-16; //xtol
    PARAMETERS[2] = 0.9; //gtol
    PARAMETERS[3] = 1.0e-20; //stpmin
    PARAMETERS[4] = 1.0e+20; //stpmax
    PARAMETERS[5] = 1.0e-9; // ||g||/max(1,||x||)
    PARAMETERS[6] = 1.0e-10; // ||g||

    INFO[0] = 20; //max_fev_in_linesearch
    INFO[1] = 0; //total_num_fev
    INFO[2] = 0; //iter
    INFO[3] = 0; //update strategy. 0: standard lbfgs, 1: m1qn3; 
    INFO[4] = 100000; // max iterations
    INFO[5] = 1; //1: print message, 0: do nothing
    INFO[6] = 10; // T: update interval of Hessian
    INFO[7] = 0; // 0: without hessian, 1: with accurate hessian
    INFO[8] = 15; // icfs parameter
    INFO[9] = 0; // 0: linesearch 1: modified linesearch (do not set 1 in pratice !)
    INFO[10] = 0; // 0: Disable preconditioned CG 1: preconditioned CG
    INFO[11] = 1; // different methods for choosing beta in CG.
    INFO[12] = 1; //internal usage. 0: only update diag in USER_DEFINED_HLBFGS_UPDATE_H
    INFO[13] = 0; // 0: standard lbfgs update, 1: Biggs's update, 2: Yuan's update; 3: Zhang and Xu's update
}
//////////////////////////////////////////////////////////////////////////
void HLBFGS_MESSAGE(bool print, int id, const double PARAMETERS[])
{
    if (!print)
    {
        return;
    }
    switch (id)
    {
    case 0:
        std::cout << "Please check your input parameters !\n";
        break;
    case 1:
        std::cout <<"Linesearch is failed !\n";
        break;
    case 2:
        std::cout <<"Convergence : ||g||/max(1,||x||) <= " << PARAMETERS[5] << std::endl;
        break;
    case 3:
        std::cout <<"Convergence : ||g|| <=  " << PARAMETERS[6] << std::endl;
        break;
    case 4:
        std::cout <<"Convergence: linesearch cannot improve anymore \n";
        break;
    case 5:
        std::cout <<"Exceeds max iteration \n";
        break;
    default:
        break;
    }
}
//////////////////////////////////////////////////////////////////////////
void HLBFGS_UPDATE_Hessian(int N, int M, double *q, double *s, double *y,  int cur_pos, double *diag, int INFO[])
{
    if (M <= 0 || INFO[2] == 0)
    {
        return;
    }

    static int inc = 1;

    int start = cur_pos*N;

    double *y_start = &y[start];
    double *s_start = &s[start];

    double ys = DDOT(&N, y_start, &inc, s_start, &inc);

    if (INFO[3] == 0)
    {
        double yy = DDOT(&N, y_start, &inc, y_start, &inc);
        double factor = ys/yy;
        if (INFO[12] == 1)
        {
            DSCAL(&N, &factor, q, &inc);
        }
        else
        {
            diag[0] = factor;
        }

    }
    else if (INFO[3] == 1)
    {

        //m1qn3 update
        double dyy = 0;
        double dinvss = 0;
        for (int i = 0; i < N; i++)
        {
            dinvss += s_start[i]*s_start[i]/diag[i];
            dyy += diag[i]*y_start[i]*y_start[i];
        }
        for (int i = 0; i < N; i++)
        {
            diag[i] = 1.0/( dyy/(ys*diag[i]) + y_start[i]*y_start[i]/ys
                - dyy*s_start[i]*s_start[i]/(ys*dinvss*diag[i]*diag[i])   );
        }
        if (INFO[12] == 1)
        {
            for (int i = 0; i < N; i++)
            {
                q[i] *= diag[i];
            }
        }

    }
}
//////////////////////////////////////////////////////////////////////////
void HLBFGS_UPDATE_First_Step(int N, int M, double *q, double *s, double *y,
                              double *rho, double *alpha, int bound, int cur_pos, int iter)
{
    if (M <= 0)
    {
        return;
    }

    int start;
    double tmp;
    static int inc = 1;

    for (int i = bound; i >= 0; i--)
    {
        start = iter<=M? cur_pos-bound+i:(cur_pos-(bound-i)+M)%M;
        alpha[i] = rho[start] * DDOT(&N, q, &inc, &s[start*N], &inc);
        tmp = -alpha[i];
        DAXPY(&N, &tmp, &y[start*N], &inc, q, &inc);
    }
}

//////////////////////////////////////////////////////////////////////////
void HLBFGS_UPDATE_Second_Step(int N, int M, double *q, double *s, double *y,
                               double *rho, double *alpha, int bound, int cur_pos, int iter)
{
    if (M <= 0)
    {
        return;
    }

    int start;
    double tmp;
    static int inc = 1;

    for (int i = 0; i <= bound; i++)
    {
        start = iter<=M? i:(cur_pos+1+i)%M;
        tmp = alpha[i]-rho[start]*DDOT(&N, &y[start*N], &inc, q, &inc);
        DAXPY(&N, &tmp, &s[start*N], &inc, q, &inc);
    }
}
//////////////////////////////////////////////////////////////////////////
void HLBFGS_BUILD_HESSIAN_INFO(HESSIAN_MATRIX& m_hessian, int INFO[])
{
    m_hessian.l_info.p = INFO[8];

    if (m_hessian.l_info.lrow_ind)
    {
        delete[] m_hessian.l_info.lrow_ind;
    }
    m_hessian.l_info.lrow_ind = new int[m_hessian.nnz+m_hessian.n*m_hessian.l_info.p];
    if (m_hessian.l_info.l)
    {
        delete[] m_hessian.l_info.l;
    }
    m_hessian.l_info.l = new double[m_hessian.nnz+m_hessian.n*m_hessian.l_info.p];

    m_hessian.l_info.icfs_alpha = 0;

    dicfs_(&m_hessian.n, &m_hessian.nnz, m_hessian.values,
        m_hessian.diag, m_hessian.colptr, m_hessian.rowind,
        m_hessian.l_info.l, m_hessian.l_info.ldiag, m_hessian.l_info.lcol_ptr, m_hessian.l_info.lrow_ind,
        &m_hessian.l_info.p, &m_hessian.l_info.icfs_alpha, m_hessian.l_info.iwa, m_hessian.l_info.wa1, m_hessian.l_info.wa2);
}
//////////////////////////////////////////////////////////////////////////
void CONJUGATE_GRADIENT_UPDATE(int N, double *q, double *prev_q_update, double *prev_q_first_stage, int INFO[])
{
    static int inc = 1;
    //determine beta
    double cg_beta = 1.0;
    if (INFO[11]  == 1)
    {
        if (INFO[2] == 0)
        {
            memcpy(prev_q_first_stage, q, sizeof(double)*N);
            memcpy(prev_q_update, q, sizeof(double)*N);
            return;
        }
        else
        {
            cg_beta = DDOT(&N, q, &inc, q, &inc);
            cg_beta /= std::fabs( cg_beta - DDOT(&N, q, &inc, prev_q_first_stage, &inc));
            memcpy(prev_q_first_stage, q, sizeof(double)*N);
        }
    }
    else
    {
        if (INFO[2] == 0)
        {
            memcpy(prev_q_update, q, sizeof(double)*N);
            return;
        }
    }
    //determine new q

    const double minus_one = -1.0;
    if (cg_beta != 1.0)
        DSCAL(&N, &cg_beta, prev_q_update, &inc);

    DAXPY(&N, &minus_one, prev_q_update, &inc, q, &inc);
    double quad_a = DDOT(&N, q, &inc, q, &inc);
    double quad_b = DDOT(&N, q, &inc, prev_q_update, &inc);
    double cg_lambda = -quad_b / quad_a;
    if (cg_lambda > 1)
        cg_lambda = 1;
    else if (cg_lambda < 0)
        cg_lambda = 0;

    static double one = 1.0;
    DSCAL(&N, &cg_lambda, q, &inc);
    DAXPY(&N, &one, prev_q_update, &inc, q, &inc);

    memcpy(prev_q_update, q, sizeof(double)*N);
}
//////////////////////////////////////////////////////////////////////////
void HLBFGS(int N, int M, double *x,
            void EVALFUNC(int,double*,double*,double*,double*),
            void EVALFUNC_H(int,double*,double*,double*,double*,HESSIAN_MATRIX&),
            void USER_DEFINED_HLBFGS_UPDATE_H(int,int,double*,double*,double*,int,double*, int[]),
            void NEWITERATION(int,int,double*,double*,double*,double*),
            double PARAMETERS[], int INFO[] )
{
    int T = INFO[6];
    if ( N < 1 || M < 0 || T < -1 || INFO[4] < 1)
    {
        HLBFGS_MESSAGE(INFO[5]!=0, 0, PARAMETERS);
        return;
    }
    //allocate mem
    double *q = new double[N];
    double *g = new double[N];
    double *alpha = M<=0? 0: new double[M];
    double *rho = M<=0? 0: new double[M];
    double *s = M<=0? 0: new double[M*N];
    double *y = M<=0? 0: new double[M*N];
    double *prev_x = new double[N];
    double *prev_g = new double[N];
    double *diag = 0;
    double *wa = new double[N];
    double update_alpha = 1;
    HESSIAN_MATRIX m_hessian(N);
    if (INFO[3] == 1)
    {
        diag = new double[N];
        for (int i = 0; i < N; i++)
        {
            diag[i] = 1.0;
        }
    }
    double *prev_q_first_stage = 0;
    double *prev_q_update = 0;
    double scale = 0.0;
    double cg_dginit = 0;
    if (INFO[10] == 1)
    {
        if (INFO[11] == 1)
            prev_q_first_stage = new double[N];

        prev_q_update = new double[N];
    }

    //initialize
    static int inc = 1;
    INFO[1] = 0;
    INFO[2] = 0;
    double f = 0;
    int maxfev = INFO[0], bound = 0, nfev = 0, cur_pos = 0, start = 0;
    //line search parameters
    double stp, ftol=PARAMETERS[0], xtol=PARAMETERS[1], gtol=PARAMETERS[2],
        stpmin = PARAMETERS[3], stpmax = PARAMETERS[4];
    int info, keep[20];
    double gnorm, rkeep[40];
    memset(rkeep, 0, sizeof(double)*40);
    memset(keep, 0, sizeof(int)*20);

    m_hessian.l_info.allocate_mem(N);
    char task1='N';
    char task2='T';
    double prev_f;

    //////////////////////////////////////////////////////////////////////////
    do
    {
        if ( INFO[7] == 1 &&  ( (T==0) || ( INFO[2] % T == 0) ) )
        {
            //std::cout << "Generate Hessian\n";
            EVALFUNC_H(N, x, INFO[2]==0?0:prev_x, &f,  g, m_hessian);
            HLBFGS_BUILD_HESSIAN_INFO(m_hessian, INFO);
        }
        else if (INFO[2] == 0)
        {
            EVALFUNC(N, x, 0, &f, g);
            INFO[1]++;
        }

        if (INFO[2] > 0 && M > 0)
        {
            //compute s and y
            start = cur_pos*N;
            for (int i = 0; i < N; i++)
            {
                s[start+i] = x[i] - prev_x[i];
                y[start+i] = g[i] - prev_g[i];
            }
            rho[cur_pos] = 1.0/DDOT(&N, &y[start], &inc, &s[start], &inc);
            if (INFO[13] == 1)
            {
                update_alpha = 1.0 / (rho[cur_pos] * 6 * ( prev_f - f + DDOT(&N, g, &inc, &s[start], &inc) ) - 2.0);
            }
            else if (INFO[13] == 2)
            {
                update_alpha = 1.0 / (rho[cur_pos] * 2 * ( prev_f - f + DDOT(&N, g, &inc, &s[start], &inc) ) );
            }
            else if (INFO[13] == 3)
            {
                update_alpha = 1.0 / (1 + rho[cur_pos] * (6 *(prev_f - f) + 3*(DDOT(&N, g, &inc, &s[start], &inc)+DDOT(&N, prev_g, &inc, &s[start], &inc)) ) );
            }
            if (INFO[13] != 0)
            {
                if (update_alpha < 0.01)
                {
                    update_alpha = 0.01;
                }
                else if (update_alpha > 100)
                {
                    update_alpha = 100;
                }
                rho[cur_pos] *= update_alpha;
            }
        }

        for (int i = 0; i < N; i++)
        {
            q[i] = -g[i];
        }

        if (INFO[2] > 0 && M > 0)
        {
            bound = INFO[2] > M ? M-1:INFO[2]-1;
            HLBFGS_UPDATE_First_Step(N, M, q, s, y, rho, alpha, bound, cur_pos, INFO[2]);
        }

        if (INFO[10] == 0)
        {
            if (INFO[7] == 1)
            {
                dstrsol_(&N,m_hessian.l_info.l,m_hessian.l_info.ldiag,m_hessian.l_info.lcol_ptr,m_hessian.l_info.lrow_ind,q,&task1);
                dstrsol_(&N,m_hessian.l_info.l,m_hessian.l_info.ldiag,m_hessian.l_info.lcol_ptr,m_hessian.l_info.lrow_ind,q,&task2);
            }
            else
            {
                USER_DEFINED_HLBFGS_UPDATE_H(N, M, q, s, y, cur_pos, diag, INFO);
            }
        }
        else
        {
            if (INFO[7] == 1)
            {
                dstrsol_(&N,m_hessian.l_info.l,m_hessian.l_info.ldiag,m_hessian.l_info.lcol_ptr,m_hessian.l_info.lrow_ind,q,&task1);
                CONJUGATE_GRADIENT_UPDATE(N, q, prev_q_update, prev_q_first_stage, INFO);
                cg_dginit = -DDOT(&N, q, &inc, q, &inc);
                dstrsol_(&N,m_hessian.l_info.l,m_hessian.l_info.ldiag,m_hessian.l_info.lcol_ptr,m_hessian.l_info.lrow_ind,q,&task2);
            }
            else
            {
                INFO[12] = 0;
                USER_DEFINED_HLBFGS_UPDATE_H(N, M, q, s, y, cur_pos, INFO[3]==0 ? (&scale):diag, INFO);
                if (INFO[3] == 0)
                {
                    if (M > 0 && INFO[2] > 0 && scale != 1.0)
                    {
                        scale = std::sqrt(scale);
                        DSCAL(&N, &scale, q, &inc);
                    }
                    CONJUGATE_GRADIENT_UPDATE(N, q, prev_q_update, prev_q_first_stage, INFO);
                    cg_dginit = -DDOT(&N, q, &inc, q, &inc);
                    if (M > 0 && INFO[2] > 0 && scale != 1.0)
                        DSCAL(&N, &scale, q, &inc);
                }
                else
                {
                    if (M > 0 && INFO[2] > 0)
                    {
                        //use prev_g as temporary array
                        for (int i = 0; i < N; i++)
                        {
                            prev_g[i] = std::sqrt(diag[i]);
                            q[i] *= prev_g[i];
                        }
                    }
                    CONJUGATE_GRADIENT_UPDATE(N, q, prev_q_update, prev_q_first_stage, INFO);
                    cg_dginit = -DDOT(&N, q, &inc, q, &inc);
                    if (M > 0 && INFO[2] > 0)
                    {
                        for (int i = 0; i < N; i++)
                        {
                            q[i] *= prev_g[i];
                        }
                    }

                }
                INFO[12] = 1;
            }
        }



        if (INFO[2] > 0 &&  M > 0)
        {
            HLBFGS_UPDATE_Second_Step(N, M, q, s, y, rho, alpha, bound, cur_pos, INFO[2]);

            cur_pos = (cur_pos+1)%M;
        }

        //store g and x
        memcpy(prev_x, x, sizeof(double)*N);
        memcpy(prev_g, g, sizeof(double)*N);
        prev_f = f;
        //linesearch, find new x
        bool blinesearch = true;
        if (INFO[2] == 0)
        {
            gnorm = DNRM2(&N, g, &inc);
            //if(gnorm > 1)
            stp = 1.0/gnorm;
            //else
            //	stp = 1;
        }
        else
        {
            stp = 1;
        }

        info = 0;

        do
        {
            MCSRCH(&N, x, &f, g, q, &stp, &ftol, &gtol, &xtol, &stpmin, &stpmax, &maxfev, &info, &nfev, wa, keep, rkeep, INFO[10] == 0?0:(&cg_dginit));
            blinesearch =(info == -1);
            if (blinesearch)
            {
                EVALFUNC(N, x, prev_x, &f, g);
                INFO[1]++;
            }

            if (INFO[9] == 1 && prev_f > f) //modify line search to avoid too many function calls
            {
                info = 1;
                break;
            }

        }
        while (blinesearch);

        gnorm = DNRM2(&N, g, &inc);
        INFO[2]++;
        NEWITERATION(INFO[2], INFO[1], x, &f, g, &gnorm);
        double xnorm =DNRM2(&N, x, &inc);
        xnorm = 1>xnorm?1:xnorm;
        rkeep[2] = gnorm;
        rkeep[8] = xnorm;

        if (info != 1)
        {
            HLBFGS_MESSAGE(INFO[5]!=0, 1, PARAMETERS);
            break;
        }
        if (gnorm/xnorm <= PARAMETERS[5])
        {
            HLBFGS_MESSAGE(INFO[5]!=0, 2, PARAMETERS);
            break;
        }
        if (gnorm < PARAMETERS[6])
        {
            HLBFGS_MESSAGE(INFO[5]!=0, 3, PARAMETERS);
            break;
        }
        if (stp < stpmin || stp > stpmax)
        {
            HLBFGS_MESSAGE(INFO[5]!=0, 4, PARAMETERS);
            break;
        }
        if (INFO[2] > INFO[4])
        {
            HLBFGS_MESSAGE(INFO[5]!=0, 5, PARAMETERS);
            break;
        }

    }
    while (true);

    //free mem
    delete[] q;
    delete[] g;
    if (M > 0)
    {
        delete[] alpha;
        delete[] rho;
        delete[] s;
        delete[] y;
    }
    delete[] prev_x;
    delete[] prev_g;
    delete[] wa;
    if (diag)
        delete[] diag;
    if (prev_q_first_stage)
        delete[] prev_q_first_stage;
    if (prev_q_update)
        delete[] prev_q_update;
}

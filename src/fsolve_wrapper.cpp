#include "fsolve.hpp"
#include "fsolve_wrapper.hpp"

void fsolve_wrapper(std::function<void(Real& x, Real& fx)> func, Real ref, Real& guess, Real& func_val){
    double* x;
    double* fx;

    // set the fx and x 
    fx = new double[1];
    x  = new double[1];
    x[0] = guess;

    int n = 1;
    double tol = 0.00001;
    int lwa;
    double *wa;

    lwa = (n * (3.0 * n + 13.0)) / 2.0;
    wa  = new double[lwa];

    // capturing guess by value, func by reference
    auto f1 = [ref, &func](int n, double x[], double fx[]){
        Real xx, ffx;
        xx = (Real) x[0];
        func(xx, ffx);

        fx[0] = (double) ffx - (double) ref;

        return;
    };

    f1(n, x, fx);

    int info = fsolve(f1, n, x, fx, tol, wa, lwa);

    // set the final values 
    guess = x[0];
    func_val = fx[0];

    std::cout << "\n";
    std::cout << "  Returned value of INFO = " << info << "\n";

    delete[] fx;
    delete[] wa;
    delete[] x;
}

void fsolve_wrapper(std::function<void(std::vector<Real>& x, std::vector<Real>& fx)> func, std::vector<Real> ref, std::vector<Real>& guess, std::vector<Real>& func_val){
    int n = guess.size();

    double* x;
    double* fx;

    // set the fx and x 
    fx = new double[n];
    x  = new double[n];

    for (int i=0;i<n;i++){
        x[i] = guess[i];
    }

    double tol = 0.00001;
    int lwa;
    double *wa;

    lwa = (n * (3.0 * n + 13.0)) / 2.0;
    wa  = new double[lwa];

    auto f = [ref, &func](int n, double x[], double fx[]){
        std::vector<Real> xx(n,0);
        std::vector<Real> ffx(n,0);

        for (int i=0;i<n;i++){
            xx[i] = (Real) x[i];
        }

        func(xx, ffx);

        for (int i=0;i<n;i++){
            fx[i] = (double) ffx[i] - (double) ref[i];
        }

        return;
    };

    f(n, x, fx);

    int info = fsolve(f, n, x, fx, tol, wa, lwa);

    // set the final values 
    for (int i=0;i<n;i++){
        guess[i] = x[i];
        func_val[i] = fx[i];
    }

    std::cout << "\n";
    std::cout << "  Returned value of INFO = " << info << "\n";

    delete[] fx;
    delete[] wa;
    delete[] x;
}

std::pair<bool,Real> SimpleSolve(std::function<Real(Real)> &func,
                                   const Real xbegin,
                                   const Real xend)
    {
    Real a=xbegin;
    Real b=xend;
    const Real epsilon=1e-6;
    const Real delta=1e-8;
    Real f1=func(a);
    Real f2=func(b);

    if (std::fabs(f1)< delta){
        return std::make_pair(true, a);
    }

    if (std::fabs(f2) < delta){
        return std::make_pair(true,b);
    }
    std::cout << "start iterating." << std::endl;

    int it=0;
    while(true)
        {
        ++it;
        Real c=(b+a)/2;

        if((b-a)<epsilon*2.0)
            return std::make_pair(true,c);

        if(std::fabs(func(c))<delta)
            return std::make_pair(true,c);
        ((func(a)*func(c))<0) ? b=c : a=c ;
        if(it>1000)
            {
            return std::make_pair(false,c);
            }
        }
}
#include <functional>
#include "tools/CommonTypes.h"
#include "Eigen/unsupported/Eigen/NonLinearOptimization"
#include <cmath>

using Real = CommonTypes::Real;
using Real2= CommonTypes::Real2;

void fsolve_wrapper(std::function<void(Real& x, Real& fx)> func, Real ref, Real& x, Real& func_val);
void fsolve_wrapper(std::function<void(std::vector<Real>& x, std::vector<Real>& fx)> func, std::vector<Real> ref, std::vector<Real>& x, std::vector<Real>& func_val);
std::pair<bool,Real> SimpleSolve(std::function<Real(Real)> &func,
                                   const Real xbegin,
                                   const Real  xend);
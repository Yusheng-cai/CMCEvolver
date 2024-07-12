#include "AFP_shapes.h"
#include <cmath>

Real3 AFP_shape::Numericaldrdu(Real u, Real v){
    std::vector<Real> ulist = {u - 0.01f, u + 0.01f};
    std::vector<Real3> values;

    for (auto uu : ulist){
        values.push_back(calculatePos(uu, v));
    }

    Real3 drdu = 1.0 / 0.02 * (values[1] - values[0]);

    return drdu;
}

Real3 AFP_shape::Numericaldrdv(Real u, Real v){
    std::vector<Real> vlist = {v-0.01f, v+0.01f};
    std::vector<Real3> values;

    for (auto vv : vlist){
        values.push_back(calculatePos(u, vv));
    }

    Real3 drdv = 1.0 / 0.02 * (values[1] - values[0]);

    return drdv;
}

Real AFP_shape::CalculateUGivenTargetLength(Real init_u, Real v, Real target_len, Real init_u_step, Real tolerance, Real stepsize){
    std::vector<Real> length_vec;
    std::vector<Real> step_vec;
    Real3 origin,next_p,diff;
    Real dist;

    origin = this->calculatePos(init_u, v);
    length_vec.push_back(0);
    step_vec.push_back(0);

    // guess 2 
    next_p= this->calculatePos(init_u + init_u_step, v);
    diff  = next_p - origin;
    dist  = std::sqrt(LinAlg3x3::DotProduct(diff,diff));
    length_vec.push_back(dist);
    step_vec.push_back(init_u_step);

    // start iterating
    int ind = 0;
    Real ret_u;
    while (true){
        Real deriv = (step_vec[ind+1] - step_vec[ind]) / (length_vec[ind+1] - length_vec[ind]);
        Real next_u_step = step_vec[ind+1] + stepsize * deriv * (target_len - length_vec[ind+1]); 
        next_p     = this->calculatePos(init_u + next_u_step, v);

        diff       = next_p - origin;
        dist       = std::sqrt(LinAlg3x3::DotProduct(diff,diff));
        length_vec.push_back(dist);
        step_vec.push_back(next_u_step);

        if (std::abs(dist - target_len) < tolerance){
            ret_u = next_u_step + init_u;
            break;
        }

        ind++;
    }

    return ret_u;
}


Eigen::MatrixXd AFP_shape::NumericalJacobian(Real u, Real v){
    // jacobian matrix 
    Real3 drdu = Numericaldrdu(u,v);
    Real3 drdv = Numericaldrdv(u,v);

    Eigen::MatrixXd jac(3,2);
    jac << drdu[0], drdv[0], drdu[1], drdv[1], drdu[2], drdv[2];

    return jac;
}

Eigen::MatrixXd AFP_shape::InvNumericalJacobian(Real u, Real v){
    auto jac = NumericalJacobian(u,v);
    auto invjac = jac.completeOrthogonalDecomposition().pseudoInverse();

    return invjac;
}

Eigen::MatrixXd AFP_shape::InvJacobian(Real u, Real v, bool useNumerical){
    if (useNumerical){
        return InvNumericalJacobian(u,v);
    }
    else{
        return InvAnalyticalJacobian(u,v);
    }
}

Eigen::MatrixXd AFP_shape::AnalyticalJacobian(Real u, Real v){
    Real3 drdu = Analyticaldrdu(u,v);
    Real3 drdv = Analyticaldrdv(u,v);

    Eigen::MatrixXd jac(3,2);
    jac << drdu[0], drdv[0], drdu[1], drdv[1], drdu[2], drdv[2];

    return jac;
}

Eigen::MatrixXd AFP_shape::InvAnalyticalJacobian(Real u, Real v){
    auto jac = AnalyticalJacobian(u,v);
    auto invjac = jac.completeOrthogonalDecomposition().pseudoInverse();

    return invjac;
}

bool AFP_shape::CalculateNumericalNormalAndTangent(Real3& point, Real3& tangent, Real3& normal, int xdir, int ydir, int zdir){
    Real u = CalculateU(point);
    Real v = CalculateV(point);

    return CalculateNumericalNormalAndTangent(u,v,tangent, normal, xdir, ydir, zdir);
} 

bool AFP_shape::CalculateAnalyticalNormalAndTangent(Real3& point, Real3& tangent, Real3& normal, int xdir, int ydir, int zdir){
    Real u = CalculateU(point);
    Real v = CalculateV(point);

    return CalculateAnalyticalNormalAndTangent(u,v,tangent,normal,xdir,ydir,zdir);
}

Real3 AFP_shape::drdu(Real u, Real v, bool useNumerical){
    if (useNumerical){
        return Numericaldrdu(u,v);
    }
    else{
        return Analyticaldrdu(u,v);
    }
}

Real3 AFP_shape::drdv(Real u, Real v, bool useNumerical){
    if (useNumerical){
        return Numericaldrdv(u,v);
    }
    else{
        return Analyticaldrdv(u,v);
    }
}

bool AFP_shape::CalculateNumericalNormalAndTangent(Real u, Real v, Real3& tangent, Real3& normal, int xdir, int ydir, int zdir){
    Real3 drdu = Numericaldrdu(u,v);
    Real3 drdv = Numericaldrdv(u,v);

    Real3 norm = LinAlg3x3::CrossProduct(drdu, drdv);
    LinAlg3x3::normalize(norm);
    normal[xdir] = -norm[0];
    normal[ydir] = -norm[1];
    normal[zdir] = -norm[2];

    Real3 tang = drdv;
    LinAlg3x3::normalize(tang);
    tangent[xdir] = -tang[0];
    tangent[ydir] = -tang[1];
    tangent[zdir] = -tang[2];

    return true;
}

bool AFP_shape::CalculateAnalyticalNormalAndTangent(Real u, Real v, Real3& tangent, Real3& normal, int xdir, int ydir, int zdir){
    // // // initialize the position matrix 
    std::vector<std::vector<Real3>> pos(3, std::vector<Real3>(3,{0,0,0}));

    Real3 drdu = Analyticaldrdu(u,v);
    Real3 drdv = Analyticaldrdv(u,v);

    Real3 norm = LinAlg3x3::CrossProduct(drdu, drdv);
    LinAlg3x3::normalize(norm);
    normal[xdir] = norm[0];
    normal[ydir] = norm[1];
    normal[zdir] = norm[2];

    Real3 tang = drdv;
    LinAlg3x3::normalize(tang);
    tangent[xdir] = -tang[0];
    tangent[ydir] = -tang[1];
    tangent[zdir] = -tang[2];

    return true;
}

Real AFP_shape::shift_u_in_range(Real u){
    if (u < 0){
        u += 2 * Constants::PI;
    }

    return u;
}

                                            /// Sphere //// 

Sphere::Sphere(const ParameterPack& pack) : AFP_shape(pack){
    pack.ReadNumber("radius", ParameterPack::KeyType::Optional, radius_);
    pack.ReadArrayNumber("center", ParameterPack::KeyType::Required, center_);
}

Sphere::Real3 Sphere::calculatePos(Real u, Real v){
    Real3 ret;

    ret[0] = radius_ * std::sin(v) * std::cos(u) + center_[0];
    ret[1] = radius_ * std::sin(v) * std::sin(u) + center_[1];
    ret[2] = radius_ * std::cos(v);

    return ret;
}

Real Sphere::CalculateValue(Real3 position, int xdir, int ydir, int zdir){
    Real rsq = radius_ * radius_;
    Real val = std::pow(position[xdir] - center_[0],2) / rsq + std::pow(position[ydir] - center_[1],2) / rsq + \
    std::pow(position[zdir] - 0,2); 

    return val;
}

Real Sphere::CalculateAreaZ(Real z){
    Real phi = std::acos(z);
    Real r   = radius_ * std::sin(phi);

    return Constants::PI * r * r;
}

Real Sphere::CalculatePeriZ(Real z){
    Real phi = std::acos(z);
    Real r   = radius_ * std::sin(phi);

    return 2.0f * Constants::PI * r;
}

Real Sphere::CalculateV(Real3 pos, int xdir, int ydir, int zdir){
    return std::acos(pos[zdir] / radius_);
}

Real Sphere::CalculateU(Real3 pos, int xdir, int ydir, int zdir){
    Real dx = pos[xdir] - center_[0];
    Real dy = pos[ydir] - center_[1];

    // we strictly follow a 0->2pi U
    Real U  = std::atan2(dy, dx);
    U = shift_u_in_range(U);

    return U;
}

Real Sphere::CalculateU(Real3 pos, Real3 box, int xdir, int ydir, int zdir){
    Real dx = pos[xdir] - center_[0];
    Real dy = pos[ydir] - center_[1];

    if (dx > 0.5 * box[0]) {dx -= box[0];}
    else if (dx < -0.5 * box[0]) {dx += box[0];}

    if (dy > 0.5 * box[1]) {dy -= box[1];}
    else if (dy < -0.5 * box[1]) {dy += box[1];}

    Real U = std::atan2(dy, dx);
    U = shift_u_in_range(U);

    return U;
}

Real3 Sphere::Analyticaldrdu(Real u, Real v){
    Real3 ret;
    ret[0] = -radius_ * std::sin(v) * std::sin(u); 
    ret[1] = radius_  * std::sin(v) * std::cos(u);
    ret[2] = 0;

    return ret;
}

Real3 Sphere::Analyticaldrdv(Real u, Real v){
    Real3 ret;
    ret[0] = radius_ * std::cos(v) * std::cos(u);
    ret[1] = radius_ * std::cos(v) * std::sin(u);
    ret[2] = -radius_ * std::sin(v);

    return ret;
}


                                    ///// bulging sphere ////////
BulgingSphere::BulgingSphere(const ParameterPack& pack) : AFP_shape(pack){
    pack.ReadNumber("radius", ParameterPack::KeyType::Optional, radius_);
    pack.ReadArrayNumber("center", ParameterPack::KeyType::Required, center_);
    pack.ReadNumber("theta_phi_factor", ParameterPack::KeyType::Optional, theta_phi_factor_);
    pack.ReadNumber("bulge_factor", ParameterPack::KeyType::Optional, bulge_factor_);
}

BulgingSphere::Real BulgingSphere::rBulging(Real v){
    return radius_ * std::sin(v) + bulge_factor_ * std::sin(theta_phi_factor_ * v);
}

BulgingSphere::Real BulgingSphere::dradius_dv(Real u, Real v){
    return radius_ * std::cos(v) + bulge_factor_ * std::cos(theta_phi_factor_ * v) * theta_phi_factor_;
}

BulgingSphere::Real BulgingSphere::dradius_du(Real u, Real v){
    return 0.0f;
}

BulgingSphere::Real3 BulgingSphere::calculatePos(Real u, Real v){
    Real3 ret;

    Real r = rBulging(v);
    ret[0] = r * std::cos(u) + center_[0];
    ret[1] = r * std::sin(u) + center_[1];
    ret[2] = radius_ * std::cos(v);

    return ret;
}

BulgingSphere::Real BulgingSphere::CalculateAreaZ(Real z){
    Real v   = std::acos(z);
    Real r   = rBulging(v);

    return Constants::PI * r * r;
}

BulgingSphere::Real BulgingSphere::CalculatePeriZ(Real z){
    Real v   = std::acos(z);
    Real r   = rBulging(v);

    return 2.0f * Constants::PI * r;
}

BulgingSphere::Real BulgingSphere::CalculateValue(Real3 position, int xdir, int ydir, int zdir){
    Real v   = CalculateV(position, xdir, ydir, zdir);
    Real r   = rBulging(v);

    Real rsq = r * r;
    Real val = std::pow(position[xdir] - center_[0],2) / rsq + std::pow(position[ydir] - center_[1],2) / rsq;

    return val;
}

BulgingSphere::Real BulgingSphere::CalculateV(Real3 pos, int xdir, int ydir, int zdir){
    return std::acos(pos[zdir] / radius_);
}

BulgingSphere::Real BulgingSphere::CalculateU(Real3 pos, int xdir, int ydir, int zdir){
    Real dx = pos[xdir] - center_[0];
    Real dy = pos[ydir] - center_[1];

    // we strictly follow a 0->2pi U
    Real U  = std::atan2(dy, dx);
    U = shift_u_in_range(U);

    return U;
}

BulgingSphere::Real BulgingSphere::CalculateU(Real3 pos, Real3 box, int xdir, int ydir, int zdir){
    Real dx = pos[xdir] - center_[0];
    Real dy = pos[ydir] - center_[1];

    if (dx > 0.5 * box[0]) {dx -= box[0];}
    else if (dx < -0.5 * box[0]) {dx += box[0];}

    if (dy > 0.5 * box[1]) {dy -= box[1];}
    else if (dy < -0.5 * box[1]) {dy += box[1];}

    Real U = std::atan2(dy, dx);
    U = shift_u_in_range(U);

    return U;
}

BulgingSphere::Real3 BulgingSphere::Analyticaldrdu(Real u, Real v){
    Real3 ret;
    Real dra_du = BulgingSphere::dradius_du(u,v);
    Real r      = BulgingSphere::rBulging(v);

    ret[0] = dra_du * std::cos(u) - r * std::sin(u); 
    ret[1] = dra_du * std::cos(u) + r * std::cos(u);
    ret[2] = 0;

    return ret;
}

BulgingSphere::Real3 BulgingSphere::Analyticaldrdv(Real u, Real v){
    Real3 ret;

    Real dra_dv = BulgingSphere::dradius_dv(u,v);
    ret[0] = dra_dv * std::cos(u);
    ret[1] = dra_dv * std::sin(u);
    ret[2] = -radius_ * std::sin(v);

    return ret;
}



                                    ///// Super Egg ////////////

SuperEgg::SuperEgg(const ParameterPack& pack) : AFP_shape(pack){
    pack.ReadNumber("a", ParameterPack::KeyType::Required, a_);
    pack.ReadNumber("b", ParameterPack::KeyType::Required, b_);
    pack.ReadNumber("n", ParameterPack::KeyType::Optional, n_);
    pack.ReadNumber("zmax", ParameterPack::KeyType::Optional, zmax_);
    pack.ReadNumber("a_taper", ParameterPack::KeyType::Optional, a_taper_);
    pack.ReadNumber("b_taper", ParameterPack::KeyType::Optional, b_taper_);
    pack.ReadNumber("a_alpha", ParameterPack::KeyType::Optional, a_alpha_);
    pack.ReadNumber("b_alpha", ParameterPack::KeyType::Optional, b_alpha_);
    pack.ReadArrayNumber("center", ParameterPack::KeyType::Required, center_);
    pack.ReadNumber("internal_refinement", ParameterPack::KeyType::Optional, internal_refinement_);


    Real du = 2 * Constants::PI / internal_refinement_;
    for (int i=0;i<internal_refinement_;i++){
        internal_ulist_.push_back(du * i);
    }
}

Real SuperEgg::aBulging(Real z){
    return a_ * (1 - std::pow((a_taper_ * z / zmax_),2)) + a_alpha_ * std::sin(z / zmax_ * Constants::PI);
}

Real SuperEgg::da_dz(Real z){
    return a_ * a_taper_ * z / zmax_ * a_taper_ / zmax_ + a_alpha_ * std::cos(z / zmax_ * Constants::PI) * Constants::PI / zmax_;
}

Real SuperEgg::bBulging(Real z){
    return b_ * (1 - std::pow((b_taper_ * z / zmax_),2)) + b_alpha_ * std::sin(z / zmax_ * Constants::PI);
}

Real SuperEgg::db_dz(Real z){
    return b_ * (b_taper_ * z / zmax_);
}


SuperEgg::Real3 SuperEgg::calculatePos(Real u, Real v){
    Real3 ret;
    ret[2] = zmax_ * std::cos(v);
    ret[0] = aBulging(ret[2]) * Algorithm::sgn(std::cos(u)) * std::pow(std::abs(std::cos(u)), 2.0/n_) + center_[0];
    ret[1] = bBulging(ret[2]) * Algorithm::sgn(std::sin(u)) * std::pow(std::abs(std::sin(u)), 2.0/n_) + center_[1];

    return ret;
}

SuperEgg::Real SuperEgg::CalculateValue(Real3 position, int xdir, int ydir, int zdir){
    Real a,b;
    a = aBulging(position[zdir]);
    b = bBulging(position[zdir]);

    return std::pow(std::abs(position[xdir] - center_[0]) / a, n_) + std::pow(std::abs(position[ydir] - center_[1]) / b, n_);
}


Real SuperEgg::CalculateV(Real3 point, int xdir, int ydir, int zdir){
    Real v = std::acos((point[zdir]) / zmax_);

    return v;
}

Real SuperEgg::CalculateU(Real3 point, Real3 box, int xdir, int ydir, int zdir){
    // first get a good guess of what u and v are
    Real dx,dy,dz;
    Real3 shifted_pos;
    Real shifted_x, shifted_y;

    dx = point[xdir] - center_[0];
    dy = point[ydir] - center_[1];

    if (dx > box[xdir] * 0.5) {dx = dx - box[xdir];}
    else if (dx < -box[xdir] * 0.5) {dx = dx + box[xdir];}

    if (dy > box[ydir] * 0.5) {dy = dy - box[ydir];}
    else if (dy < -box[ydir] * 0.5) {dy = dy + box[ydir];}

    shifted_x = center_[0] + dx;
    shifted_y = center_[1] + dy;
    shifted_pos[0] = shifted_x; shifted_pos[1] = shifted_y;
    shifted_pos[2] = point[zdir];

    Real u = CalculateU(shifted_pos, xdir, ydir, zdir);
    u = shift_u_in_range(u);

    return u;
}


Real SuperEgg::CalculateU(Real3 point, int xdir, int ydir, int zdir){
    // first get a good guess of what u and v are
    Real u = std::atan2(point[ydir] - center_[1], point[xdir] - center_[0]);
    Real v = CalculateV(point, xdir, ydir, zdir);

    // have a list of u solutions
    std::vector<Real> u_solutions;
    std::vector<Real> err_list;
    std::vector<Real3> pos_s;

    // solve y first 
    FunctorY funcy(point, this, xdir, ydir, zdir);
    Eigen::HybridNonLinearSolver<FunctorY> NL_solver(funcy);
    Eigen::VectorXd y1,y2;
    Real3 y1_pos, y2_pos;
    Real err_y1, err_y2;
    y1.resize(1), y2.resize(1);
    y1[0] = Algorithm::sgn(u) * Constants::PI  - Algorithm::sgn(u) * 0.1;
    y2[0] = 0 + Algorithm::sgn(u) * 0.1;
    auto infoy = NL_solver.solveNumericalDiff(y1);
    auto infoy2= NL_solver.solveNumericalDiff(y2);

    // keep track of solutions
    u_solutions.push_back(y1[0]);
    u_solutions.push_back(y2[0]);
    y1_pos = calculatePos(y1[0], v), y2_pos = calculatePos(y2[0], v);
    pos_s.push_back(y1_pos); pos_s.push_back(y2_pos);
    err_y1 = std::pow(y1_pos[0] - point[xdir],2) + std::pow(y1_pos[1] - point[ydir],2);
    err_y2 = std::pow(y2_pos[0] - point[xdir],2) + std::pow(y2_pos[1] - point[ydir],2);
    err_list.push_back(err_y1);
    err_list.push_back(err_y2);

    // solve x 
    FunctorX funcx(point,this, xdir, ydir, zdir);
    Eigen::HybridNonLinearSolver<FunctorX> NL_solverX(funcx);
    Eigen::VectorXd x1,x2;
    Real3 x1_pos, x2_pos;
    Real err_x1, err_x2;
    x1.resize(1), x2.resize(1);
    x1[0] = Algorithm::sgn(u) * Constants::PI  - Algorithm::sgn(u) * 0.1;
    x2[0] = 0 + Algorithm::sgn(u) * 0.1;
    auto infox = NL_solverX.solveNumericalDiff(x1);
    auto infox2= NL_solverX.solveNumericalDiff(x2);

    // keep track of solutions
    u_solutions.push_back(x1[0]);
    u_solutions.push_back(x2[0]);
    x1_pos = calculatePos(x1[0],v), x2_pos = calculatePos(x2[0],v);
    pos_s.push_back(x1_pos); pos_s.push_back(x2_pos);
    err_x1 = std::pow(x1_pos[0] - point[xdir],2) + std::pow(x1_pos[1] - point[ydir],2);
    err_x2 = std::pow(x2_pos[0] - point[xdir],2) + std::pow(x2_pos[1] - point[ydir],2);
    err_list.push_back(err_x1);
    err_list.push_back(err_x2);

    // set u 
    int index = Algorithm::argmin(err_list);

    if (u_solutions[index] < -2*Constants::PI || u_solutions[index] > 2*Constants::PI){
        return std::fmod(u_solutions[index], 2*Constants::PI);
    }
    else{
        return u_solutions[index];
    }
}

Real SuperEgg::CalculateAreaZ(Real z){
    // green's theorem says A = 0.5 * \int x dy - y dx 
    Real A = 0.0f;

    Real3 p= {0,0,z};
    Real v = CalculateV(p, 0,1,2);

    std::vector<Real3> pos_list;
    for (int i=0;i<internal_refinement_;i++){
        pos_list.push_back(calculatePos(internal_ulist_[i], v));
    }

    for (int i=0;i<internal_refinement_;i++){
        Real3 diff = pos_list[(i+1) % internal_refinement_] - pos_list[i];
        A += pos_list[i][0] * diff[1] - pos_list[i][1] * diff[0];
    }

    return 0.5 * A;
}

Real SuperEgg::CalculatePeriZ(Real z){
    // we calculate perimeter by summing up the little chuncks
    Real3 p= {0,0,z};
    Real v = CalculateV(p, 0,1,2);

    std::vector<Real3> pos_list;
    for (int i=0;i<internal_refinement_;i++){
        pos_list.push_back(calculatePos(internal_ulist_[i], v));
    }

    Real P = 0.0f;

    for (int i=0;i<internal_refinement_;i++){
        Real3 diff = pos_list[(i+1) % internal_refinement_] - pos_list[i];
        P += LinAlg3x3::norm(diff);
    }

    return P;
}

Real3 SuperEgg::Analyticaldrdu(Real u, Real v){
    Real3 drdu;
    Real dydu, dxdu;
    Real z = zmax_ * std::cos(v);

    if (u >= 0){
        dydu = bBulging(z) * 2.0 / (float)n_ * std::pow(std::sin(u), 2.0/(float)n_ - 1) * std::cos(u);
        dxdu = aBulging(z) * 2.0 / (float)n_ * std::pow(std::cos(u), 2.0/(float)n_ - 1) * (-std::sin(u));
    }
    else{
        dydu = bBulging(z) * 2.0 / (float)n_ * std::pow(-std::sin(u), 2.0/(float)n_-1) * (-std::cos(u));
        dxdu = aBulging(z) * 2.0 / (float)n_ * std::pow(-std::cos(u), 2.0/(float)n_-1) * (std::sin(u));
    }

    drdu[0] = dxdu;
    drdu[1] = dydu; 
    drdu[2] = 0.0;

    return drdu;
}

Real3 SuperEgg::Analyticaldrdv(Real u, Real v){
    return {0,0,zmax_ * std::cos(v)};
}
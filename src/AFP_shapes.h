#pragma once

#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "tools/Constants.h"
#include "LinAlgTools.h"
#include "tools/Algorithm.h"
#include "fsolve_wrapper.hpp"
#include "Eigen/Core"


#include <vector>

class AFP_shape {
    public:
        using Real3 = CommonTypes::Real3;
        using Real  = CommonTypes::Real;
        typedef Eigen::Matrix<Real, 3, 2> Matrix3x2;
        typedef Eigen::Matrix<Real, 2, 3> Matrix2x3;

        AFP_shape(const ParameterPack& pack){};

        // calculate numerical normal and tangent given a point
        virtual bool CalculateNumericalNormalAndTangent(Real3& point, Real3& tangent, Real3& normal, int xdir=0, int ydir=1, int zdir=2); 
        // calculate analytical normal and tangent given a point
        virtual bool CalculateAnalyticalNormalAndTangent(Real3& point, Real3& tangent, Real3& normal, int xdir=0, int ydir=1, int zdir=2);
        // calculate numerical normal and tangent given u,v
        virtual bool CalculateNumericalNormalAndTangent(Real u, Real v, Real3& tangent, Real3& normal, int xdir=0, int ydir=1, int zdir=2);
        // calculate analytical normal and tangent given u,v
        virtual bool CalculateAnalyticalNormalAndTangent(Real u, Real v, Real3& tangent, Real3& normal, int xdir=0, int ydir=1, int zdir=2);
        // calculate numerical drdu/drdv --> depending on the calculate Pos function
        virtual Real3 Numericaldrdu(Real u, Real v);
        virtual Real3 Numericaldrdv(Real u, Real v);
        virtual Real3 drdu(Real u, Real v, bool useNumerical=true);
        virtual Real3 drdv(Real u, Real v, bool useNumerical=true);
        // calculate the numerical jacobian
        virtual Eigen::MatrixXd NumericalJacobian(Real u, Real v);
        virtual Eigen::MatrixXd InvNumericalJacobian(Real u, Real v);
        virtual Eigen::MatrixXd AnalyticalJacobian(Real u, Real v);
        virtual Eigen::MatrixXd InvAnalyticalJacobian(Real u, Real v);
        virtual Eigen::MatrixXd InvJacobian(Real u, Real v, bool useNumerical=true);

                                // pure virtual functions //
        // calculate the function value given position
        virtual Real CalculateValue(Real3 position, int xdir=0, int ydir=1, int zdir=2)=0;
        // calculate V
        virtual Real CalculateV(Real3 pos, int xdir=0, int ydir=1, int zdir=2)=0;
        // calculate U
        virtual Real CalculateU(Real3 pos, int xdir=0, int ydir=1, int zdir=2)=0;
        // calculate U for pbc
        virtual Real CalculateU(Real3 pos, Real3 box, int xdir=0, int ydir=1, int zdir=2)=0;
        virtual Real3 calculatePos(Real u, Real v)=0;
        virtual Real3 Analyticaldrdu(Real u, Real v)=0;
        virtual Real3 Analyticaldrdv(Real u, Real v)=0;

        virtual Real CalculateAreaZ(Real z)=0;
        virtual Real CalculatePeriZ(Real z)=0;

        Real shift_u_in_range(Real u);

        Real CalculateUGivenTargetLength(Real init_u, Real v, Real target_len, Real init_u_step=0.001, Real tolerance=1e-6, Real stepsize=0.1);

    protected:
        Real2 center_;
};

class Sphere : public AFP_shape{
    public:
        Sphere(const ParameterPack& pack);

        // calculate position given u,v
        virtual Real CalculateValue(Real3 position, int xdir=0, int ydir=1, int zdir=2) override;
        virtual Real CalculateV(Real3 pos, int xdir=0, int ydir=1, int zdir=2) override;
        virtual Real CalculateU(Real3 pos, int xdir=0, int ydir=1, int zdir=2) override;
        virtual Real CalculateU(Real3 pos, Real3 box, int xdir=0, int ydir=1, int zdir=2) override;
        virtual Real3 calculatePos(Real u, Real v) override;
        virtual Real3 Analyticaldrdu(Real u, Real v) override;
        virtual Real3 Analyticaldrdv(Real u, Real v) override;
        virtual Real CalculateAreaZ(Real z) override;
        virtual Real CalculatePeriZ(Real z) override;


    private:
        Real radius_;
};

class BulgingSphere : public AFP_shape{
    public:
        BulgingSphere(const ParameterPack& pack);

        // calculate position given u,v
        Real rBulging(Real v);
        Real dradius_dv(Real u, Real v);
        Real dradius_du(Real u, Real v);
        virtual Real CalculateValue(Real3 position, int xdir=0, int ydir=1, int zdir=2) override;
        virtual Real CalculateV(Real3 pos, int xdir=0, int ydir=1, int zdir=2) override;
        virtual Real CalculateU(Real3 pos, int xdir=0, int ydir=1, int zdir=2) override;
        virtual Real CalculateU(Real3 pos, Real3 box, int xdir=0, int ydir=1, int zdir=2) override;
        virtual Real3 calculatePos(Real u, Real v) override;
        virtual Real3 Analyticaldrdu(Real u, Real v) override;
        virtual Real3 Analyticaldrdv(Real u, Real v) override;
        virtual Real CalculateAreaZ(Real z) override;
        virtual Real CalculatePeriZ(Real z) override;

    private:
        Real radius_=1.0;
        Real bulge_factor_;
        Real theta_phi_factor_=2.0;
};

class SuperEgg : public AFP_shape{
    public:
        using Real3 = CommonTypes::Real3;
        using Real  = CommonTypes::Real;
        using Real2 = CommonTypes::Real2;

        SuperEgg(const ParameterPack& pack);

        Real aBulging(Real z);
        Real bBulging(Real z);
        Real da_dz(Real z);
        Real db_dz(Real z);

        virtual Real CalculateValue(Real3 position, int xdir=0, int ydir=1, int zdir=2) override;
        virtual Real CalculateV(Real3 point, int xdir=0, int ydir=1, int zdir=2) override;
        virtual Real CalculateU(Real3 pos, int xdir=0, int ydir=1, int zdir=2) override;
        virtual Real CalculateU(Real3 pos, Real3 box, int xdir=0, int ydir=1, int zdir=2);
        virtual Real3 calculatePos(Real u, Real v) override;
        virtual Real3 Analyticaldrdu(Real u, Real v) override;
        virtual Real3 Analyticaldrdv(Real u, Real v) override;
        virtual Real CalculateAreaZ(Real z) override;
        virtual Real CalculatePeriZ(Real z) override;

        Real getn() {return (Real)n_;}
        Real2 getcenter() {return center_;}

        struct FunctorY{
            Real3 point_;
            SuperEgg* se_;
            int xdir_, ydir_, zdir_;
            FunctorY(Real3 point, SuperEgg* se, int xdir, int ydir, int zdir) : point_(point), se_(se), xdir_(xdir), ydir_(ydir), zdir_(zdir) {};

            int operator()(Eigen::VectorXd &b, Eigen::VectorXd &fvec){
                fvec[0] = se_->bBulging(point_[zdir_] + se_->offset_height) * Algorithm::sgn(std::sin(b[0])) * std::pow(std::abs(std::sin(b[0])), 2.0/se_->getn()) \
                - (point_[ydir_] - se_->getcenter()[1]);

                return 0;
            }
            
            int operator()(Eigen::VectorXd &b, Eigen::VectorXd &fvec) const{
                fvec[0] = se_->bBulging(point_[zdir_] + se_->offset_height) * Algorithm::sgn(std::sin(b[0])) * std::pow(std::abs(std::sin(b[0])), 2.0/se_->getn()) \
                - (point_[ydir_] - se_->getcenter()[1]);

                return 0;
            }
        };

        struct FunctorX{
            Real3 point_;
            SuperEgg* se_;
            int xdir_, ydir_, zdir_;
            FunctorX(Real3 point, SuperEgg* se, int xdir, int ydir, int zdir) : point_(point), se_(se), xdir_(xdir), ydir_(ydir), zdir_(zdir) {};

            int operator()(Eigen::VectorXd &b, Eigen::VectorXd &fvec){
                fvec[0] = se_->aBulging(point_[zdir_] + se_->offset_height) * Algorithm::sgn(std::cos(b[0])) * std::pow(std::abs(std::cos(b[0])), 2.0/se_->getn()) \
                - (point_[xdir_] - se_->getcenter()[0]);

                return 0;
            }
            
            int operator()(Eigen::VectorXd &b, Eigen::VectorXd &fvec) const{
                fvec[0] = se_->aBulging(point_[zdir_] + se_->offset_height) * Algorithm::sgn(std::cos(b[0])) * std::pow(std::abs(std::cos(b[0])), 2.0/se_->getn()) \
                - (point_[xdir_] - se_->getcenter()[0]);

                return 0;
            }
        };

    private:
        Real a_;
        Real b_;
        Real n_=4;
        Real n_taper_=2.0;
        Real zmax_=1.5;
        Real a_taper_=1.0;
        Real b_taper_=0.5;
        Real a_alpha_=0.0;
        Real b_alpha_=0.0;
        Real offset_height=0.0;
        Real2 center_;
        int internal_refinement_=100000;

        std::vector<Real> internal_ulist_;
};

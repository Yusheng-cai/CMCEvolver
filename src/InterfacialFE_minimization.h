#pragma once 

#include "MeshRefineStrategy.h"
#include "Curvature.h"
#include "tools/CommonTypes.h"
#include "MeshCurvatureflow.h"
#include "AFP_shapes.h"
#include "Eigen/Core"
#include "tools/CommonOperations.h"

#include <vector>
#include <string>
#include <array>
#include <memory>
#include <chrono>
#include <algorithm>
#include <numeric>

// We are trying to solve the general equation of 
// dxi/dt = (k0 - kmean) * ni
// However, our kmean could be what we calculate 

// class InterfacialFE_Solver{
//     using INT3 = CommonTypes::index3;

//     private:
//         std::vector<INT3> face_;
//     public:
//         InterfacialFE_Solver(std::vector<INT3>& faces);
//         double operator()(Const Eigen::VectorXd& x, Eigen::VectorXd& grad);
// }

class InterfacialFE_minimization : public MeshRefineStrategy{
    public:
        InterfacialFE_minimization(MeshRefineStrategyInput& input);

        void update_Mesh();

        virtual void refine(Mesh& mesh) override;
        void refineBoundary(Mesh& mesh, AFP_shape* shape);

        const std::vector<Real>& getFE() {return FE_;}

        // function to obtain the gamma
        Real CalculateGamma(Real temperature);
        Real CalculateMu(Real temperature);

        Real getL() {return L_;}
        Real getK() {return rho_ * (L_ + mu_) / (2*gamma_);}
        Real getarea() {return area_list_[area_list_.size()-1];}
        Real getvolume() {return volume_list_[volume_list_.size()-1];}
        Real getVnbs() {return Vnbs_list_[Vnbs_list_.size()-1];}
        Real getAnbs() {return Anbs_list_[Anbs_list_.size()-1];}
        Real getVunderneath() {return Vunderneath_;}
        Real getVnbs_underneath() {return Vnbs_underneath_;}
        Mesh getMeshFlatContact() {return m_flatContact_;}

        void setL(Real L) {L_=L;}
        void setL2(Real L2) {L2_=L2;}
        void setK(Real k) {L_ = ConvertKToL(k);}
        void setSquaredGradients(Mesh& m);
        Real ConvertKToL(Real k) {return k * (2 * gamma_) / rho_ - mu_;}

        Real getrho() {return rho_;}
        Real getmu() {return mu_;}
        Real getgamma() {return gamma_;}

        Real getL2() {return L2_;}
        const std::vector<Real>& getL2_list() {return L2_list_;}
        const std::vector<Real>& getarea_list() {return area_list_;}
        const std::vector<Real>& getvolume_list() {return volume_list_;}
        const std::vector<Real>& getVnbs_list() {return Vnbs_list_;}
        const std::vector<Real>& getAnbs_list() {return Anbs_list_;}

    private:
        // triangle areas 
        std::vector<Real> TriangleAreas_;
        std::vector<std::vector<int>> MapVertexToFace_;
        std::map<INT2, std::vector<int>> MapEdgeToOpposingVerts_;

        int numVerts_, maxStep_=1e5, maxBoundaryStep_=1e5, maxL2step_=1e5;
        Real tol_=1e-5, boundarytol_=1e-5, L2tol_=5e-5;
        Real stepsize_=1e-1, k0_, boundarystepsize_=0.5, L2_stepsize_=10.0f;

        std::vector<vertex> newVertices_;
        std::vector<Real> TotalArea_;

        std::vector<Real3> dAdpi_, dVdpi_;

        Real temperature_=298.0f;
        Real L_=0.0f;
        Real L2_=0.0f;
        Real gamma_; // kJ/nm2
        Real mu_;
        Real rho_=5.0333e-23; // mol/nm3
        Real dgamma_gamma_;
        Real zstar_;
        Real zstar_deviation_=0.003;

        bool debug_=false, useNumerical_=true;

        int print_every=100, optimize_every=1e10, boundary_optimize_every=300;

        std::vector<Real> FE_;

        bool MaxStepCriteria=true, boundaryMaxStepCriteria_=true;

        // list of things to keep track of
        std::vector<Real> volume_list_, area_list_, Vnbs_list_, Anbs_list_, L2_list_;
        Real a_,V_, Vnbs_, Anbs_, Vunderneath_, Vnbs_underneath_;

        Mesh m_flatContact_;

        std::vector<Real3> squared_gradients_;
        Real3 epsilon_={1e-10,1e-10,1e-10};
};
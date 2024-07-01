#pragma once 

#include "MeshRefineStrategy.h"
#include "Curvature.h"
#include "tools/CommonTypes.h"
#include "MeshCurvatureflow.h"
#include "ShortEdgeRemoval.h"
#include "ObtuseTriangleRemoval.h"

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

class CurvatureEvolution : public MeshRefineStrategy
{
    public:
        using curveptr = std::unique_ptr<Curvature>;
        using edgeRemovePtr = std::unique_ptr<ShortEdgeRemoval>;
        using angleRemovePtr= std::unique_ptr<ObtuseTriangleRemoval>;
        using flowptr  = std::unique_ptr<MeshCurvatureflow>;
        using Matrix   = CommonTypes::Matrix;
        using Real2    = CommonTypes::Real2;

        CurvatureEvolution(MeshRefineStrategyInput& input);

        virtual void refine(Mesh& mesh) override;

        void init();
        void update();
        void CleanMesh();
        void setMeanCurvature(Real meanC) {meanCurvature_ = meanC;}

    private:
        std::string curvaturetype_;
        curveptr curvatureCalc_;
        flowptr  curvatureflow_;

        // The indices for the outside vertices to not be kept constant 
        std::vector<int> VertexIndices_;

        std::string fixed_index_file_;
        std::vector<int> fixed_index_;
        std::vector<int> isfixed_;

        Real StepSize_;
        Real meanCurvature_;
        int num_vertices_;

        Real tol_ = 1e-4;

        // keep track of number of iteration
        int iteration_=0;
        int maxStep = 1e6;

        // whether or not we want to print the error 
        Real err_;
        int skip_=1;

        // are we doing fairing 
        bool fairing_=false, cleanMesh_=false;
        std::string fairing_iteration_, fairing_step_;
        Real edgeCutoffLength_;

        // short edge removal
        edgeRemovePtr edgeRemover_;
        angleRemovePtr angleRemover_;

        // keep track of the original positions
        std::vector<Real3> original_pos_;

        // debug
        bool debug_=false;
};
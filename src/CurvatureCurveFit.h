#pragma once
#include "Graph.h"
#include "Curvature.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "Graph.h"
#include "LinAlgTools.h"

#include <vector>
#include <array>

class CurvatureCurveFit : public Curvature
{
    public:
        using Matrix = CommonTypes::Matrix;
        using Real2  = CommonTypes::Real2;
        CurvatureCurveFit(CurvatureInput& input);

        virtual void calculate(Mesh& mesh) override;

        // printing function specifically for curvefit
        void printSecondFundamentalForm(std::string name);
    
    private:
        int NumParameters_ = 3;
        int NumNeighbors_ = 1;

        std::vector<Real3> ff2_parameter_;
        bool extend_boundary_=false;
        int boundary_neighbors_=0;
};
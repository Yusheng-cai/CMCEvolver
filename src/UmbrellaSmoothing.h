#pragma once

#include "MeshRefineStrategy.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "Mesh.h"
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "Eigen/QR"
#include "Eigen/LU"
#include "Eigen/SparseLU"


#include <vector>
#include <array>
#include <string>
#include <iostream>

// umbrella smoothing 
/*
    umbrella smoothing attempts to solve the following differential equation 

    dX/dt = lambda L(X)

    where L is the umbrella operator, namely L(xi) = 1/m \sum_{j=1}^{N1(i)} xj - xi where N1(i) signifies the 
    1-ring neighbors of vertex i

    The iterative scheme is as follow Xn+1 = (I + \lambda dt * L)Xn
*/

class UmbrellaSmoothing : public MeshRefineStrategy
{
    public:
        using triplet = Eigen::Triplet<Real,int>;
        using Sparse_mat = Eigen::SparseMatrix<Real>;
        using INT2    = CommonTypes::index2;

        UmbrellaSmoothing(MeshRefineStrategyInput& input);

        virtual void refine(Mesh& mesh) override;

        void refineStepImplicit();
        void prepareImplicitMatrix();
    
    private:
        int numIterations_;
        
        Real lambdadt_=0.5;

        std::vector<vertex> newVertices_;
        std::vector<vertex> oldVertices_;

        // eigen triplet is a data structure that is useful for sparse matrix storage (i,j,value)
        std::vector<triplet> triplets_;
        Eigen::SparseMatrix<Real> L_;
        std::vector<Real3> Lfactors_;

        // initialize the solver
        Eigen::BiCGSTAB<Sparse_mat> solver_;

        bool scale_ = false;

        Real initialVolume_;
};
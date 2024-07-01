#include "UmbrellaSmoothing.h"

namespace MeshRefineStrategyFactory
{
    registry_<UmbrellaSmoothing> registerUmbrella("umbrella");
}

UmbrellaSmoothing::UmbrellaSmoothing(MeshRefineStrategyInput& input)
:MeshRefineStrategy(input)
{
    input.pack.ReadNumber("iterations", ParameterPack::KeyType::Required, numIterations_);
    input.pack.ReadNumber("lambdadt", ParameterPack::KeyType::Optional, lambdadt_);
    input.pack.Readbool("scale", ParameterPack::KeyType::Optional, scale_);
}

void UmbrellaSmoothing::refine(Mesh& mesh)
{
    mesh_ = &mesh;

    // calculate boundary vertices
    MeshTools::CalculateVertexNeighbors(*mesh_, neighborIndices_);
    MeshTools::MapEdgeToFace(*mesh_, MapEdgeToFace_, MapVertexToEdge_);
    MeshTools::CalculateBoundaryVertices(*mesh_, MapEdgeToFace_,boundaryIndicator_);

    // get the vertices
    const auto& vertices = mesh_->getvertices();

    // fill the old vertices with the original vertices
    oldVertices_.insert(oldVertices_.end(), vertices.begin(), vertices.end());

    // resize the new vertices 
    newVertices_.insert(newVertices_.end(), vertices.begin(), vertices.end());

    // prepare the implicit matrix, it doesn't change for umbrella smoothing
    prepareImplicitMatrix();

    // find volume if needed 
    if ( scale_){
        initialVolume_ = mesh_->calculateVolume();
        std::cout << "Initial volume = " << initialVolume_ << std::endl;
    }

    for (int i=0;i<numIterations_;i++){
        std::cout << "Step " << i << std::endl;
        refineStepImplicit();
    }

    auto& vert = mesh_->accessvertices();
    vert.clear();
    vert.insert(vert.end(), newVertices_.begin(), newVertices_.end());
    mesh_->CalcVertexNormals();
}

void UmbrellaSmoothing::refineStepImplicit()
{
    const auto& vertices = mesh_->getvertices();

    // obtain the rhs of the equation to be solved 
    std::vector<Eigen::VectorXd> rhs;
    std::vector<Eigen::VectorXd> xyz;

    rhs.resize(3, Eigen::VectorXd::Zero(vertices.size()));
    xyz.resize(3, Eigen::VectorXd::Zero(vertices.size()));

    #pragma omp parallel for
    for(int i = 0; i < vertices.size(); ++i)
    {
        rhs[0][i] = Lfactors_[i][0] * lambdadt_ + vertices[i].position_[0];
        rhs[1][i] = Lfactors_[i][1] * lambdadt_ + vertices[i].position_[1];
        rhs[2][i] = Lfactors_[i][2] * lambdadt_ + vertices[i].position_[2];
    }

    for (int i=0;i<3;i++){
        xyz[i] = solver_.solveWithGuess(rhs[i],rhs[i]);
    }

    for (int i=0;i<vertices.size();i++){
        for (int j=0;j<3;j++){
            newVertices_[i].position_[j] = xyz[j][i];
        }
    }

    // update the vertices 
    auto& vert = mesh_->accessvertices();
    vert.clear();
    vert.insert(vert.end(), newVertices_.begin(), newVertices_.end());

    if (scale_){
        Real vol = mesh_->calculateVolume();
        Real scale = std::pow(initialVolume_/vol, 1.0/3.0);
        mesh_->scaleVertices(scale);
    }
}


void UmbrellaSmoothing::prepareImplicitMatrix()
{
    const auto& vertices = mesh_->getvertices();
    L_.resize(vertices.size(), vertices.size());

    // initialize the triplets for EIGEN sparse matrix 
    triplets_.clear();

    // resize Lfactors
    Lfactors_.resize(vertices.size());

    // iterate over vertices 
    #pragma omp parallel
    {
        std::vector<triplet> localtriplets;

        #pragma omp for 
        for (int i=0;i<vertices.size();i++)
        {
            Real3 localLfactor = {{0,0,0}}; 
            if (! MeshTools::IsBoundary(i, boundaryIndicator_))
            {
                int numneighbors = neighborIndices_[i].size();
                Real factor = 1.0/numneighbors;

                for (int j=0;j<numneighbors;j++)
                {
                    localtriplets.push_back(triplet(i, neighborIndices_[i][j], factor));

                    // perform adding the lfactor only if mesh is periodic 
                    if (mesh_->isPeriodic()){
                        Real3 shiftVec;
                        mesh_ -> CalculateShift(vertices[neighborIndices_[i][j]].position_, vertices[i].position_, shiftVec);
                        localLfactor = localLfactor + shiftVec * factor;
                    }
                }
                localtriplets.push_back(triplet(i,i, -1));
            }
            Lfactors_[i] = localLfactor;
        }

        #pragma omp critical
        {
            triplets_.insert(triplets_.end(), localtriplets.begin(), localtriplets.end());
        }
    }

    // set and factor the L matrix
    L_.setFromTriplets(triplets_.begin(), triplets_.end());
    Eigen::SparseMatrix<Real> I = Eigen::MatrixXd::Identity(vertices.size(), vertices.size()).sparseView();
    L_ = I - lambdadt_*L_;
    solver_.compute(L_);
}
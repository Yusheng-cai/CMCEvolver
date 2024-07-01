#include "MeshCurvatureflow.h"

namespace MeshRefineStrategyFactory
{
    registry_<MeshCurvatureflow> registercurveflow("curvatureflow");
}

MeshCurvatureflow::MeshCurvatureflow(MeshRefineStrategyInput& input)
: MeshRefineStrategy(input)
{
    input.pack.ReadNumber("iterations", ParameterPack::KeyType::Required, numIterations_);
    input.pack.ReadNumber("lambdadt", ParameterPack::KeyType::Optional, lambdadt_);
    input.pack.Readbool("scale", ParameterPack::KeyType::Optional, scale_);
    input.pack.Readbool("Decimate", ParameterPack::KeyType::Optional, decimate_);
    input.pack.Readbool("shiftCOM", ParameterPack::KeyType::Optional, shift_COM_);

    pack_.ReadNumber("NumBoundarySmoothing", ParameterPack::KeyType::Optional, StopBoundarySmoothing_);

    // set Eigen to be using the correct number of threads 
    Eigen::initParallel();
    Eigen::setNbThreads(0);
}

void MeshCurvatureflow::updateMesh(){
    // find COM if necessary
    if(shift_COM_){
        COM_ = mesh_->CalculateCOM();
    }

    // obtain map from edge index {minIndex, maxIndex} to the face index 
    // obtain map from vertex index {index} to the Edge Index {minIndex, maxIndex}
    MeshTools::MapEdgeToFace(*mesh_, MapEdgeToFace_, MapVertexToEdge_);

    // Calculate the boundary vertices --> vertices which that has an edge shared by only 1 face
    MeshTools::CalculateBoundaryVertices(*mesh_, MapEdgeToFace_, boundaryIndicator_);

    // Calculate the vertex neighbors
    MeshTools::CalculateVertexNeighbors(*mesh_, neighborIndices_);

    // Map from vertex indices to face indices 
    MeshTools::MapVerticesToFaces(*mesh_, MapVertexToFace_);

    // Map from edges to their opposing vertex indices 
    MeshTools::MapEdgeToOpposingVertices(*mesh_, MapEdgeToFace_, MapEdgeToOpposingVerts_);

    // find number of vertices 
    const auto& v = mesh_->getvertices();
    numVerts_ = v.size();
    newVertices_.clear();newVertices_.resize(numVerts_);
    TotalArea_.clear();TotalArea_.resize(numVerts_);

    // obtain the rhs of the equation to be solved 
    rhs_.clear();xyz_.clear();
    rhs_.resize(3, Eigen::VectorXd::Zero(numVerts_));
    xyz_.resize(3, Eigen::VectorXd::Zero(numVerts_));
}

void MeshCurvatureflow::refine(Mesh& mesh){
    mesh_ = &mesh;

    updateMesh();

    // if we are scaling, we need to calculate the volume
    if (scale_){initialVolume_ = std::abs(mesh_->calculateVolume()); std::cout << "Initial volume = " << initialVolume_ << "\n";}
    

    // start refining 
    for (int i=1;i<=numIterations_;i++){
        // every iteration, decimate the degenerate triangles if necessary
        if (decimate_){
            if (MeshTools::decimateDegenerateTriangle(*mesh_)){
                updateMesh();
            }
        }

        std::cout << "Iteration " << i << std::endl;
        if (i > StopBoundarySmoothing_){
            SmoothBoundary_=false;
        }
        refineImplicitStep();
    }

    mesh_->CalcVertexNormals();
}

bool MeshCurvatureflow::CalculateCotangentWeights(int i, Real3& Lfactor, std::vector<Real>& cotFactors)
{
    // clear the input 
    cotFactors.clear();
    Lfactor.fill(0);

    // get the neighbor indices 
    std::vector<int> neighbors = neighborIndices_[i];
    int neighborSize = neighbors.size();

    // obtain triangles and vertices from mesh
    const auto& triangles = mesh_->gettriangles();
    const auto& vertices  = mesh_->getvertices();

    // iterate over the neighbor indices
    for (int j=0;j<neighborSize;j++){
        // make an edge between this vertex and its neighbor
        INT2 edge = MeshTools::makeEdge(i,neighbors[j]);

        // map this edge to the faces that it corresponds to 
        std::vector<int> faces;
        bool FoundEdge = Algorithm::FindInMap(MapEdgeToFace_, edge, faces);
        ASSERT(FoundEdge, "The edge " << edge <<  " is not found.");

        // map this edge to the opposing vertices indices 
        std::vector<int> OpposingVerts;
        bool FoundVerts = Algorithm::FindInMap(MapEdgeToOpposingVerts_, edge, OpposingVerts);
        ASSERT(FoundVerts, "The edge " << edge << " is not found.");

        // factor here is cot(\alpha) + cot(\beta) --> where \alpha and \beta are the opposing angles shared by an edge 
        Real cotOpposingAnglesSum = 0.0;
        int index1 = edge[0];
        int index2 = edge[1];

        // calculate the cosine theta with respect to opposing points 
        for (int OpposingIdx : OpposingVerts)
        {
            Real3 vec1, vec2;
            Real vec1sq, vec2sq;

            // calculate the vertex distance and vector 
            mesh_->getVertexDistance(vertices[index1], vertices[OpposingIdx], vec1, vec1sq);
            mesh_->getVertexDistance(vertices[index2], vertices[OpposingIdx], vec2, vec2sq);

            // find the costine angle between
            Real costheta  = LinAlg3x3::findCosangle(vec1, vec2);
            Real sin2theta = 1 - costheta*costheta;

            // if sin2theta is too close to 0 or even less than 0 then we just return false 
            if (sin2theta < epsilon_){return false;}

            Real sintheta = std::sqrt(sin2theta);
            cotOpposingAnglesSum += costheta/sintheta;
        }

        if (mesh_->isPeriodic()){
            Real3 ShiftVec;

            // xj - xi 
            mesh_ -> CalculateShift(vertices[neighbors[j]].position_, vertices[i].position_,ShiftVec);

            Lfactor = Lfactor + cotOpposingAnglesSum * ShiftVec;
        }

        cotFactors.push_back(cotOpposingAnglesSum);
    }

    return true;
}

void MeshCurvatureflow::refineImplicitStep()
{
    // obtain the normals of the surfaces 
    mesh_->CalcVertexNormals();

    // Get sum of triangles for each of the vertices  
    std::vector<Real3> normals;
    MeshTools::CalculateTriangleAreasAndFaceNormals(*mesh_, TriangleAreas_, normals);
    std::fill(TotalArea_.begin(), TotalArea_.end(), 0.0);

    #pragma omp parallel for
    for (int i=0;i<MapVertexToFace_.size();i++){
        for (int j=0;j<MapVertexToFace_[i].size();j++){
            TotalArea_[i] += 1.0 / 3.0 * TriangleAreas_[MapVertexToFace_[i][j]];
        }
    }

    // obtain the implicit matrices
    L_  = CalculateImplicitMatrix();

    // populate the solved matrices 
    const auto& vertices = mesh_->getvertices();
    #pragma omp parallel for
    for(int i = 0; i < numVerts_; ++i){
        // we must scale the rhs by Lfactors_ as in (xj - xi + L)
        rhs_[0][i] = lambdadt_ * Lfactors_[i][0] + vertices[i].position_[0];
        rhs_[1][i] = lambdadt_ * Lfactors_[i][1] + vertices[i].position_[1];
        rhs_[2][i] = lambdadt_ * Lfactors_[i][2] + vertices[i].position_[2];
    }

    // solve the sparse system of equations
    Sparse_LU sparsesolver;
    sparsesolver.compute(L_);
    ASSERT((sparsesolver.info() == Eigen::Success), "Compute step failed.");
    for (int i=0;i<3;i++){
        xyz_[i] = sparsesolver.solve(rhs_[i]);
        ASSERT((sparsesolver.info() == Eigen::Success), "The solver failed.");
    }

    // copy over to new vertices data 
    #pragma omp parallel for
    for (int i=0;i<numVerts_;i++){
        for (int j=0;j<3;j++){
            newVertices_[i].position_[j] = xyz_[j][i];
        }
    }

    // update the mesh vertices 
    auto& vert = mesh_->accessvertices();
    vert.clear();
    vert.insert(vert.end(), newVertices_.begin(), newVertices_.end());

    // calculate the volume if needed 
    if (scale_){
        Real vol = std::abs(mesh_->calculateVolume());
        Real scale = std::pow(initialVolume_/vol, 1.0/3.0);

        // this already performs update 
        mesh_->scaleVertices(scale);
        Real newvol = std::abs(mesh_->calculateVolume());

        // shift com
        if (shift_COM_){
            Real3 COM_before = mesh_->CalculateCOM();
            mesh_->ShiftCOMWithRespectTo(COM_);
            Real3 COM_after = mesh_->CalculateCOM();

        }
        std::cout << "new volume = " << newvol << "\n";
    }
}


Eigen::SparseMatrix<MeshCurvatureflow::Real> MeshCurvatureflow::CalculateImplicitMatrix()
{
    const auto& vertices = mesh_->getvertices();

    Lfactors_.clear();
    Lfactors_.resize(vertices.size(), {{0,0,0}});
    triplets_.clear();

    #pragma omp parallel
    {
        std::vector<triplet> triplet_local;
        #pragma omp for
        for (int i=0;i<vertices.size();i++){
            if (! MeshTools::IsBoundary(i, boundaryIndicator_)){
                int numneighbors = neighborIndices_[i].size();

                // get the weights of the neighbor indices 
                Real3 Lfactor = {};
                std::vector<Real> CotWeights;

                // success is if there is not NAN values 
                bool SUCCESS = CalculateCotangentWeights(i, Lfactor, CotWeights);
                Real w = 0.0;

                // Skip the calculation is Total area is 0 --> less than some threshold --> 1e-5
                if (SUCCESS){
                    Real factor = 1.0/(4.0 * TotalArea_[i]);

                    for (int j=0;j<numneighbors;j++){
                        w += CotWeights[j];
                    }

                    for (int j=0;j<numneighbors;j++){
                        triplet_local.push_back(triplet(i, neighborIndices_[i][j], factor * CotWeights[j]));
                    }

                    Lfactors_[i] = Lfactor * factor;
                    triplet_local.push_back(triplet(i,i, -factor * w));
                }
            }
        }

        #pragma omp critical
        {
            triplets_.insert(triplets_.end(), triplet_local.begin(), triplet_local.end());
        }
    }


    // add boundary vertices to the triplets     
    if (SmoothBoundary_){
        std::cout << "Smoothing boundary" << "\n";
        smoothBoundaryVertices();
    }

    // set the L matrix 
    Eigen::SparseMatrix<Real> L;
    L.setZero();
    L.resize(numVerts_, numVerts_);
    L.setFromTriplets(triplets_.begin(), triplets_.end());
    Eigen::SparseMatrix<Real> I = Eigen::MatrixXd::Identity(numVerts_, numVerts_).sparseView();
    L = I - lambdadt_*L;

    return L;
}

void MeshCurvatureflow::smoothBoundaryVertices()
{
    std::vector<int> boundaryIndices;
    const auto& verts = mesh_->getvertices();

    for (int i=0;i<numVerts_;i++){
        if (MeshTools::IsBoundary(i, boundaryIndicator_)){
            boundaryIndices.push_back(i);
        }
    }

    for (int i=0;i<boundaryIndices.size();i++){
        // boundary index 
        int bindex = boundaryIndices[i];

        // find its neighbors 
        std::vector<int> neighbors = neighborIndices_[bindex];

        // the neighbors of this boundary which also happens to be boundaries 
        std::vector<int> bNeighbor;
        for (int n : neighbors){
            if (MeshTools::IsBoundary(n, boundaryIndicator_)){
                bNeighbor.push_back(n);
            }
        }
        ASSERT((bNeighbor.size() == 2), "Boundary vertices can at most be connected to 2 other boundary vertex while it is connected to " \
                                        << bNeighbor.size());

        // smoothing 
        Real E=0;        
        std::vector<Real> factor(bNeighbor.size());
        Real factor_sum=0.0;
        Real3 Lfactor = {};       
        for (int j=0;j<bNeighbor.size();j++){
            int n = bNeighbor[j];
            Real3 distvec;
            Real eij;
            mesh_->getVertexDistance(verts[n].position_, verts[bindex].position_, distvec, eij);
            E += eij;

            if (mesh_->isPeriodic()){
                Real3 shift;
                mesh_->CalculateShift(verts[n].position_, verts[bindex].position_, shift);
                Lfactor = Lfactor + shift * 2.0/eij;
            }
            factor[j] = 2.0/eij;
            factor_sum += factor[j];
        }
        factor = factor / E;
        Lfactor = Lfactor / E;

        for (int j=0;j<bNeighbor.size();j++){
            triplets_.push_back(triplet(bindex, bNeighbor[j], factor[j]));
        }
        triplets_.push_back(triplet(bindex, bindex, -factor_sum/E));

        if (mesh_->isPeriodic()){
            Lfactors_[bindex] = Lfactors_[bindex] + Lfactor;
        }
    }
}
#include "CurvatureEvolution.h"

namespace MeshRefineStrategyFactory
{
    registry_<CurvatureEvolution> registercurve("CurvatureEvolution");
}

CurvatureEvolution::CurvatureEvolution(MeshRefineStrategyInput& input)
: MeshRefineStrategy(input)
{
    // initialize curvature 
    auto cPack = pack_.findParamPack("Curvature", ParameterPack::KeyType::Required);
    cPack->ReadString("type", ParameterPack::KeyType::Required, curvaturetype_);
    CurvatureInput in = {const_cast<ParameterPack&>(*cPack)};
    curvatureCalc_ = curveptr(CurvatureRegistry::Factory::instance().create(curvaturetype_,in));

    // read in the necessary ones for update 
    pack_.ReadNumber("stepsize", ParameterPack::KeyType::Required, StepSize_);
    pack_.ReadNumber("tolerance", ParameterPack::KeyType::Optional, tol_);
    pack_.ReadNumber("k0", ParameterPack::KeyType::Required, meanCurvature_);
    pack_.ReadNumber("maxstep", ParameterPack::KeyType::Optional, maxStep);
    pack_.Readbool("fairing", ParameterPack::KeyType::Optional, fairing_);
    pack_.Readbool("cleanMesh", ParameterPack::KeyType::Optional, cleanMesh_);
    pack_.Readbool("debug", ParameterPack::KeyType::Optional, debug_);

    fixed_index_.clear();
    if (pack_.ReadString("fixed_index_file", ParameterPack::KeyType::Optional, fixed_index_file_)){
        StringTools::ReadTabulatedData(fixed_index_file_, 0, fixed_index_);
    }

    // whether or not we are doing fairing 
    if (fairing_){
        pack_.ReadString("fairing_iteration", ParameterPack::KeyType::Required, fairing_iteration_);
        pack_.ReadString("fairing_step", ParameterPack::KeyType::Required, fairing_step_);

        ParameterPack fairing_pack;
        fairing_pack.insert("iterations", fairing_iteration_);
        fairing_pack.insert("lambdadt", fairing_step_);
        fairing_pack.insert("name", "f");
        fairing_pack.insert("Decimate", "false");
        MeshRefineStrategyInput input = {fairing_pack};

        curvatureflow_ = flowptr(new MeshCurvatureflow(input));
    }

    // whether or not we are cleaning 
    if (cleanMesh_){
        pack_.ReadNumber("edgeLengthCutoff", ParameterPack::KeyType::Required, edgeCutoffLength_);
    }
}


void CurvatureEvolution::init(){
    const auto& vertices = mesh_->getvertices();
    num_vertices_ = vertices.size();
    original_pos_.clear();
    original_pos_.resize(num_vertices_);
    for (int i=0;i<original_pos_.size();i++){
        original_pos_[i] = vertices[i].position_;
    }

    // find the boundary vertices
    MeshTools::MapEdgeToFace(*mesh_, MapEdgeToFace_, MapVertexToEdge_);
    MeshTools::CalculateBoundaryVertices(*mesh_, MapEdgeToFace_, boundaryIndicator_);

    isfixed_.clear();
    isfixed_.resize(vertices.size(),0);
    for (int i=0;i<fixed_index_.size();i++){
        int ind = fixed_index_[i];
        isfixed_[ind] = 1;
    }

    VertexIndices_.clear();
    for (int i=0;i<vertices.size();i++){
        if (! MeshTools::IsBoundary(i, boundaryIndicator_)){
            VertexIndices_.push_back(i);
        }
    }
}

void CurvatureEvolution::update(){
    const auto& vertices = mesh_->getvertices();

    // find the boundary vertices
    MeshTools::MapEdgeToFace(*mesh_, MapEdgeToFace_, MapVertexToEdge_);
    MeshTools::CalculateBoundaryVertices(*mesh_, MapEdgeToFace_, boundaryIndicator_);

    VertexIndices_.clear();
    for (int i=0;i<vertices.size();i++){
        if (! MeshTools::IsBoundary(i, boundaryIndicator_)){
            VertexIndices_.push_back(i);
        }
    }
}

void CurvatureEvolution::CleanMesh(){
    // set boundary or fixed vertices to be of high importance
    std::vector<int> importance(mesh_->getNumVertices(),0);
    for (int i=0;i<importance.size();i++){
        if ((MeshTools::IsBoundary(i, boundaryIndicator_))){
            importance[i] = -1;
        }
        // -1 importance means fixed
        if (isfixed_[i]){
            importance[i] = -1;
        }
    }

    // set the importance of the edges
    edgeRemover_ = edgeRemovePtr(new ShortEdgeRemoval(*mesh_));
    edgeRemover_->set_importance(importance);

    // cleant he mesh by collapsing edges 
    int num_collapsed = edgeRemover_->calculate(edgeCutoffLength_);
    
    std::cout << "num collapse = " << num_collapsed << "\n";

    if (num_collapsed != 0){
        MeshTools::MapEdgeToFace(*mesh_, MapEdgeToFace_, MapVertexToEdge_, false);
        MeshTools::CalculateBoundaryVertices(*mesh_, MapEdgeToFace_, boundaryIndicator_);
        int numB=0;
        for (int i=0;i<boundaryIndicator_.size();i++){
            if (boundaryIndicator_[i]){
                numB += 1;
            }
        }
        MeshTools::RemoveIsolatedFaces(*mesh_);
        MeshTools::RemoveDuplicatedFaces(*mesh_);
        std::vector<int> MapOldToNew = MeshTools::RemoveIsolatedVertices(*mesh_);

        MeshTools::MapEdgeToFace(*mesh_, MapEdgeToFace_, MapVertexToEdge_);
        MeshTools::CalculateBoundaryVertices(*mesh_, MapEdgeToFace_, boundaryIndicator_);
        int numBafter=0;
        for (int i=0;i<boundaryIndicator_.size();i++){
            if (boundaryIndicator_[i]){
                numBafter += 1;
            }
        }

        ASSERT((numB==numBafter), "No new boundary should be created.");

        isfixed_.clear();
        isfixed_.resize(mesh_->getNumVertices(),0);
        std::vector<int> newfixed_index;
        for (int i=0;i<fixed_index_.size();i++){
            int ind = MapOldToNew[fixed_index_[i]];
            if (ind != -100){
                newfixed_index.push_back(ind);
                isfixed_[ind] = 1;
            }
        }
        fixed_index_ = newfixed_index;
        std::cout << "Fixed index = " << newfixed_index.size() << "\n";

        update();
    }
}

void CurvatureEvolution::refine(Mesh& mesh){
    // get a copy of the mesh (pointer)
    mesh_ = &mesh;
    mesh_->CalcVertexNormals();

    // find the vertices to be updated
    init();

    // Start refining , initialize err and iteration
    err_ = std::numeric_limits<Real>::max();
    iteration_=1;
    Real last_max_err=std::numeric_limits<Real>::max();

    while (err_ >= tol_){
        // set max error numeric min
        Real maxerr = -std::numeric_limits<Real>::max();
        int maxIndex;

        // first let's calculate the curvatures 
        curvatureCalc_ -> calculate(*mesh_);

        // then let's obtain the curvatures for each vertices 
        const auto& curvatures = curvatureCalc_ -> getAvgCurvaturePerVertex();

        // get vertics 
        auto& vertices = mesh_->accessvertices();

        // update the vertices 
        Real avgE=0.0;
        #pragma omp parallel
        {
            Real e=-std::numeric_limits<Real>::max();
            int max_I;
            Real avgElocal = 0.0;

            #pragma omp for
            for (int j=0;j<VertexIndices_.size();j++){
                int index = VertexIndices_[j];

                // find the difference in curvature 
                Real diffkappa = meanCurvature_ - curvatures[index];

                // add to local average error
                avgElocal += std::abs(diffkappa);

                if (std::abs(diffkappa) > e){e = std::abs(diffkappa); max_I=j;}

                // update the vertex positions 
                vertices[index].position_ = vertices[index].position_ + StepSize_ * vertices[index].normals_ * diffkappa;
            }

            #pragma omp critical
            if (e > maxerr)
            {
                maxerr = e;
                maxIndex = max_I;
            }

            #pragma omp critical
            {
                avgE += avgElocal;
            }

        }

        avgE = avgE / VertexIndices_.size();
        err_ = std::abs(avgE);

        // if ((err_- last_max_err) > 1e-5){
        //     StepSize_ /= 2;
        // }

        last_max_err = err_;

        if (StepSize_ < 1e-5){
            std::cout << "Stepsize has become too small, stepsize = " << StepSize_ << ". Exiting...." << "\n";
            break;
        }

        std::cout << "Max error is " << maxerr << std::endl;
        std::cout << "average error is " << err_ << "\n";

        // update the normals as well 
        mesh_->CalcVertexNormals();

        // clean mesh before we perform the calculation
        if (cleanMesh_){
            CleanMesh();
        }

        // fair the mesh if necessary
        if (fairing_){
            if (err_ > 10){
                curvatureflow_->refine(*mesh_);
            }
        }

        // break the calculation if iteration is already maxed out
        if (iteration_ > maxStep){
            std::cout << "Evolution finished premature at iteration " << iteration_-1 << "\n";
            break;
        }

        if (debug_){
            std::string ite = StringTools::TypeToString(iteration_);
            MeshTools::writePLY("Iteration_" + ite +".ply", *mesh_);
            std::ofstream ofs;
            ofs.open("Curvature_" + ite + ".out");
            for (int i=0;i<curvatures.size();i++){
                ofs << curvatures[i] << "\n";
            }
            ofs.close();
        }

        // update the iterations
        iteration_ ++;
    }
}
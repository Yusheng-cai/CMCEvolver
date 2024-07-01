#include "CurvatureFDM.h"

namespace CurvatureRegistry
{
    registry<CurvatureFDM> registerFDM("FDM");
}

CurvatureFDM::CurvatureFDM(CurvatureInput& input)
:Curvature(input)
{
}

void CurvatureFDM::calculate(Mesh& mesh)
{
    initialize(mesh);

    MeshTools::CalculateVertexNeighbors(mesh, neighbor_indices_);

    // vertices 
    const auto& vertices = mesh.getvertices();
    int nv = vertices.size();

    // Now we calculate the curvature
    #pragma omp parallel for
    for (int i=0;i<nv;i++){
        auto& v1 = vertices[i];
        auto& neighbors = neighbor_indices_[i];
        int numneighbors= neighbors.size();
        Real GaussCurve = 1.0;

        for (int n : neighbors){ 
            auto& v2 = vertices[n];

            Real curve_Val;
            Real3 diff = v2.position_ - v1.position_;
            Real3 diffn = v2.normals_ - v1.normals_;


            Real diffnnorm = LinAlg3x3::norm(diffn); 
            Real sq_diff   = LinAlg3x3::DotProduct(diff, diff);
            Real NDotPos   = LinAlg3x3::DotProduct(diffn, diff); 

            curve_Val = NDotPos/sq_diff;
            avgCurvaturePerVertex_[i] += std::abs(curve_Val);
            GaussCurve *= std::abs(curve_Val);
        }
        GaussCurvaturePerVertex_[i] = GaussCurve;
    }

    #pragma omp parallel for
    for (int i=0;i<nv;i++){ 
        int neighborSize = neighbor_indices_[i].size();
        avgCurvaturePerVertex_[i] = avgCurvaturePerVertex_[i]/neighborSize;
        GaussCurvaturePerVertex_[i] = std::pow(GaussCurvaturePerVertex_[i], 1.0/neighborSize);
    } 

    CalculateFaceCurvature(mesh, avgCurvaturePerVertex_, GaussCurvaturePerVertex_, CurvaturePerVertex_,FaceCurvature_);
}

void CurvatureFDM::printCurvature(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# Average Gaussian\n";
    for (int i=0;i<avgCurvaturePerVertex_.size();i++){
        ofs << avgCurvaturePerVertex_[i] << " " << GaussCurvaturePerVertex_[i] << "\n";
    }
    ofs.close();
}
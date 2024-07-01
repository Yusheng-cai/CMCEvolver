#include "CurvatureJetFit.h"

namespace CurvatureRegistry
{
    registry<CurvatureJetFit> registerJet("jetfit");
}

CurvatureJetFit::CurvatureJetFit(CurvatureInput& input)
:Curvature(input)
{
    foundnumneighrs_ = input.pack.ReadNumber("neighbors", ParameterPack::KeyType::Optional, numneighbors_);
    input.pack.ReadNumber("degree", ParameterPack::KeyType::Optional, degree_);
    input.pack.ReadNumber("MongeCoefficient", ParameterPack::KeyType::Optional, MongeCoefficient_);

    outputs_.registerOutputFunc("neighborIndices", [this](std::string name) -> void {this -> printNeighbors(name);});
    outputs_.registerOutputFunc("coefficients", [this](std::string name) -> void {this -> printCoefficientPerVertex(name);});
    outputs_.registerOutputFunc("PCAeigenvector", [this](std::string name) -> void {this -> printPCAeigenvector(name);});

    numpoints_ = (degree_ + 1)*(degree_+2)/2;
}

void CurvatureJetFit::printPCAeigenvector(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    for (int i=0;i<PCAeigenvector_.size();i++){
        for (int j=0;j<3;j++){
            for (int k=0;k<3;k++){
                ofs_ << PCAeigenvector_[i][j][k] << " ";
            }
        }
        ofs_ << "\n";
    }
    ofs_.close();
}

void CurvatureJetFit::calculate(Mesh& mesh)
{
    initialize(mesh);

    // First calculate the vertex neighbors for each of the vertex -->
    // simply just the triangles that each vertex is connected to
    std::vector<std::vector<int>> VertexNeighbors;
    MeshTools::CalculateVertexNeighbors(mesh, VertexNeighbors);

    // start calculating number of neighbors 
    const auto& vertices = mesh.getvertices();
    int nv = vertices.size();

    if (foundnumneighrs_){
        std::cout << "Using number of vertices." << std::endl;
        Graph::BFS_kring_neighbor(VertexNeighbors, numneighbors_, NeighborIndicesNVertex_);
    }
    else{
        std::cout << "Using number of neighbors." << std::endl;
        Graph::getNNearbyIndices(VertexNeighbors,numpoints_-1,NeighborIndicesNVertex_);
    }

    coefficientsPerVertex_.resize(nv);
    PCAeigenvector_.resize(nv);

    #pragma omp parallel
    {
        MongeViaJetFitting jetfitterLocal;

        #pragma omp for
        for (int i=0;i<nv;i++){
            std::vector<Dpoint> vec;
            Real3 thisPos = vertices[i].position_;
            Real3 thisNormal = vertices[i].normals_;
            Dpoint point(thisPos[0], thisPos[1], thisPos[2]);
            vec.push_back(point);

            for (int j=0;j<NeighborIndicesNVertex_[i].size();j++){
                int neighborIndex = NeighborIndicesNVertex_[i][j];
                Real3 vertpos;

                vertpos = mesh.getShiftedVertexPosition(vertices[neighborIndex], vertices[i]);

                Dpoint point(vertpos[0], vertpos[1], vertpos[2]);
                vec.push_back(point);
            }

            ASSERT((vec.size() >= numpoints_), "The neighbors for indices " << i << " " << vec.size() << " is not enough for fitting for degree " << degree_ <<\
            " which has to be at least " << numpoints_);

            auto mform = jetfitterLocal(vec.begin(), vec.end(), degree_, MongeCoefficient_);
            for (int j=0;j<3;j++){
                for (int k=0;k<3;k++){
                    PCAeigenvector_[i][j][k] = jetfitterLocal.pca_basis(j).second[k];
                }
            }

            DVector norm(-thisNormal[0], -thisNormal[1], -thisNormal[2]);

            // fix normal
            mform.comply_wrt_given_normal(norm);

            std::vector<Real> coeff;
            for (int i=0;i<mform.coefficients().size();i++){
                coeff.push_back(mform.coefficients()[i]);
            }

            coefficientsPerVertex_[i] = coeff;
            CurvaturePerVertex_[i][0] = mform.principal_curvatures(0);
            CurvaturePerVertex_[i][1] = mform.principal_curvatures(1);

            DVector vec1 = mform.maximal_principal_direction();
            DVector vec2 = mform.minimal_principal_direction();
            Real3 v;
            Real3 v2;

            for (int j=0;j<3;j++){
                v[j] = vec1[j];
                v2[j] = vec2[j];
            }
            principalDir1_[i] = v;
            principalDir2_[i] = v2;
        }

        #pragma omp for
        for (int i=0;i<CurvaturePerVertex_.size();i++)
        {
            Real avg=0.0;
            Real gauss = 1.0;
            for (int j=0;j<2;j++){
                gauss *= CurvaturePerVertex_[i][j];
                avg   += CurvaturePerVertex_[i][j]/2.0;
            }
            avgCurvaturePerVertex_[i] = avg;
            GaussCurvaturePerVertex_[i] = gauss;
        }
    }

    CalculateFaceCurvature(mesh, avgCurvaturePerVertex_, GaussCurvaturePerVertex_, CurvaturePerVertex_, FaceCurvature_);
}

void CurvatureJetFit::printNeighbors(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    for (int i=0;i<NeighborIndicesNVertex_.size();i++){
        int neighborsize = NeighborIndicesNVertex_[i].size();
        // print self first 
        ofs_ << i << " ";

        for (int j=0;j<neighborsize;j++){
            ofs_ << NeighborIndicesNVertex_[i][j] << " ";
        }
        ofs_ << "\n";
    }
    ofs_.close();
}

void CurvatureJetFit::printCoefficientPerVertex(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    ofs_ << "# ";
    for (int i=0;i<coefficientsPerVertex_[0].size();i++){
        ofs_ << "coefficient" << i+1 << " ";
    }
    ofs_ << "\n";
   

    for (int i=0;i<coefficientsPerVertex_.size();i++){
        int sizePerVertex = coefficientsPerVertex_[i].size();
        for (int j=0;j<sizePerVertex;j++){
            if (j != sizePerVertex-1){
                ofs_ << coefficientsPerVertex_[i][j] << " ";
            }
            else{
                ofs_ << coefficientsPerVertex_[i][j] << "\n";
            }
        }
    }

    ofs_.close();
}
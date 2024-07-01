#include "CurvatureCurveFit.h"

namespace CurvatureRegistry
{
    registry<CurvatureCurveFit> registerCurvefit("curvefit");
}

CurvatureCurveFit::CurvatureCurveFit(CurvatureInput& input)
:Curvature(input){
    input.pack.ReadNumber("neighbors", ParameterPack::KeyType::Optional, NumNeighbors_);
    input.pack.Readbool("extend_boundary", ParameterPack::KeyType::Optional, extend_boundary_);
    if (extend_boundary_){
        input.pack.ReadNumber("boundary_neighbors", ParameterPack::KeyType::Optional, boundary_neighbors_);
    }
}

void CurvatureCurveFit::calculate(Mesh& mesh)
{
    // initialize the mesh 
    initialize(mesh);

    const auto& vertices = mesh.getvertices();
    int nv = vertices.size();

    // resize fundamental form vector which contains (L,M,N)
    ff2_parameter_.resize(nv);

    // calculate vertex neighbors
    std::vector<std::vector<int>> VertexNeighbors;
    MeshTools::CalculateVertexNeighbors(mesh, VertexNeighbors);

    // map edge to face 
    std::map<INT2, std::vector<int>> EToF;
    std::vector<bool> boundary_indicator;
    MeshTools::MapEdgeToFace(mesh, EToF);
    MeshTools::CalculateBoundaryVertices(mesh, EToF, boundary_indicator);

    // calculate neighbor vertices that are N vertex away 
    std::vector<std::vector<int>> NeighborIndicesNVertex;
    Graph::getNearbyIndicesNVertexAway(VertexNeighbors, NumNeighbors_, NeighborIndicesNVertex);
    // if (! extend_boundary_){Graph::BFS_kring_neighbor(VertexNeighbors, NumNeighbors_, NeighborIndicesNVertex);}
    // else{Graph::BFS_kring_neighbor_boundary(VertexNeighbors, boundary_indicator, NumNeighbors_, boundary_neighbors_, NeighborIndicesNVertex);}

    // the reference direction of the normal vector is the z vector
    Real3 referenceDir = {{0,0,1}};

    #pragma omp parallel for
    for (int i=0;i<nv;i++){
        auto& v = vertices[i];
        auto& neighbors = NeighborIndicesNVertex[i];

        if (boundary_indicator[i]){
            int inc_neighbors = boundary_neighbors_ + 1;
            while (neighbors.size() < 3){
                Graph::BFS_kring_neighbor(VertexNeighbors, inc_neighbors, i, neighbors);
                inc_neighbors += 1;
            }
        }
        else{
            ASSERT((neighbors.size()>=3), "The number of points to be fit must be larger or equal to 3, \
            index = " << i << " number of neighbors = " << neighbors.size());
        }

        // How to rotate normal vector to be z vector
        //ASSERT((std::abs(LinAlg3x3::DotProduct(v.normals_ , referenceDir) + 1) > 1e-5), "Normals too close to reference");
        Matrix rotationMat = LinAlg3x3::GetRotationMatrix(v.normals_, referenceDir); 
        Matrix revrotMat   = LinAlg3x3::GetRotationMatrix(referenceDir, v.normals_);
        Real3 vector       = LinAlg3x3::MatrixDotVector(rotationMat, v.normals_);

        // initialize the 
        Eigen::Matrix3d mat = Eigen::Matrix3d::Zero();
        Eigen::Vector3d b   = Eigen::Vector3d::Zero();

        // Fitting equation to the form of b1*x^2 + b2 * y^2 + b3 * x * y
        for (int j=0;j<neighbors.size();j++){
            int neighborId = neighbors[j];
            Real3 diff;
            Real diffsq;

            mesh.getVertexDistance(vertices[neighborId], v, diff, diffsq);

            // Rotate neighbor position to the frame of reference of interest
            Real3 neighborRotatedPos = LinAlg3x3::MatrixDotVector(rotationMat, diff);
            Real xpos = neighborRotatedPos[0];
            Real ypos = neighborRotatedPos[1];
            Real zpos = neighborRotatedPos[2];

            mat(0,0) += std::pow(xpos, 4.0);
            mat(0,1) += std::pow(xpos, 2.0) * std::pow(ypos, 2.0);
            mat(0,2) += std::pow(xpos, 3.0) * ypos;
            mat(1,1) += std::pow(ypos, 4.0);
            mat(1,2) += xpos * std::pow(ypos, 3.0);

            b[0] += std::pow(xpos, 2.0) * zpos;
            b[1] += std::pow(ypos, 2.0) * zpos;
            b[2] += xpos * ypos * zpos;
        }

        mat(1,0) = mat(0,1);
        mat(2,0) = mat(0,2);
        mat(2,1) = mat(1,2);
        mat(2,2) = mat(0,1);


        // Solve the least squares fitting problem 
        Eigen::Vector3d ans = mat.bdcSvd(Eigen::ComputeFullU|Eigen::ComputeFullV).solve(b);
        Eigen::EigenSolver<Eigen::Matrix2d> eigensolver;

        // Ans = [ b1, b2, b12 ]
        Eigen::Matrix2d SecondFundamentalMat;
        SecondFundamentalMat(0,0) = -2*ans[0];
        SecondFundamentalMat(0,1) = ans[2];
        SecondFundamentalMat(1,0) = ans[2];
        SecondFundamentalMat(1,1) = -2*ans[1];
        eigensolver.compute(SecondFundamentalMat);

        // ff2
        Real3 ff2;
        for (int j=0;j<3;j++){
            ff2[j] = ans[j];
        }
        ff2_parameter_[i] = {{2*ff2[0], ff2[1], 2*ff2[2]}};

        // Compute the eigenvalues 
        Eigen::Vector2d eigenvalues = eigensolver.eigenvalues().real();
        Eigen::Matrix2d eigenvectors= eigensolver.eigenvectors().real();

        // sort the eigenvalues from greater to smaller 
        Real2 eigenvals;
        eigenvals[0] = eigenvalues[0];
        eigenvals[1] = eigenvalues[1];
        bool Switch=false;

        if (eigenvals[0] < eigenvals[1]){
            Switch=true;
        }

        std::sort(eigenvals.begin(), eigenvals.end(), std::greater<Real>());
        CurvaturePerVertex_[i] = eigenvals;

        // rotate the eigenvectors 
        Real3 eigvec1;
        eigvec1.fill(0);
        Real3 eigvec2;
        eigvec2.fill(0);

        eigvec1[0] = eigenvectors(0,0);
        eigvec1[1] = eigenvectors(1,0);
        eigvec2[0] = eigenvectors(0,1);
        eigvec2[1] = eigenvectors(1,1);

        std::array<Real3, 2> eigenvectorPair;

        Real3 eigvec1Rot = LinAlg3x3::MatrixDotVector(revrotMat, eigvec1);
        Real3 eigvec2Rot = LinAlg3x3::MatrixDotVector(revrotMat, eigvec2);

        // eigenvectors in real space
        if (Switch){
            principalDir1_[i] = eigvec2Rot;
            principalDir2_[i] = eigvec1Rot;
        }
        else{
            principalDir1_[i] = eigvec1Rot;
            principalDir2_[i] = eigvec2Rot;
        }
    }

    for (int i=0;i<CurvaturePerVertex_.size();i++){
        Real avg = 0.0;
        Real gauss = 1.0;
        for (int j=0;j<2;j++){
            avg += CurvaturePerVertex_[i][j]/2.0;
            gauss *= CurvaturePerVertex_[i][j];
        }

        avgCurvaturePerVertex_[i] = avg;
        GaussCurvaturePerVertex_[i] = gauss;
    }

    CalculateFaceCurvature(mesh, avgCurvaturePerVertex_, GaussCurvaturePerVertex_, CurvaturePerVertex_, FaceCurvature_);
}

void CurvatureCurveFit::printSecondFundamentalForm(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# VertexNumber L M N\n";

    for (int i=0;i<ff2_parameter_.size();i++){
        ofs << i << " " << ff2_parameter_[i][0] << " " << ff2_parameter_[i][1] << " " << ff2_parameter_[i][2] << "\n";
    }

    ofs.close();
}
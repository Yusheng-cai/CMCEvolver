#include "CurvatureTensor.h"

namespace CurvatureRegistry
{
   registry<CurvatureTensor> registerCurvatureTensor("tensor"); 
}

CurvatureTensor::CurvatureTensor(CurvatureInput& input)
:Curvature(input)
{
    outputs_.registerOutputFunc("FF2", [this](std::string name) -> void {this -> printFF2(name);});
}

void CurvatureTensor::calculate(Mesh& mesh)
{
    mesh.CalcPerVertexDir();
    mesh.CalculateCornerArea();

    const auto& triangles = mesh.gettriangles();
    const auto& vertices  = mesh.getvertices();
    const auto& pervertexdir1 = mesh.getPerVertexDir1();
    const auto& pervertexdir2 = mesh.getPerVertexDir2();
    const auto& triangleArea = mesh.getTriangleArea();
    const auto& cornerArea   = mesh.getCornerAreas();

    curvatureTensorPerTriangle_.resize(triangles.size());

    // Fill the curvature Tensor Per vertex with zero arrays
    curvatureTensorPerVertex_.resize(vertices.size(), {});

    TotalAreaPerVertex_.resize(vertices.size(), 0.0);
    CurvaturePerVertex_tot.resize(vertices.size());

    for (int i=0;i<triangles.size();i++)
    {
        auto& t = triangles[i].triangleindices_; 
        std::array<Real3,3> edges;
        Real edge0sq, edge1sq, edge2sq;

        mesh.getVertexDistance(vertices[t[2]], vertices[t[1]], edges[0], edge0sq);
        mesh.getVertexDistance(vertices[t[0]], vertices[t[2]], edges[1], edge1sq);
        mesh.getVertexDistance(vertices[t[1]], vertices[t[0]], edges[2], edge2sq);

        Real3 U = edges[0];
        LinAlg3x3::normalize(U);
        Real3 N = LinAlg3x3::CrossProduct(edges[0], edges[1]);
        LinAlg3x3::normalize(N);
        Real3 V = LinAlg3x3::CrossProduct(N,U);
        LinAlg3x3::normalize(V); 

        // initialize A and B matrix to be solved
        Eigen::Matrix3d A;
        Eigen::Vector3d b;
        b.fill(0);
        A.fill(0);

        std::vector<std::array<Real,2>> eJVec;
        std::vector<std::array<Real,2>> nJVec;
 
        for (int j=0;j<3;j++)
        {
            Real ejU = LinAlg3x3::DotProduct(edges[j], U);
            Real ejV = LinAlg3x3::DotProduct(edges[j], V); 
            eJVec.push_back({{ejU, ejV}});

            A(0,0) += ejU*ejU;
            A(0,1) += ejU*ejV;
            A(2,2) += ejV*ejV;

            // convert ids.
            int id1 = j - 1;
            if ( id1 < 0){id1 += 3;}
            int id2 = j - 2;
            if ( id2 < 0){id2 += 3;}

            Real3 diffN = vertices[t[id1]].normals_ - vertices[t[id2]].normals_;
            Real diffNU = LinAlg3x3::DotProduct(diffN, U);
            Real diffNV = LinAlg3x3::DotProduct(diffN, V);

            nJVec.push_back({{diffNU, diffNV}});
            
            b[0] += diffNU*ejU;
            b[1] += diffNU*ejV + diffNV*ejU;
            b[2] += diffNV*ejV;
        }
        A(1,1) = A(0,0) + A(2,2);
        A(1,2) = A(0,1);
        A(2,1) = A(0,1);
        A(1,0) = A(0,1);

        Eigen::EigenSolver<Eigen::Matrix2d> eigensolver;
        Eigen::Vector3d soln = A.bdcSvd(Eigen::ComputeFullU|Eigen::ComputeFullV).solve(b);
        Real3 ans;
        Eigen::Matrix2d triangleMat;
        triangleMat(0,0) = soln[0];
        triangleMat(0,1) = soln[1];
        triangleMat(1,0) = soln[1];
        triangleMat(1,1) = soln[2];


        eigensolver.compute(triangleMat);

        for (int j=0;j<3;j++)
        {
            ans[j] = soln[j];
        }

        Eigen::Vector2d eig = eigensolver.eigenvalues().real();

        curvatureTensorPerTriangle_[i] = ans;
        
        for (int j=0;j<3;j++)
        {
            int id_ = t[j];
            Real3 u = pervertexdir1[id_];
            Real3 v = pervertexdir2[id_];
            Real3 newcurv = projectCurvature(u,v, U, V,ans);

            Eigen::Matrix2d tempMat;
            Eigen::EigenSolver<Eigen::Matrix2d> eigensolver;
            tempMat(0,0) = newcurv[0];
            tempMat(0,1) = newcurv[1];
            tempMat(1,0) = newcurv[1];
            tempMat(1,1) = newcurv[2];

            eigensolver.compute(tempMat);
            Eigen::Vector2d eig = eigensolver.eigenvalues().real();

            CurvaturePerVertex_tot[id_].push_back(newcurv);
            
            for (int k = 0;k<3;k++)
            {
                curvatureTensorPerVertex_[id_][k] += newcurv[k] * cornerArea[i][j];
            }
            TotalAreaPerVertex_[id_] += triangleArea[i];
        }
    } 

    calculatePrincipalCurvatures(mesh);

    //CalculateFaceCurvature(mesh, avgCurvaturePerVertex_, FaceCurvature_);
}

void CurvatureTensor::calculatePrincipalCurvatures(Mesh& mesh)
{
    initialize(mesh);

    const auto& vertices  = mesh.getvertices();
    const auto& pointArea = mesh.getTriangleArea();

    // Use z as a reference direection
    Real3 referenceDir = {{0,0,1}};
    CurvaturePerVertex_.resize(vertices.size());
    avgCurvaturePerVertex_.resize(vertices.size() ,0.0);
    GaussCurvaturePerVertex_.resize(vertices.size(),0.0);

    principalDir1_.resize(vertices.size());
    principalDir2_.resize(vertices.size());

    for (int i=0;i<vertices.size();i++)
    {
        // find rotation matrix that rotates the normal from [0,0,1] to itself 
        auto& normal = vertices[i].normals_;
        Matrix rotationMat = LinAlg3x3::GetRotationMatrix(referenceDir, normal);

        for (int j=0;j<3;j++)
        {
            curvatureTensorPerVertex_[i][j] = curvatureTensorPerVertex_[i][j]/pointArea[i];
            #ifdef MY_DEBUG
            std::cout << "curvatureTensor = " << curvatureTensorPerVertex_[i][j] << std::endl;
            std::cout << "pointArea " << i << " = " << pointArea[i] << std::endl;
            #endif 
        }

        Eigen::Matrix2d mat;
        Eigen::EigenSolver<Eigen::Matrix2d> eigensolver;

        mat(0,0) = curvatureTensorPerVertex_[i][0];
        mat(0,1) = curvatureTensorPerVertex_[i][1];
        mat(1,0) = curvatureTensorPerVertex_[i][1];
        mat(1,1) = curvatureTensorPerVertex_[i][2];

        eigensolver.compute(mat);
        Eigen::Vector2d eigenvalues = eigensolver.eigenvalues().real();

        if (eigenvalues[0] > eigenvalues[1])
        {
            CurvaturePerVertex_[i][0] = eigenvalues[0];
            CurvaturePerVertex_[i][1] = eigenvalues[1];
        }
        else
        {
            CurvaturePerVertex_[i][0] = eigenvalues[1];
            CurvaturePerVertex_[i][1] = eigenvalues[0];
        }

        Eigen::Matrix2d eigenvectors = eigensolver.eigenvectors().real();

        Real3 eigvec1;
        eigvec1.fill(0);
        Real3 eigvec2;
        eigvec2.fill(0);

        eigvec1[0] = eigenvectors(0,0);
        eigvec1[1] = eigenvectors(1,0);
        eigvec2[0] = eigenvectors(0,1);
        eigvec2[1] = eigenvectors(1,1);

        std::array<Real3, 2> eigenvectorPair;

        Real3 eigvec1Rot = LinAlg3x3::MatrixDotVector(rotationMat, eigvec1);
        Real3 eigvec2Rot = LinAlg3x3::MatrixDotVector(rotationMat, eigvec2);

        if (eigenvalues[0] >= eigenvalues[1])
        {
            principalDir1_[i] = eigvec1Rot;
            principalDir2_[i] = eigvec2Rot;
        }
        else
        {
            principalDir2_[i] = eigvec1Rot;
            principalDir1_[i] = eigvec2Rot;
        }
    }

    for (int i=0;i<CurvaturePerVertex_.size();i++)
    {
        Real avg=0.0;
        Real gauss=1.0;
        for (int j=0;j<2;j++)
        {
            avg += CurvaturePerVertex_[i][j]/2.0;
            gauss *= CurvaturePerVertex_[i][j];
        }
        avgCurvaturePerVertex_[i] = avg;
        GaussCurvaturePerVertex_[i] = gauss;
    }
}

void CurvatureTensor::printFF2(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    for (int i=0;i<curvatureTensorPerVertex_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs_ << curvatureTensorPerVertex_[i][j] << " ";
        }
        ofs_ << "\n";
    }
    ofs_.close();
}

CurvatureTensor::Real3 CurvatureTensor::projectCurvature(const Real3& oldu, const Real3& oldv, const Real3& refu, const Real3& refv,const Real3& curvature)
{
    Real3 newu;
    Real3 newv;
 

    Real3 oldN = LinAlg3x3::CrossProduct(oldu, oldv);
    LinAlg3x3::normalize(oldN);


    Real3 refN = LinAlg3x3::CrossProduct(refu, refv);
    LinAlg3x3::normalize(refN);

    LinAlg3x3::RotateBasisSet(oldN, refN, oldu, oldv, newu, newv);

	Real u1 = LinAlg3x3::DotProduct(newu, refu);
    Real v1 = LinAlg3x3::DotProduct(newu, refv);
    Real u2 = LinAlg3x3::DotProduct(newv, refu);
    Real v2 = LinAlg3x3::DotProduct(newv, refv);

    Real3 newcurvature_;
	newcurvature_[0]  = curvature[0] * u1*u1 + curvature[1] * (2.0  * u1*v1) + curvature[2] * v1*v1;
	newcurvature_[1]  = curvature[0] * u1*u2 + curvature[1] * (u1*v2 + u2*v1) + curvature[2] * v1*v2;
	newcurvature_[2]  = curvature[0] * u2*u2 + curvature[1] * (2.0  * u2*v2) + curvature[2] * v2*v2;

    return newcurvature_;
}
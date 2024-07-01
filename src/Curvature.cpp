#include "Curvature.h"

Curvature::Curvature(CurvatureInput& input)
:pack_(input.pack){
    input.pack.ReadVectorString("outputs", ParameterPack::KeyType::Optional, OutputNames_);
    input.pack.ReadVectorString("outputFiles", ParameterPack::KeyType::Optional, OutputFileNames_);
    input.pack.ReadString("name", ParameterPack::KeyType::Optional, name_);

    outputs_.registerOutputFunc("curvature", [this](std::string name) -> void { this -> printCurvature(name);});
    outputs_.registerOutputFunc("principaldir", [this](std::string name) -> void {this -> printPrincipalDir(name);});
    outputs_.registerOutputFunc("FaceCurvature", [this](std::string name)-> void {this -> printFaceCurvature(name);});

    ASSERT(( OutputNames_.size() == OutputFileNames_.size()), "The number of outputs does not agree with number of output files.");
}

void Curvature::initialize(Mesh& mesh){
    const auto& vertices = mesh.getvertices();
    int nv = vertices.size();

    CurvaturePerVertex_.clear();
    CurvaturePerVertex_.resize(nv);

    avgCurvaturePerVertex_.clear();
    avgCurvaturePerVertex_.resize(nv ,0.0);

    GaussCurvaturePerVertex_.clear();
    GaussCurvaturePerVertex_.resize(nv,0.0);

    principalDir1_.clear();
    principalDir1_.resize(nv);

    principalDir2_.clear();
    principalDir2_.resize(nv);
}

void Curvature::printOutput(bool bootstrap, int numTimes){

    if (bootstrap){
        for (int i=0;i<OutputNames_.size();i++){
            std::string name = StringTools::AppendIndexToFileName(OutputFileNames_[i], numTimes);
            outputs_.getOutputFuncByName(OutputNames_[i])(name);
        }
    }
    else{
        for (int i=0;i<OutputNames_.size();i++){
            outputs_.getOutputFuncByName(OutputNames_[i])(OutputFileNames_[i]);
        }
    }
}

void Curvature::printCurvature(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# k1 k2 Average Gaussian" << "\n";
    for(int i=0;i<CurvaturePerVertex_.size();i++){
        ofs << CurvaturePerVertex_[i][0] << " " << CurvaturePerVertex_[i][1] << " " << avgCurvaturePerVertex_[i] << " " \
        << GaussCurvaturePerVertex_[i] << "\n";
    }
    ofs.close();
}

void Curvature::printFaceCurvature(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# TriangleIndex  MeanCurvature  GaussCurvature\n";
    for (int i=0;i<FaceCurvature_.size();i++)
    {
        ofs << i+1 << " " << FaceCurvature_[i] << "\n";
    }
    ofs.close();
}

void Curvature::printPrincipalDir(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    ofs << "#p1x p1y p1z p2x p2y p2z\n";

    for (int i=0;i<principalDir1_.size();i++){
        for (int j=0;j<3;j++){
            ofs << principalDir1_[i][j] << " ";
        }
        for (int j=0;j<3;j++){
            if (j<2){
                ofs << principalDir2_[i][j] << " ";
            }
            else{
                ofs << principalDir2_[i][j] << "\n";
            }
        }
    }

    ofs.close();
}

void Curvature::CalculateFaceCurvature(Mesh& mesh, const std::vector<Real>& VertexCurvature, \
                                                   const std::vector<Real>& VertexGaussCurvature,\
                                                   const std::vector<Real2>& CurvaturePerVertex, \
                                                   std::vector<std::array<Real,4>>& FaceCurvature)
{
    const auto& triangle = mesh.gettriangles();
    const auto& vertices = mesh.getvertices();

    FaceCurvature.clear();
    FaceCurvature.resize(triangle.size());

    #pragma omp parallel for
    for (int i=0;i<triangle.size();i++){
        auto& t = triangle[i].triangleindices_;
        Real FaceMeanC = 0.0;
        Real FaceGaussC= 0.0;
        Real FaceMeank1= 0.0;
        Real FaceMeank2= 0.0;

        for (int j=0;j<3;j++){
            FaceMeanC += VertexCurvature[t[j]];
            FaceGaussC+= VertexGaussCurvature[t[j]];
            FaceMeank1+= CurvaturePerVertex[t[j]][0];
            FaceMeank2+= CurvaturePerVertex[t[j]][1];
        }

        FaceMeanC = FaceMeanC / 3.0;
        FaceGaussC= FaceGaussC / 3.0;
        FaceMeank1= FaceMeank1 / 3.0;
        FaceMeank2= FaceMeank2 / 3.0;

        FaceCurvature[i][0] = FaceMeanC;
        FaceCurvature[i][1] = FaceGaussC;
        FaceCurvature[i][2] = FaceMeank1;
        FaceCurvature[i][3] = FaceMeank2;
    }
}
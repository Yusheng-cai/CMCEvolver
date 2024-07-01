#include "CurvatureQuadraticFit.hpp"

namespace CurvatureRegistry
{
    registry<CurvatureQuadraticFit> registerQuadraticfit("quadraticfit");
}

CurvatureQuadraticFit::CurvatureQuadraticFit(CurvatureInput& input)
: Curvature(input){
    pack_.Readbool("MonteCarlo", ParameterPack::KeyType::Optional, MonteCarlo_);
    if (MonteCarlo_){
        pack_.ReadNumber("MonteCarloN", ParameterPack::KeyType::Required, MonteCarloN_);
    }   
    pack_.ReadNumber("neighbors", ParameterPack::KeyType::Required, N_ring_);
}

void CurvatureQuadraticFit::calculate(Mesh& m){
    srand(time(NULL));

    // set the mesh pointer
    m_ = &m;

    std::vector<std::vector<int>> neighborIndices;
    std::vector<std::vector<int>> NearbyNeighbors;
    std::vector<Real> Areas;
    std::vector<Real3> normals;

    // find neighbors of each vertex
    MeshTools::CalculateVertexNeighbors(m, neighborIndices);

    // calculate vertex normals 
    m.CalcVertexNormals();

    // find k-ring neighbor using bfs
    Graph::BFS_kring_neighbor(neighborIndices, N_ring_, NearbyNeighbors);

    const auto& v = m.getvertices();
    CurvaturePerVertex_.resize(v.size());
    avgCurvaturePerVertex_.resize(v.size());
    GaussCurvaturePerVertex_.resize(v.size());

    for (int i=0;i<v.size();i++){
        ASSERT((NearbyNeighbors[i].size() >= 6), "Must have at least 6 points for curvature quadratic fit while neighbors of vertex " << i << " is " << neighborIndices[i].size());
        std::vector<int> neighbors = NearbyNeighbors[i];

        // make sure all the points have the same orientation as the middle one
        std::vector<int> new_neighbors;
        applyProjOnPlane(v, i, neighbors, new_neighbors);
        if (new_neighbors.size() >= 6 && new_neighbors.size() < neighbors.size()){
            neighbors = new_neighbors;
        }

        // projecting
        Real3 avg_normal = v[i].normals_;

        std::vector<int> sampled_points;
        MonteCarloSample(neighbors, sampled_points);
        neighbors = sampled_points;

        // compute reference frame --> based on the average normal in the frame
        std::vector<Real3> ref;
        computeReferenceFrame(i, neighborIndices[i], avg_normal, ref);

        // perform quadratic fit
        QuadraticFit(i, ref, neighbors);
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
    
    CalculateFaceCurvature(m, avgCurvaturePerVertex_, GaussCurvaturePerVertex_, CurvaturePerVertex_, FaceCurvature_);
}

void CurvatureQuadraticFit::QuadraticFit(int i, const std::vector<Real3>& ref, const std::vector<int>& neighbors){
    std::vector<Real3> shifted_points;
    const auto& v = m_->getvertices();

    for (auto r : ref){
        std::cout << r << "\n";
    }
    for (int j=0;j<neighbors.size();j++){
        Real3 cp = v[neighbors[j]].position_;

        Real3 vTang;
        Real temp;
        m_->getVertexDistance(cp, v[i].position_, vTang, temp);

        Real x,y,z;
        x = LinAlg3x3::DotProduct(vTang, ref[0]);
        y = LinAlg3x3::DotProduct(vTang, ref[1]);
        z = LinAlg3x3::DotProduct(vTang, ref[2]);
        Real3 shift = {x,y,z};

        shifted_points.push_back(shift);
    }

    Eigen::MatrixXd A(neighbors.size(),5);
    Eigen::MatrixXd b(neighbors.size(),1);
    Eigen::MatrixXd sol(5,1);

    for (int c=0;c<neighbors.size();c++){
        Real u = shifted_points[c][0];
        Real v = shifted_points[c][1];
        Real n = shifted_points[c][2];

        A(c,0) = u * u;
        A(c,1) = u * v;
        A(c,2) = v * v;
        A(c,3) = u;
        A(c,4) = v;

        b(c) = n;
    }

    sol = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    Real a,bb,c,d,e; 
    a = sol(0); 
    bb = sol(1); 
    c = sol(2); 
    d = sol(3);
    e = sol(4);

    Real E = 1.0  + d * d;
    Real F = d * e;
    Real G = 1.0 + e * e;

    Real3 n = {{-d,-e,1.0}};
    LinAlg3x3::normalize(n);

    Real L = 2.0 * a * n[2];
    Real M = bb * n[2];
    Real N = 2.0 * c * n[2];

    Eigen::Matrix2d m;
    m << L * G - M * F, M*E-L*F, M*E-L*F, N*E-M*F;
    m = m / (E*G - F*F);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(m);

    Eigen::Vector2d c_Val = eig.eigenvalues();
    c_Val = -c_Val;
    Real2 c_V = {{c_Val[0], c_Val[1]}};
    Algorithm::sort(c_V);

    // write to curvature per vertex 
    CurvaturePerVertex_[i] = c_V;
}

void CurvatureQuadraticFit::MonteCarloSample(const std::vector<int>& neighbors, std::vector<int>& sampled_points)
{
    sampled_points.clear();
    if (MonteCarloN_ > neighbors.size()){
        sampled_points = neighbors;

        return;
    }

    Real p =  (Real)MonteCarloN_ / (Real)neighbors.size();

    for (int i=0;i<neighbors.size();i++){
        Real r = (Real)rand() / (Real)RAND_MAX;

        if (r < p){
            sampled_points.push_back(neighbors[i]);
        }
    }
}

void CurvatureQuadraticFit::applyProjOnPlane(const std::vector<vertex>& v, int j, const std::vector<int>& neighbors, std::vector<int>& new_neighbors){
    Real3 ppn = v[j].normals_;
    new_neighbors.clear();

    for (int i=0;i<neighbors.size();i++){
        int neighbor_idx = neighbors[i];
        if (LinAlg3x3::DotProduct(ppn, v[neighbor_idx].normals_) < 0){
            new_neighbors.push_back(neighbor_idx);
        }
    }
}

void CurvatureQuadraticFit::computeReferenceFrame(int i, const std::vector<int>& neighbors, const Real3& normal, std::vector<Real3>& ref){
    ref.resize(3);

    // subtract v from first of neighbors 
    const auto& v = m_->getvertices();

    Real3 longest_v = v[neighbors[0]].position_;

    longest_v = Project(v[i].position_, longest_v, normal) - v[i].position_;
    LinAlg3x3::normalize(longest_v);

    // find the other 2 axis
    Real3 y_axis = LinAlg3x3::CrossProduct(normal, longest_v);
    LinAlg3x3::normalize(y_axis);
    ref[0] = longest_v;
    ref[1] = y_axis;
    ref[2] = normal;
}

CurvatureQuadraticFit::Real3 CurvatureQuadraticFit::Project(const Real3& v, const Real3& vp, const Real3& normal){
    Real3 diff;
    Real diffsq;
    m_-> getVertexDistance(vp, v, diff, diffsq);

    Real dot = LinAlg3x3::DotProduct(diff, normal);
    Real3 shift = normal * dot;

    return vp - shift;
}
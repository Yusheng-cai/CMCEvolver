#include "Curvature.h"
#include "Mesh.h"
#include "Graph.h"
#include "tools/CommonTypes.h"

class CurvatureQuadraticFit : public Curvature
{
    public:
        using INT2 = CommonTypes::index2;

        CurvatureQuadraticFit(CurvatureInput& input);

        virtual void calculate(Mesh& m) override;

        void MonteCarloSample(const std::vector<int>& neighbors, std::vector<int>& sampled_points);

        // check whether or not all the normals pointing in the same direction as the center one
        void applyProjOnPlane(const std::vector<vertex>& v, int j, const std::vector<int>& neighbors, std::vector<int>& new_neighbors);

        // projection
        Real3 Project(const Real3& v, const Real3& vp, const Real3& normals);

        // compute the reference frame
        void computeReferenceFrame(int i, const std::vector<int>& neighbors, const Real3& normal, std::vector<Real3>& ref);

        // quadratic fit
        void QuadraticFit(int i, const std::vector<Real3>& ref, const std::vector<int>& neighbors);

    private:
        int N_ring_, MonteCarloN_;
        bool MonteCarlo_=true;

        Mesh* m_=nullptr;

        int f_=0;
};
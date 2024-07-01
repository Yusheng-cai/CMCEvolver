#include "Curvature.h"

class CurvatureFDM: public Curvature
{
    public:
        CurvatureFDM(CurvatureInput& input);

        virtual void calculate(Mesh& mesh);
        virtual ~CurvatureFDM(){};

        virtual void printCurvature(std::string name) override;

    private:
        std::vector<std::vector<int>> neighbor_indices_;
};
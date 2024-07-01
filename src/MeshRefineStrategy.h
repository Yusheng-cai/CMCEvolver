#pragma once
#include "tools/Assert.h"
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"
#include "tools/GenericFactory.h"
#include "tools/OutputFunction.h"

#include <vector>
#include <array>
#include <map>


// forward declare Mesh
class Mesh;

struct MeshRefineStrategyInput{
    ParameterPack& pack;
};

class MeshRefineStrategy{
    public:
        using Real = CommonTypes::Real;
        using Real3 = CommonTypes::Real3;
        using INT2 = CommonTypes::index2;

        MeshRefineStrategy(MeshRefineStrategyInput& input);

        virtual ~MeshRefineStrategy(){};

        virtual void refine(Mesh& mesh) = 0;

        virtual void setFixed(const std::vector<int>& fixed_indices){fixed_indices_=fixed_indices;}

        virtual std::string getName() {return name_;}

    protected:
        Mesh* mesh_ = nullptr;

        ParameterPack& pack_;
        Output output_;
        std::string name_ = "refine";

        // neighbor Indices
        std::vector<std::vector<int>> neighborIndices_;
        std::map<INT2, std::vector<int>> MapEdgeToFace_;
        std::vector<std::vector<INT2>> MapVertexToEdge_;
        std::vector<bool> boundaryIndicator_;

        std::vector<int> fixed_indices_;
        std::vector<bool> isfixed_;
};

namespace MeshRefineStrategyFactory
{
    using key = std::string;
    using base = MeshRefineStrategy;

    using Factory = GenericFactory<base, key, MeshRefineStrategyInput&>;

    template<typename D>
    using registry_ = RegisterInFactory<base, D, key, MeshRefineStrategyInput&>;
};
#pragma once 

#include "tools/Assert.h"
#include "Curvature.h"
#include "MeshRefineStrategy.h"

#include <string>
#include <map>
#include <vector>

class Registry
{
    public:
        Registry()  =default;

        void registerCurvature(std::string name, const Curvature& curve);
        Curvature& getCurvature(std::string name);

        void registerMeshRefinement(std::string name, const MeshRefineStrategy& meshrefine);
        MeshRefineStrategy& getMeshRefineStrat(std::string name);

    private:
        std::map<std::string, Curvature*> MapNameToCurvature_;
        std::map<std::string, MeshRefineStrategy*> MapNameToMeshRefinement_;
};
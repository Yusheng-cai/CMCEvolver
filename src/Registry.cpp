#include "Registry.h"

void Registry::registerCurvature(std::string name, const Curvature& curve)
{
    auto it = MapNameToCurvature_.find(name);
    ASSERT((it == MapNameToCurvature_.end()), "The name for curvature " << name << " is already registered.");

    MapNameToCurvature_.insert(std::make_pair(name, const_cast<Curvature*>(&curve)));
}

Curvature& Registry::getCurvature(std::string name)
{
    auto it = MapNameToCurvature_.find(name);

    ASSERT((it != MapNameToCurvature_.end()), "The name for curvature " << name << " is not registered.");

    return *(it->second);
}

void Registry::registerMeshRefinement(std::string name, const MeshRefineStrategy& refineStrat)
{
    auto it = MapNameToMeshRefinement_.find(name);
    ASSERT((it == MapNameToMeshRefinement_.end()), "The name for mesh refinement " << name << " is already registered.");

    MapNameToMeshRefinement_.insert(std::make_pair(name, const_cast<MeshRefineStrategy*>(&refineStrat)));
}

MeshRefineStrategy& Registry::getMeshRefineStrat(std::string name)
{
    auto it = MapNameToMeshRefinement_.find(name);

    ASSERT((it != MapNameToMeshRefinement_.end()), "The name for mesh refinement " << name << " is not registered.");

    return *(it->second);
}
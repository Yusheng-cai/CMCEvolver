#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <list>

namespace Graph
{
    using index3 = CommonTypes::index3;
    using Real   = CommonTypes::Real;
    using Real3  = CommonTypes::Real3;

    // Need the input of nearby Indices of the vertex connected to it 
    void getNearbyIndicesNVertexAway(const std::vector<std::vector<int>>& NearbyIndices, int N, std::vector<std::vector<int>>& NerabyIndicesNVertex);

    void getNearbyIndicesNVertexAway(const std::vector<std::vector<int>>& NearbyIndices, int N, int index, std::vector<int>& NearbyIndicesNVertex);

    void getNNearbyIndices(const std::vector<std::vector<int>>& NearbyIndices, int N, std::vector<std::vector<int>>& NearbyIndicesN);

    // BFS k-ring neighbor
    void BFS_kring_neighbor(const std::vector<std::vector<int>>& NearbyIndices, int N, std::vector<std::vector<int>>& NNearbyNeighbors);

    // BFS k-ring neighbor for a single vertex 
    void BFS_kring_neighbor(const std::vector<std::vector<int>>& NearbyIndices, int N, int index, std::vector<int>& NNearbyNeighbors);

    // BFS k-ring neighbor but special for boundary vertices 
    void BFS_kring_neighbor_boundary(const std::vector<std::vector<int>>& NearbyIndices, const std::vector<bool>& boundaryIndicator, \
                                     int N, int Nboundary, std::vector<std::vector<int>>& NNearbyNeighbors);
};
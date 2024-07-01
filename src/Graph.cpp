#include "Graph.h"

void Graph::getNNearbyIndices(const std::vector<std::vector<int>>& NearbyIndices, int N, std::vector<std::vector<int>>& NearbyIndicesN)
{
    int NumberOfVertices = NearbyIndices.size();

    NearbyIndicesN.clear();
    NearbyIndicesN.resize(NumberOfVertices);

    #pragma omp parallel for
    for(int i=0;i<NumberOfVertices;i++)
    {
        int initial_count = 0;
        int final_count = 0;
        int nearbySize = NearbyIndices[i].size();
        const std::vector<int>& currentNearbyIndices = NearbyIndices[i];
        std::vector<int> currentNearbyIndicesNVertex;
        currentNearbyIndicesNVertex.resize(0);
        
        while (currentNearbyIndicesNVertex.size() < N){
            if (currentNearbyIndicesNVertex.size() == 0){
                for (int j=0;j<currentNearbyIndices.size();j++){
                    bool found = (std::find(currentNearbyIndicesNVertex.begin(), currentNearbyIndicesNVertex.end(), currentNearbyIndices[j])) != currentNearbyIndicesNVertex.end();
                    if ( ! found && currentNearbyIndices[j] != i){
                        currentNearbyIndicesNVertex.push_back(currentNearbyIndices[j]);
                    }
                }

                final_count = currentNearbyIndicesNVertex.size();
            }            
            else{
                std::vector<int> temp;
                temp.resize(0);
                for (int k=initial_count;k<final_count;k++){
                    int id = currentNearbyIndicesNVertex[k];

                    for (int m=0;m<NearbyIndices[id].size();m++){
                        bool found = (std::find(currentNearbyIndicesNVertex.begin(), currentNearbyIndicesNVertex.end(),NearbyIndices[id][m])) != currentNearbyIndicesNVertex.end();
                        if (! found && NearbyIndices[id][m] != i){
                            currentNearbyIndicesNVertex.push_back(NearbyIndices[id][m]);
                        }
                    }
                }

                initial_count = final_count;
                final_count = currentNearbyIndicesNVertex.size();
            }
        }

        int size = currentNearbyIndicesNVertex.size();
        std::vector<int> temp_;

        if (size > N){
            temp_.insert(temp_.end(),currentNearbyIndicesNVertex.begin(), currentNearbyIndicesNVertex.begin() + N);
            std::sort(temp_.begin(), temp_.end());
            NearbyIndicesN[i] = temp_;
        }
        else
        {
            std::sort(currentNearbyIndicesNVertex.begin(), currentNearbyIndicesNVertex.end());
            NearbyIndicesN[i] = currentNearbyIndicesNVertex;
        }
    }

}

void Graph::getNearbyIndicesNVertexAway(const std::vector<std::vector<int>>& NearbyIndices, int N, int index, std::vector<int>& NearbyIndicesNVertex){
    int initial_count = 0;
    int final_count = 0;
    int nearbySize = NearbyIndices[index].size();
    const std::vector<int>& currentNearbyIndices = NearbyIndices[index];
    NearbyIndicesNVertex.clear();

    for (int j=0;j<N;j++){
        if (NearbyIndicesNVertex.size() == 0){
            NearbyIndicesNVertex.insert(NearbyIndicesNVertex.end(),\
            currentNearbyIndices.begin(), currentNearbyIndices.end());

            final_count = NearbyIndicesNVertex.size();
        }            
        else{
            std::vector<int> temp;
            for (int k=initial_count;k<final_count;k++){
                int id = NearbyIndicesNVertex[k];

                for (int m=0;m<NearbyIndices[id].size();m++){
                    if (NearbyIndices[id][m] != index){
                        temp.push_back(NearbyIndices[id][m]);
                    }
                }
            }

            NearbyIndicesNVertex.insert(NearbyIndicesNVertex.end(), temp.begin(), temp.end());
            initial_count = final_count;
            final_count = NearbyIndicesNVertex.size();
        }
    }


    // sort the vector
    std::sort(NearbyIndicesNVertex.begin(), NearbyIndicesNVertex.end());

    // The non-unique terms will be undefined ***
    auto it = std::unique(NearbyIndicesNVertex.begin(), NearbyIndicesNVertex.end());

    // removes the non-unique terms
    NearbyIndicesNVertex.resize(std::distance(NearbyIndicesNVertex.begin(), it));
}

void Graph::getNearbyIndicesNVertexAway(const std::vector<std::vector<int>>& NearbyIndices, int N, std::vector<std::vector<int>>& NearbyIndicesNVertex)
{
    int NumberOfVertices = NearbyIndices.size();

    NearbyIndicesNVertex.clear();
    NearbyIndicesNVertex.resize(NumberOfVertices);

    #pragma omp parallel for
    for(int i=0;i<NumberOfVertices;i++){
        int initial_count = 0;
        int final_count = 0;
        int nearbySize = NearbyIndices[i].size();
        const std::vector<int>& currentNearbyIndices = NearbyIndices[i];
        std::vector<int> currentNearbyIndicesNVertex;
        currentNearbyIndicesNVertex.resize(0);

        for (int j=0;j<N;j++){
            if (currentNearbyIndicesNVertex.size() == 0){
                currentNearbyIndicesNVertex.insert(currentNearbyIndicesNVertex.end(),\
                currentNearbyIndices.begin(), currentNearbyIndices.end());

                final_count = currentNearbyIndicesNVertex.size();
            }            
            else{
                std::vector<int> temp;
                temp.resize(0);
                for (int k=initial_count;k<final_count;k++){
                    int id = currentNearbyIndicesNVertex[k];

                    for (int m=0;m<NearbyIndices[id].size();m++){
                        if (NearbyIndices[id][m] != i){
                            temp.push_back(NearbyIndices[id][m]);
                        }
                    }
                }

                currentNearbyIndicesNVertex.insert(currentNearbyIndicesNVertex.end(), temp.begin(), temp.end());
                initial_count = final_count;
                final_count = currentNearbyIndicesNVertex.size();
            }
        }


        // sort the vector
        std::sort(currentNearbyIndicesNVertex.begin(), currentNearbyIndicesNVertex.end());

        // The non-unique terms will be undefined ***
        auto it = std::unique(currentNearbyIndicesNVertex.begin(), currentNearbyIndicesNVertex.end());

        // removes the non-unique terms
        currentNearbyIndicesNVertex.resize(std::distance(currentNearbyIndicesNVertex.begin(), it));

        NearbyIndicesNVertex[i] = currentNearbyIndicesNVertex;
    }
}

void Graph::BFS_kring_neighbor(const std::vector<std::vector<int>>& NearbyIndices, int N, int i, std::vector<int>& NNearbyNeighbors){
    int numv = NearbyIndices.size();
    std::list<std::pair<int,int>> queue;
    std::vector<bool> visited(numv, false);
    queue.push_back(std::pair<int,int>(i,0));
    NNearbyNeighbors.clear();

    while (! queue.empty()){
        int toVisit = queue.front().first;
        int distance = queue.front().second;

        // removes the first element
        queue.pop_front();
        NNearbyNeighbors.push_back(toVisit);
        if (distance < N){
            for (int j=0;j<NearbyIndices[toVisit].size();j++){
                int neighbor = NearbyIndices[toVisit][j];
                if (! visited[neighbor]){
                    queue.push_back(std::pair<int,int>(neighbor, distance+1));
                    visited[neighbor] = true;
                }
            }
        }
    }
}

void Graph::BFS_kring_neighbor(const std::vector<std::vector<int>>& NearbyIndices, int N, std::vector<std::vector<int>>& NNearbyNeighbors)
{
    int numv = NearbyIndices.size();
    NNearbyNeighbors.clear(); NNearbyNeighbors.resize(numv);

    #pragma omp parallel for
    for (int i=0;i<numv;i++){
        std::vector<int> vv;
        BFS_kring_neighbor(NearbyIndices, N, i, vv);
        NNearbyNeighbors[i] = vv;
    }
}

void Graph::BFS_kring_neighbor_boundary(const std::vector<std::vector<int>>& NearbyIndices, const std::vector<bool>& boundaryIndicator, \
                                        int N, int Nboundary, std::vector<std::vector<int>>& NNearbyNeighbors)
{
    int numv = NearbyIndices.size();
    NNearbyNeighbors.clear(); NNearbyNeighbors.resize(numv);

    #pragma omp parallel for
    for (int i=0;i<numv;i++){
        std::vector<int> vv;
        if (boundaryIndicator[i]){
            BFS_kring_neighbor(NearbyIndices, Nboundary, i, vv);
        }
        else{
            BFS_kring_neighbor(NearbyIndices, N, i, vv);
        }

        NNearbyNeighbors[i] = vv;
    }
}
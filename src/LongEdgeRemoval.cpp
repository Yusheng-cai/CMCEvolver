/* This file is part of PyMesh. Copyright (c) 2015 by Qingnan Zhou */
#include "LongEdgeRemoval.h"

void LongEdgeRemoval::split(std::list<Real3>& end_points,
        std::list<Real3>::iterator begin_itr,
        Real length, Real threshold) {
    if (length <= threshold) return;
    if (!std::isfinite(length)) return;

    auto end_itr = std::next(begin_itr);
    ASSERT((end_itr != end_points.end()), "wrong");
    const auto& v1 = *begin_itr;
    const auto& v2 = *end_itr;
    Real3 distVec;
    Real dist;
    mesh_.getVertexDistance(v2, v1, distVec, dist);

    Real3 mid = v1 + 0.5 * distVec;
    Real half_length = length * 0.5;

    auto mid_itr = end_points.insert(end_itr, mid);

    split(end_points, begin_itr, half_length, threshold);
    split(end_points, mid_itr, half_length, threshold);
}

void LongEdgeRemoval::run(Real max_length, bool recursive) {
    int num_faces = m_faces_.size();
    m_ori_faces_ = Algorithm::arange(0, num_faces, 1);

    int num_added_triangles = 0;
    do {
        init_edge_map();
        split_long_edges(max_length);
        num_added_triangles = retriangulate();
        std::cout << "num_added_triangles = " << num_added_triangles << "\n";
        mesh_.SetVerticesAndTriangles(m_vertices_, m_faces_);
    } while (recursive && num_added_triangles != 0);

    mesh_.SetVerticesAndTriangles(m_vertices_, m_faces_);
}

void LongEdgeRemoval::init_edge_map() {
    m_edge_map_.clear();
    int num_faces = m_faces_.size();

    MeshTools::MapEdgeToFace(mesh_, m_edge_map_, false);
}

void LongEdgeRemoval::split_long_edges(Real max_length) {
    int num_ori_vertices = m_vertices_.size();
    int vertex_count = num_ori_vertices;
    std::vector<vertex> new_vertices;
    for (auto& itr : m_edge_map_) {
        const auto& edge = itr.first;
        std::list<Real3> chain;
        chain.push_back(m_vertices_[edge[0]].position_);
        chain.push_back(m_vertices_[edge[1]].position_);
        Real edge_length;
        Real3 edge_diff;
        mesh_.getVertexDistance(m_vertices_[edge[0]], m_vertices_[edge[1]], edge_diff, edge_length);
        split(chain, chain.begin(), edge_length, max_length);

        std::vector<int> vertex_indices;
        vertex_indices.push_back(edge[0]);
        const auto begin = std::next(chain.begin());
        const auto end = std::prev(chain.end());
        for (auto it = begin; it != end; it++) {
            vertex v;
            v.position_ = *it;
            new_vertices.push_back(v);
            vertex_indices.push_back(vertex_count);
            vertex_count++;
        }
        vertex_indices.push_back(edge[1]);

        itr.second = vertex_indices;
    }

    if (new_vertices.size() > 0) {
        m_vertices_.insert(m_vertices_.end(), new_vertices.begin(), new_vertices.end());
    }
}

int LongEdgeRemoval::retriangulate() {
    const int num_faces = m_faces_.size();
    ASSERT((num_faces == m_ori_faces_.size()), "wrong");

    std::vector<INT3> faces;
    std::vector<triangle> tris;
    std::vector<int> face_marks;
    for (int i=0; i<num_faces; i++) {
        int v0_idx, v1_idx, v2_idx;
        auto chain = get_vertex_chain_around_triangle(
                i, v0_idx, v1_idx, v2_idx);
        if (chain.size() == 3) {
            faces.push_back(m_faces_[i].triangleindices_);
        } else {
            triangulate_chain(faces, chain, v0_idx, v1_idx, v2_idx);
        }
        face_marks.push_back(faces.size());
    }

    tris.resize(faces.size());
    for (int i=0;i<faces.size();i++){
        triangle t;
        t.triangleindices_ = faces[i];
        tris[i] = t;
    }

    m_faces_ = tris;
    const int num_refined_faces = m_faces_.size();

    std::vector<int> ori_faces(num_refined_faces);
    int counter = 0;
    int fid = 0;
    for (const auto mark : face_marks) {
        for (; counter<mark; counter++) {
            ori_faces[counter] = m_ori_faces_[fid];
        }
        fid++;
    }
    m_ori_faces_.swap(ori_faces);

    return num_refined_faces - num_faces;
}

void LongEdgeRemoval::triangulate_chain(
        std::vector<INT3>& faces,
        const std::vector<int>& chain,
        int v0_idx, int v1_idx, int v2_idx) {
    const int chain_size = chain.size();
    auto next = [&](int i) { return (i+1) % chain_size; };
    auto prev = [&](int i) { return (i+chain_size-1) % chain_size; };
    auto length = [&](int vi, int vj) {
        Real dist;
        Real3 distVec;
        mesh_.getVertexDistance(m_vertices_[vi], m_vertices_[vj], distVec, dist);

        return dist;
    };

    std::vector<INT3> visited(chain_size, {0,0,0});
    visited[v0_idx][0] = 1;
    visited[v1_idx][1] = 1;
    visited[v2_idx][2] = 1;
    std::array<std::array<int,6>,3> candidates;
    candidates = {{{v0_idx, next(v0_idx), prev(v0_idx), 0, 0, 0},\
                  {v1_idx, next(v1_idx), prev(v1_idx), 0, 0, 0},\ 
                  {v2_idx, next(v2_idx), prev(v2_idx), 0, 0, 0}}};

    std::array<std::array<Real,2>,3> candidate_lengths;
    const Real NOT_USED = std::numeric_limits<Real>::max();
    candidate_lengths = {{{length(chain[candidates[0][1]], chain[candidates[0][2]]), NOT_USED}, \
                          {length(chain[candidates[1][1]], chain[candidates[1][2]]), NOT_USED}, \
                          {length(chain[candidates[2][1]], chain[candidates[2][2]]), NOT_USED}}};


    auto index_comp = [&](int i, int j) {
        // Use greater than operator so the queue is a min heap.
        Real min_i = Algorithm::min(candidate_lengths[i]);
        Real min_j = Algorithm::min(candidate_lengths[j]);

        return min_i > min_j;
    };
    std::priority_queue<int, std::vector<int>, decltype(index_comp)>
        Q(index_comp);
    Q.push(0);
    Q.push(1);
    Q.push(2);

    while (!Q.empty()) {
        int idx = Q.top();
        Q.pop();
        int selection;
        if (candidate_lengths[idx][0] != NOT_USED &&
                candidate_lengths[idx][0] <= candidate_lengths[idx][1]) {
            selection = 0;
        } else if (candidate_lengths[idx][1] != NOT_USED &&
                candidate_lengths[idx][1] < candidate_lengths[idx][0]){
            selection = 1;
        } else {
            continue;
        }
        int base_v = candidates[idx][selection * 3 + 0];
        int right_v = candidates[idx][selection * 3 + 1];
        int left_v = candidates[idx][selection * 3 + 2];
        ASSERT((visited[base_v][idx] >= 1), "wrong");
        if (Algorithm::sum(visited[base_v]) > 1 ||
                visited[right_v][idx] > 1 ||
                visited[left_v][idx] > 1) {
            candidate_lengths[idx][selection] = NOT_USED;
            Q.push(idx);
            continue;
        }

        visited[right_v][idx] = 1;
        visited[left_v][idx] = 1;
        visited[base_v][idx] = 2;
        faces.push_back({chain[base_v], chain[right_v], chain[left_v]});

        if (Algorithm::sum(visited[right_v]) == 1) {
            int right_to_right = next(right_v);
            Real edge_len = length(chain[left_v], chain[right_to_right]);
            candidate_lengths[idx][0] = edge_len;
            candidates[idx][0] = right_v;
            candidates[idx][1] = right_to_right;
            candidates[idx][2] = left_v;
        } else {
            candidate_lengths[idx][0] = NOT_USED;
        }
        if (Algorithm::sum(visited[left_v]) == 1) {
            int left_to_left = prev(left_v);
            Real edge_len = length(chain[right_v], chain[left_to_left]);
            candidate_lengths[idx][1] = edge_len;
            candidates[idx][3] = left_v;
            candidates[idx][4] = right_v;
            candidates[idx][5] = left_to_left;
        } else {
            candidate_lengths[idx][1] = NOT_USED;
        }
        Q.push(idx);
    }

    int num_count = 0;
    std::vector<int> visited_sum(chain_size, 0);
    for (int i=0;i<chain_size;i++){
        for (int j=0;j<3;j++){
            visited_sum[i] += visited[i][j];
        }

        if (visited_sum[i] > 1){
            num_count ++;
        }
    }
    if (num_count == 3) {
        INT3 face;
        int count = 0;
        for (int i=0; i<chain_size; i++) {
            if (visited_sum[i] > 1) {
                ASSERT((count < 3), "More than 3 vertices in a tri.");
                face[count] = chain[i];
                count++;
            }
        }
        faces.push_back(face);
    }
}

std::vector<int> LongEdgeRemoval::get_vertex_chain_around_triangle(
        int fi, int& v0_idx, int& v1_idx, int& v2_idx) {
    const auto& f = m_faces_[fi];
    INT2 edge1 = MeshTools::makeEdge(f[0], f[1]);
    INT2 edge2 = MeshTools::makeEdge(f[1], f[2]);
    INT2 edge3 = MeshTools::makeEdge(f[2], f[0]);
    std::vector<int> chain_01, chain_12, chain_20;
    bool found = Algorithm::FindInMap(m_edge_map_, edge1, chain_01);
    ASSERT((found), "edge " << edge1 << " not found.");
    found = Algorithm::FindInMap(m_edge_map_, edge2, chain_12);
    ASSERT((found), "edge " << edge2 << " not found.")
    found = Algorithm::FindInMap(m_edge_map_, edge3, chain_20);
    ASSERT((found), "edge " << edge3 << " not found.")

    std::vector<int> chain;
    chain.reserve(chain_01.size() + chain_12.size() + chain_20.size() - 3);
    v0_idx = 0;
    chain.push_back(f[0]);
    if (f[0] == chain_01.front()) {
        chain.insert(chain.end(),
                std::next(chain_01.begin()),
                std::prev(chain_01.end()));
    } else {
        ASSERT((f[0] == chain_01.back()), "wrong");
        chain.insert(chain.end(),
                std::next(chain_01.rbegin()),
                std::prev(chain_01.rend()));
    }

    v1_idx = chain.size();
    chain.push_back(f[1]);
    if (f[1] == chain_12.front()) {
        chain.insert(chain.end(),
                std::next(chain_12.begin()),
                std::prev(chain_12.end()));
    } else {
        ASSERT((f[1] == chain_12.back()), "wrong");
        chain.insert(chain.end(),
                std::next(chain_12.rbegin()),
                std::prev(chain_12.rend()));
    }

    v2_idx = chain.size();
    chain.push_back(f[2]);
    if (f[2] == chain_20.front()) {
        chain.insert(chain.end(),
                std::next(chain_20.begin()),
                std::prev(chain_20.end()));
    } else {
        ASSERT((f[2] == chain_20.back()), "wrong");
        chain.insert(chain.end(),
                std::next(chain_20.rbegin()),
                std::prev(chain_20.rend()));
    }
    return chain;
}

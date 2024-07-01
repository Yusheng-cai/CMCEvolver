/* This file is part of PyMesh. Copyright (c) 2015 by Qingnan Zhou */
#include "ShortEdgeRemoval.h"


const int ShortEdgeRemoval::UNMAPPED = std::numeric_limits<int>::max();
const Real  ShortEdgeRemoval::INFINITE = std::numeric_limits<Real>::max();

ShortEdgeRemoval::ShortEdgeRemoval(Mesh& mesh) :
    m_heap_(false),
    m_vertices_(mesh.getvertices()),
    m_faces_(mesh.gettriangles()),
    m_num_collapsed_(0), 
    mesh_(mesh)
{
    m_importance_ = std::vector<int>(m_vertices_.size(), 0);
}

int ShortEdgeRemoval::calculate(Real threshold) {
    int num_collapsed = m_num_collapsed_;
    init();
    do {
        num_collapsed = m_num_collapsed_;
        collapse(threshold);
        update();
        if (num_collapsed == m_num_collapsed_) break;
    } while (get_num_faces() > 0 && min_edge_length() <= threshold);


    return m_num_collapsed_;
}

std::vector<vertex> ShortEdgeRemoval::get_vertices() const {
    const int num_vertices = get_num_vertices();
    const int num_ori_vertices = m_vertices_.size();
    const int num_new_vertices = m_new_vertices_.size();

    std::vector<vertex> vertices;
    vertices.insert(vertices.end(), m_vertices_.begin(), m_vertices_.end());
    vertices.insert(vertices.end(), m_new_vertices_.begin(), m_new_vertices_.end());

    return vertices;
}

void ShortEdgeRemoval::init() {
    init_vertex_map();
    init_edges();
    init_edge_length_heap();
    init_face_indices();
    init_vertex_face_neighbors();
}

void ShortEdgeRemoval::update() {
    update_faces();
    update_importance();
    update_vertices();
    mesh_.SetVerticesAndTriangles(m_vertices_, m_faces_);
    init_vertex_map();
    init_edges();
    init_edge_length_heap();
    init_vertex_face_neighbors();
}

void ShortEdgeRemoval::init_vertex_map() {
    const int num_vertices = get_num_vertices();
    m_vertex_map_.clear();
    m_vertex_map_.resize(num_vertices);
    std::fill(m_vertex_map_.begin(), m_vertex_map_.end(), UNMAPPED);
}

void ShortEdgeRemoval::init_face_indices() {
    const int num_faces = get_num_faces();
    m_face_indices_.resize(num_faces);

    for (int i=0; i<num_faces; i++) {
        m_face_indices_[i] = i;
    }
}

void ShortEdgeRemoval::init_vertex_face_neighbors() {
    m_vertex_face_neighbors_.clear();
    MeshTools::MapVerticesToFaces(mesh_, m_vertex_face_neighbors_);
}

void ShortEdgeRemoval::init_edges() {
    std::map<INT2, std::vector<int>> edges;
    std::vector<std::vector<INT2>> temp;
    MeshTools::MapEdgeToFace(mesh_, edges, temp, false);

    m_edges_.clear();
    for (auto it = edges.begin(); it != edges.end(); it++){
        m_edges_.push_back(it -> first);
    }
}

void ShortEdgeRemoval::init_edge_length_heap() {
    const int num_edges = m_edges_.size();

    std::vector<Real> edge_lengths(num_edges);
    for (int i=0; i<num_edges; i++) {
        if (edge_can_be_collapsed(i)) {
            edge_lengths[i] = compute_edge_length(m_edges_[i]);
        } else {
            edge_lengths[i] = INFINITE;
        }
    }
    m_heap_.init(edge_lengths);
}

void ShortEdgeRemoval::update_vertices() {
    m_vertices_ = get_vertices();
    m_new_vertices_.clear();
}

void ShortEdgeRemoval::update_faces() {
    const int num_faces = m_faces_.size();

    std::vector<triangle> faces;
    std::vector<int> face_indices;
    for (int i=0; i<num_faces; i++) {
        for (int j=0; j<3; j++) {
            int mapped_idx = m_vertex_map_[m_faces_[i][j]];
            if (mapped_idx != UNMAPPED)
                m_faces_[i][j] = mapped_idx;
        }

        // if we have degenerate triangles 
        if (m_faces_[i][0] == m_faces_[i][1] ||
            m_faces_[i][1] == m_faces_[i][2] ||
            m_faces_[i][2] == m_faces_[i][0])
            continue;
        faces.push_back(m_faces_[i]);
        face_indices.push_back(m_face_indices_[i]);
    }

    m_faces_.clear();
    m_faces_.insert(m_faces_.end(), faces.begin(), faces.end());

    m_face_indices_.clear();
    m_face_indices_.insert(m_face_indices_.end(), face_indices.begin(), face_indices.end());
}

void ShortEdgeRemoval::update_importance() {
    const int num_ori_vertices = m_vertex_map_.size();
    const int num_vertices = get_num_vertices();
    std::vector<int> extra(num_vertices-num_ori_vertices, 0);
    m_importance_.insert(m_importance_.end(), extra.begin(), extra.end());

    for (int i=0; i<num_ori_vertices; i++) {
        int mapped_idx = m_vertex_map_[i];
        if (mapped_idx != UNMAPPED) {
            int cur_val = m_importance_[mapped_idx];
            int old_val = m_importance_[i];
            assert(cur_val >= -1);
            assert(old_val >= -1);
            if (cur_val < 0 || old_val < 0) {
                m_importance_[mapped_idx] = -1;
            } else {
                m_importance_[mapped_idx] = std::max(cur_val, old_val);
            }
        }
    }
}

void ShortEdgeRemoval::collapse(Real threshold) {
    while (!m_heap_.empty()) {
        int edge_idx = m_heap_.top();
        Real edge_len = m_heap_.top_value();

        m_heap_.pop();

        if (edge_len > threshold) break;
        if (!edge_is_valid(edge_idx)) continue;
        if (!edge_can_be_collapsed(edge_idx)) continue;

        collapse_edge(edge_idx);
    }
}

bool ShortEdgeRemoval::edge_is_valid(int edge_idx) const {
    const INT2& edge = m_edges_[edge_idx];
    int v1_idx = edge[0];
    int v2_idx = edge[1];
    return m_vertex_map_[v1_idx] == UNMAPPED &&
           m_vertex_map_[v2_idx] == UNMAPPED;
}

bool ShortEdgeRemoval::edge_can_be_collapsed(int edge_idx) const {
    const INT2& edge = m_edges_[edge_idx];
    int v1_idx = edge[0];
    int v2_idx = edge[1];
    return (m_importance_[v1_idx] >= 0) || (m_importance_[v2_idx] >= 0);
}

bool ShortEdgeRemoval::collapse_would_cause_fold_over(
        int edge_idx, const ShortEdgeRemoval::Real3& v) const {
    const INT2& e = m_edges_[edge_idx];
    const int i1 = e[0];
    const int i2 = e[1];
    const auto& v1_face_neighbors = m_vertex_face_neighbors_[i1];
    const auto& v2_face_neighbors = m_vertex_face_neighbors_[i2];
    if (faces_would_flip(i1, i2, v, v1_face_neighbors)) return true;
    if (faces_would_flip(i1, i2, v, v2_face_neighbors)) return true;
    return false;
}

bool ShortEdgeRemoval::faces_would_flip(int i1, int i2,
        const Real3& v,
        const std::vector<int>& faces) const {
    auto index_of = [=](const INT3& array, int val) -> int {
        if (array[0] == val) return 0;
        if (array[1] == val) return 1;
        if (array[2] == val) return 2;
        return std::numeric_limits<int>::max();
    };
    for (const auto fi : faces) {
        const INT3& f = m_faces_[fi].triangleindices_;
        int local_i1 = index_of(f, i1);
        int local_i2 = index_of(f, i2);
        if (local_i1 < 3 && local_i2 < 3) {
            // This face contains the edge and would be eliminated after
            // collapse.
            continue;
        } else if (local_i1 < 3) {
            Real3 v_old = get_vertex(i1);
            const Real3 vo1 = get_vertex(f[(local_i1+1)%3]);
            const Real3 vo2 = get_vertex(f[(local_i1+2)%3]);
            if (face_would_flip(v_old, v, vo1, vo2)) { return true; }
        } else if (local_i2 < 3) {
            Real3 v_old = get_vertex(i2);
            const Real3 vo1 = get_vertex(f[(local_i2+1)%3]);
            const Real3 vo2 = get_vertex(f[(local_i2+2)%3]);
            if (face_would_flip(v_old, v, vo1, vo2)) { return true; }
        }
    }
    return false;
}

bool ShortEdgeRemoval::face_would_flip(
        const Real3& v_old, const Real3& v_new,
        const Real3& v_o1,  const Real3& v_o2) const {
    Real3 en1, en2, eo1, eo2;
    Real d;
    mesh_.getVertexDistance(v_o1, v_new, en1, d);
    mesh_.getVertexDistance(v_o2, v_new, en2, d);
    mesh_.getVertexDistance(v_o1, v_old, eo1, d);
    mesh_.getVertexDistance(v_o2, v_old, eo2, d);

    Real3 normal_new = LinAlg3x3::CrossProduct(en1, en2);
    Real norm_new = LinAlg3x3::norm(normal_new);
    Real3 normal_old = LinAlg3x3::CrossProduct(eo1, eo2);
    Real norm_old = LinAlg3x3::norm(normal_old);

    if (!norm_old < 1e-8) {
        // If old triangle is degenerated, things won't get worse if we
        // collapse.
        return false;
    }
    if (!norm_new < 1e-8) {
        // To not collapse if new Triangle is degenerated while old one is not.
        return true;
    }

    normal_new = normal_new / norm_new;
    normal_old = normal_old / norm_old;
    Real dot = LinAlg3x3::DotProduct(normal_new, normal_old);

    return dot < 0.5;
}

void ShortEdgeRemoval::collapse_edge(int edge_idx) {
    // get the edge
    const INT2& e = m_edges_[edge_idx];
    const int i1 = e[0];
    const int i2 = e[1];

    // find original vertices and new vertices size
    const int num_ori_vertices = m_vertices_.size();
    const int num_new_vertices = m_new_vertices_.size();

    // Obtain the vertex position and its importance 
    const Real3 v1 = get_vertex(i1);
    const Real3 v2 = get_vertex(i2);
    const int v1_importance = m_importance_[i1];
    const int v2_importance = m_importance_[i2];

    Real3 new_v;
    if (v1_importance < 0) {
        ASSERT((v2_importance >= 0), "v1 is negative so v2 importance has to be positive while it is " << v2_importance);
        new_v = v1;
    } else if (v2_importance < 0) {
        ASSERT((v1_importance>=0), "v2 wrong");
        new_v = v2;
    } else {
        if (v1_importance == v2_importance) {
            Real3 distVec;
            Real dist;
            // v2 - v1 + L
            mesh_.getVertexDistance(v2,v1, distVec, dist);
            new_v = v1 + 0.5 * distVec;
        } else if (v1_importance > v2_importance) {
            new_v = v1;
        } else {
            new_v = v2;
        }
    }

    if (collapse_would_cause_fold_over(edge_idx, new_v)) { return; }

    vertex v;
    v.position_ = new_v; 
    m_new_vertices_.push_back(v);
    int idx_mid = num_ori_vertices + num_new_vertices;
    m_vertex_map_[i1] = idx_mid;
    m_vertex_map_[i2] = idx_mid;

    m_num_collapsed_++;
}

ShortEdgeRemoval::Real3 ShortEdgeRemoval::get_vertex(int i) const {
    const int num_ori_vertices = m_vertices_.size();
    if (i<num_ori_vertices) {
        return m_vertices_[i].position_;
    } else {
        i -= num_ori_vertices;
        return m_new_vertices_[i].position_;
    }
}

Real ShortEdgeRemoval::min_edge_length() const {
    ASSERT((m_heap_.size() > 0), "Edge heap is empty!");

    return m_heap_.top_value();
}

ShortEdgeRemoval::Real ShortEdgeRemoval::compute_edge_length(const INT2& e) const {
    int i1 = e[0];
    int i2 = e[1];

    Real3 v1 = get_vertex(i1);
    Real3 v2 = get_vertex(i2);
    Real dist;
    Real3 distVec;
    mesh_.getVertexDistance(v1, v2, distVec, dist);
    return dist;
}

int ShortEdgeRemoval::get_num_vertices() const {
    const int num_ori_vertices = m_vertices_.size();
    const int num_new_vertices = m_new_vertices_.size();
    const int num_vertices = num_ori_vertices + num_new_vertices;
    return num_vertices;
}

int ShortEdgeRemoval::get_num_faces() const {
    return m_faces_.size();
}
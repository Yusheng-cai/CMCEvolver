/* This file is part of PyMesh. Copyright (c) 2015 by Qingnan Zhou */
/* Modified by Yusheng Cai to be part of the InstantaneousInterface package.*/

#include "ObtuseTriangleRemoval.h"

ObtuseTriangleRemoval::ObtuseTriangleRemoval(
        Mesh& mesh)
    : m_vertices_(mesh.getvertices()), m_faces_(mesh.gettriangles()), mesh_(mesh) { }

int ObtuseTriangleRemoval::calculate(Real max_angle_allowed, int max_iterations) {
    int total_num_split = 0;

    // convert angle to radian
    max_angle_allowed = max_angle_allowed * Constants::PI / 180;

    // start doing the iterations
    int count = 0;
    do {
        clear_intermediate_data();
        set_all_faces_as_valid();
        compute_face_angles();
        compute_opposite_vertices();
        compute_edge_face_adjacency();

        int num_split = split_obtuse_triangles(max_angle_allowed);
        total_num_split += num_split;
        count++;

        finalize_geometry();
        if (num_split == 0) break;
    } while (count < max_iterations);
    mesh_.SetVerticesAndTriangles(m_vertices_, m_faces_);

    return total_num_split;
}

void ObtuseTriangleRemoval::clear_intermediate_data() {
    m_face_angles_.clear();
    m_opp_vertices_.clear();
    m_edges_.clear();
    m_valid_.clear();
    m_edge_faces_.clear();
    m_new_vertices_.clear();
    m_new_faces_.clear();
}

void ObtuseTriangleRemoval::set_all_faces_as_valid() {
    m_valid_ = std::vector<bool>(m_faces_.size(), true);
}

void ObtuseTriangleRemoval::compute_face_angles() {
    const int num_faces = m_faces_.size();
    m_face_angles_.resize(num_faces*3);
    for (int i=0; i<num_faces; i++) {
        Real3 face_angles;
        Real3 v1 = m_vertices_[m_faces_[i][0]].position_;
        Real3 v2 = m_vertices_[m_faces_[i][1]].position_;
        Real3 v3 = m_vertices_[m_faces_[i][2]].position_;

        Real3 v2v1, v3v1, v3v2, v1v2, v1v3, v2v3;
        Real d;
        mesh_.getVertexDistance(v2, v1, v2v1, d);
        mesh_.getVertexDistance(v3, v1, v3v1, d);
        mesh_.getVertexDistance(v3, v2, v3v2, d);
        mesh_.getVertexDistance(v1, v2, v1v2, d);
        mesh_.getVertexDistance(v1, v3, v1v3, d);
        mesh_.getVertexDistance(v2, v3, v2v3, d);

        m_face_angles_[i*3]   = LinAlg3x3::findAngle(v2v1, v3v1);
        m_face_angles_[i*3+1] = LinAlg3x3::findAngle(v3v2, v1v2);
        m_face_angles_[i*3+2] = LinAlg3x3::findAngle(v1v3, v2v3);
    }
}

void ObtuseTriangleRemoval::compute_opposite_vertices() {
    const int num_faces = m_faces_.size();
    m_opp_vertices_.reserve(num_faces * 3);
    for (int i=0; i<num_faces; i++) {
        m_opp_vertices_.push_back(0);
        m_opp_vertices_.push_back(1);
        m_opp_vertices_.push_back(2);
    }
}

void ObtuseTriangleRemoval::compute_edge_face_adjacency() {
    int num_faces = m_faces_.size();
    m_edges_.clear();
    for (int i=0; i<num_faces; i++) {
        INT2 e1 = MeshTools::makeEdge(m_faces_[i][1], m_faces_[i][2]);
        INT2 e2 = MeshTools::makeEdge(m_faces_[i][0], m_faces_[i][2]);
        INT2 e3 = MeshTools::makeEdge(m_faces_[i][0], m_faces_[i][1]);

        m_edges_.push_back(e1);
        m_edges_.push_back(e2);
        m_edges_.push_back(e3);
    }

    for (int i=0; i<3*num_faces; i++) {
        INT2 edge = m_edges_[i];
        m_edge_faces_.insert(std::make_pair(edge, i));
    }
}

bool ObtuseTriangleRemoval::edge_can_be_splited(int ext_idx) const {
    bool result = true;
    const INT2& e = m_edges_[ext_idx];
    auto adj_faces = m_edge_faces_.equal_range(e);
    for (auto it = adj_faces.first;
            it != adj_faces.second; it++) {
        int f_ext_idx = it->second;
        int f_idx = f_ext_idx / 3;
        result &= m_valid_[f_idx];
    }
    
    return result;
}

int ObtuseTriangleRemoval::split_obtuse_triangles(Real max_angle) {
    int num_splited = 0;
    IndexHeap<Real> candidates(m_face_angles_);
    while (!candidates.empty()) {
        int ext_idx = candidates.top();
        candidates.pop();
        if (m_face_angles_[ext_idx] < max_angle) break;
        if (!edge_can_be_splited(ext_idx)) continue;

        split_triangle(ext_idx);
        num_splited ++;
    }
    return num_splited;
}

void ObtuseTriangleRemoval::split_triangle(int ext_idx) {
    INT2 edge = m_edges_[ext_idx];

    Real3 proj_point = project(ext_idx);
    vertex newv;
    newv.position_ = proj_point;
    m_new_vertices_.push_back(newv);
    const int num_old_v = m_vertices_.size();
    int v_mid = num_old_v + m_new_vertices_.size() - 1;

    auto range = m_edge_faces_.equal_range(edge);

    for (auto itr = range.first; itr != range.second; itr ++) {
        int f_ext_idx = itr -> second;
        int f_idx = f_ext_idx / 3;
        int opp_idx = m_opp_vertices_[f_ext_idx];
        int v0 = m_faces_[f_idx][opp_idx];
        int v1 = m_faces_[f_idx][(opp_idx+1)%3];
        int v2 = m_faces_[f_idx][(opp_idx+2)%3];

        triangle f1, f2;
        f1.triangleindices_ = {{v0, v1, v_mid}};
        f2.triangleindices_ = {{v0, v_mid, v2}};

        m_new_faces_.push_back(f1);
        m_new_faces_.push_back(f2);

        m_valid_[f_idx] = false;
    }
}

ObtuseTriangleRemoval::Real3 ObtuseTriangleRemoval::project(int ext_idx) {
    int f_idx = ext_idx / 3;
    int opp_idx = m_opp_vertices_[ext_idx];
    INT2& e = m_edges_[ext_idx];
    Real3 v1 = m_vertices_[e[0]].position_;
    Real3 v2 = m_vertices_[e[1]].position_;
    Real3 v3 = m_vertices_[m_faces_[f_idx][opp_idx]].position_;

    Real3 e1, e2;
    Real e1sq, e2sq;
    mesh_.getVertexDistance(v2, v1, e1, e1sq);
    mesh_.getVertexDistance(v3, v1, e2, e2sq);

    Real dot = LinAlg3x3::DotProduct(e1, e2);
    Real3 offset = dot / (e1sq*e1sq) * e1;
    Real3 proj = v1 + offset;
    return proj;
}

void ObtuseTriangleRemoval::finalize_geometry() {
    finalize_vertices();
    finalize_faces();
}

void ObtuseTriangleRemoval::finalize_vertices() {
    if (m_new_vertices_.size() == 0) return;
    m_vertices_.insert(m_vertices_.end(), m_new_vertices_.begin(), m_new_vertices_.end());
}

void ObtuseTriangleRemoval::finalize_faces() {
    if (m_new_faces_.size()==0) return;

    std::vector<triangle> valid_faces;
    int num_ori_f = m_faces_.size();
    for (int i=0; i<num_ori_f; i++) {
        if (m_valid_[i]) {
            valid_faces.push_back(m_faces_[i]);
        }
    }
    valid_faces.insert(valid_faces.end(),
            m_new_faces_.begin(), m_new_faces_.end());
    m_faces_.clear();
    m_faces_.insert(m_faces_.end(), valid_faces.begin(), valid_faces.end());
}

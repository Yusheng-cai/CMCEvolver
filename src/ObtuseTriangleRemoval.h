/* This file is part of PyMesh. Copyright (c) 2015 by Qingnan Zhou */
/* Modified by Yusheng Cai to be part of the InstantaneousInterface package.*/
#pragma once
#include <vector>
#include <map>
#include <iostream>

#include "tools/CommonTypes.h"
#include "tools/Algorithm.h"
#include "Mesh.h"
#include "IndexHeap.h"

class ObtuseTriangleRemoval {
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using INT3 = CommonTypes::index3;
        using INT2 = CommonTypes::index2;

        ObtuseTriangleRemoval(Mesh& mesh);

    public:
        // Angle in radian
        int calculate(Real max_angle_allowed, int max_iterations=1);
        std::vector<vertex> get_vertices() const { return m_vertices_; }
        std::vector<triangle> get_faces() const { return m_faces_; }

    private:
        void clear_intermediate_data();
        void set_all_faces_as_valid();
        void compute_face_angles();
        void compute_opposite_vertices();
        void compute_edge_face_adjacency();
        bool edge_can_be_splited(int ext_idx) const;

        int split_obtuse_triangles(Real max_angle);
        void split_triangle(int ext_idx);
        Real3 project(int ext_idx);

        void finalize_geometry();
        void finalize_vertices();
        void finalize_faces();

    private:
        std::vector<Real> m_face_angles_;
        std::vector<int> m_opp_vertices_;
        std::vector<INT2> m_edges_;
        std::vector<bool> m_valid_;
        std::multimap<INT2, int> m_edge_faces_;
        std::vector<vertex> m_new_vertices_;
        std::vector<triangle> m_new_faces_;
        Mesh& mesh_;

        std::vector<vertex> m_vertices_;
        std::vector<triangle> m_faces_;
};
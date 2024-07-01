/* This file is part of PyMesh. Copyright (c) 2015 by Qingnan Zhou */
#pragma once

#include <vector>
#include <functional>
#include <iostream>
#include <list>
#include <queue>
#include <unordered_map>
#include <map>

#include "tools/Algorithm.h"
#include "Mesh.h"
#include "tools/CommonTypes.h"

class LongEdgeRemoval {
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using INT2 = CommonTypes::index2;
        using INT3 = CommonTypes::index3;

        LongEdgeRemoval(Mesh& mesh) :
            m_vertices_(mesh.getvertices()), m_faces_(mesh.gettriangles()), mesh_(mesh) {}

    public:
        void run(Real max_length, bool recursive=true);

        void split(std::list<Real3>& end_points, \
                    std::list<Real3>::iterator begin_itr, \
                    Real length, Real threshold);

        std::vector<vertex> get_vertices() const { return m_vertices_; }
        std::vector<triangle> get_faces() const { return m_faces_; }

        std::vector<int> get_ori_faces() const { return m_ori_faces_; }

    private:
        void init_edge_map();
        void split_long_edges(Real max_length);
        int retriangulate();
        void triangulate_chain(
                std::vector<INT3>& faces,
                const std::vector<int>& chain,
                int v0_idx, int v1_idx, int v2_idx);
        std::vector<int> get_vertex_chain_around_triangle(
                int fi, int& v0_idx, int& v1_idx, int& v2_idx);

    private:
        Mesh& mesh_;

        std::vector<vertex> m_vertices_;
        std::vector<triangle> m_faces_;
        std::vector<int>  m_ori_faces_;

        std::map<INT2,std::vector<int>> m_edge_map_;
};
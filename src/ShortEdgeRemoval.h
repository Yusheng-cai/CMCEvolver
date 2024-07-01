/* This file is part of PyMesh. Copyright (c) 2015 by Qingnan Zhou */
#pragma once
#include <vector>
#include <numeric>
#include <unordered_map>
#include <map>

#include "Mesh.h"
#include "tools/CommonTypes.h"
#include "IndexHeap.h"

class ShortEdgeRemoval {
    public:
        using Real = CommonTypes::Real;
        using Real3 = CommonTypes::Real3;
        using INT3 = CommonTypes::index3;
        using INT2 = CommonTypes::index2;

        ShortEdgeRemoval(Mesh& mesh);

    public:
        /**
         * Importance is an integer value per vertex.  The vertex with higher
         * importance would keep its position during edge collapsing.  Mid-point
         * is used when collapsing edge with same importance at each end.
         *
         * Vertex with negative importance is not considered during collapsing.
         * i.e. they will keep their original location.
         */
        void set_importance(const std::vector<int>& importance) {
            m_importance_ = importance;
        }
        /**
         * Remove all edges that <= thresold
         * If thresold=0, remove all degenerated edges.
         */
        int calculate(Real threshold);
        std::vector<vertex> get_vertices() const;
        std::vector<triangle> get_faces() const { return m_faces_; }
        std::vector<int>  get_face_indices() const { return m_face_indices_; }

    private:
        void init();
        void update();
        void init_vertex_map();
        void init_face_indices();
        void init_edges();
        void init_edge_length_heap();
        void init_vertex_face_neighbors();
        void update_vertices();
        void update_faces();
        void update_importance();
        void collapse(Real threshold);
        bool edge_is_valid(int edge_idx) const;
        bool edge_can_be_collapsed(int edge_idx) const;
        bool collapse_would_cause_fold_over(int edge_idx,
                const Real3& v) const;
        bool faces_would_flip(int i1, int i2, const Real3& v,
                const std::vector<int>& faces) const;
        bool face_would_flip(const Real3& v_old, const Real3& v_new,
                const Real3& v_o1, const Real3& v_o2) const;
        void collapse_edge(int edge_idx);
        Real3 get_vertex(int i) const;
        Real min_edge_length() const;
        Real compute_edge_length(const INT2& e) const;
        int get_num_vertices() const;
        int get_num_faces() const;

    private:
        std::vector<int> m_vertex_map_;
        std::vector<INT2> m_edges_;
        IndexHeap<Real> m_heap_;

        Mesh& mesh_;

        std::vector<std::vector<int>> m_vertex_face_neighbors_;

        std::vector<vertex> m_vertices_;
        std::vector<triangle> m_faces_;
        std::vector<int>  m_face_indices_;
        std::vector<int>  m_importance_;

        std::vector<vertex> m_new_vertices_;

        int m_num_collapsed_;

    private:
        static const int UNMAPPED;
        static const Real INFINITE;
};
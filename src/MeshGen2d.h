#pragma once

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Constrained_Delaunay_triangulation_2.h"
#include "CGAL/Delaunay_mesher_2.h"
#include "CGAL/Delaunay_mesh_face_base_2.h"
#include "CGAL/Delaunay_mesh_size_criteria_2.h"
#include "CGAL/Triangulation_conformer_2.h"
#include "CGAL/IO/write_VTU.h"
#include "CGAL/lloyd_optimize_mesh_2.h"
#include "Mesh.h"

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <list>
#include <algorithm>
#include <sstream>
#include <map>
#include <limits>

class MeshGen2d
{
    public:
        using K        = CGAL::Exact_predicates_inexact_constructions_kernel;
        using Vb       = CGAL::Triangulation_vertex_base_2<K>;
        using Fb       = CGAL::Delaunay_mesh_face_base_2<K>;
        using Tds      = CGAL::Triangulation_data_structure_2<Vb, Fb>;
        using CDT      = CGAL::Constrained_Delaunay_triangulation_2<K, Tds>;
        using Criteria = CGAL::Delaunay_mesh_size_criteria_2<CDT>;
        using Mesher   = CGAL::Delaunay_mesher_2<CDT,Criteria>;
        using Vertex_handle = CDT::Vertex_handle;
        using Point    = CDT::Point;
        using Real2    = CommonTypes::Real2;
        using Real3    = CommonTypes::Real3;
        using INT2     = CommonTypes::index2;
        using Real     = CommonTypes::Real;
        using INT3     = CommonTypes::index3;
        
        MeshGen2d(std::string file, Real aspect_bound=0.125, Real size_bound=1.0);
        MeshGen2d(std::vector<Real2>& vertices, std::vector<INT2>& edges, std::vector<Real2>& seed_point, Real aspect_bound=0.125, Real size_bound=1.0,  bool periodic=true);
        void setBoxLength(Real2& box){box_length_=box; isPBC_=true;}
        void generate();

        // The input file I designed has the following features 
        // File 
        // # of points
        // # of constraints
        // # of seeds 
        // Pointsx , Pointsy
        // idx1, idx2
        void InputFileReader();

        // check the max distances between the edges 
        void checkMaxDistanceBetweenEdges();

        bool isPBC() {return isPBC_;}

        // find distance between 2 points taking into account the pbc based on what the class is provided with
        void findDistance(const Real2& p1, const Real2& p2, Real2& dist, Real& distsq);

        // update the mesh object 
        void updateMesh();

        // make periodic
        void MakePeriodic();

        const Mesh& getMesh() {return mesh_;}

    private:
        std::string FileName_;

        int numPoints_;
        int numEdges_;
        int numSeed_;

        std::vector<Real2> points_;
        std::vector<INT2> constraintIndices_;
        std::vector<Real2> seeds_;
        Real MaxDistanceBetweenEdges_;

        // initialize what is needed for the constrained delaunay triangulation
        CDT cdt_;
        std::vector<Vertex_handle> Vertex_handles_;
        std::vector<Point> ListSeeds_;

        Real aspect_bound_=0.125;
        Real size_bound_=1.0;

        // map from the vertex iterator to the index number 
        std::map<Vertex_handle, int> MapFromVertexToIndex_;

        // output vertices and faces
        int numOutputPoints_;
        int numOutputFaces_;
        int numInputPoints_;
        std::vector<Real3> OutputVertices_;
        std::vector<INT3> OutputFaces_;

        // boolean that signifies whether or not the calculation is for pbc
        bool isPBC_=false;
        Real2 box_length_;

        std::map<INT2, std::vector<int>> MapEdgeToFace_;
        std::vector<bool> boundaryIndicator_;

        Mesh mesh_;
};
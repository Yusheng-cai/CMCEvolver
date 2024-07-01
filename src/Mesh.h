#pragma once

#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "LinAlgTools.h"
#include "tools/InputParser.h"
#include "MeshRefineStrategy.h"
#include "tools/OutputFunction.h"
#include "tools/CommonOperations.h"
#include "tools/Constants.h"
#include "tools/Algorithm.h"
#include "happly.h"
#include "Graph.h"
#include "cmc_surface/fast_rdt.h"
#include "AFP_shapes.h"
#include "tools/CommandLineArguments.h"

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <memory>
#include <numeric>
#include <unordered_map>
#include <iomanip>

#define CGAL_PMP_USE_CERES_SOLVER
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/convex_hull_3.h>

// triangle structs
struct vertex{
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;

    Real3 position_ = {{0,0,0}};
    Real3 normals_ = {{0,0,0}};

    Real& operator[](int i){return position_[i];}
    Real operator[](int i) const {return position_[i];}

    Real3 operator-(const vertex& otherV) const{return position_ - otherV.position_;}
    Real3 operator+(const vertex& otherV) const{return position_ + otherV.position_;}
};

struct triangle{
    using  INT3 = CommonTypes::index3;

    INT3 triangleindices_;

    // each triangle has 3 vertices
    std::array<vertex,3> vertices_;

    int& operator[](int i){return triangleindices_[i];}
    int operator[](int i) const {return triangleindices_[i];}
};

class Mesh
{
    public:
        using Real3 = CommonTypes::Real3;
        using Real  = CommonTypes::Real;
        using refinePtr = std::unique_ptr<MeshRefineStrategy>;
        using INT2  = CommonTypes::index2;
        using INT3  = CommonTypes::index3;
        using double3 = CommonTypes::double3;

        // default constructor and destructor
        Mesh() {};
        ~Mesh(){};
        Mesh(const std::vector<Real3>& v, const std::vector<INT3>& faces);

                                        /*******************************************
                                         * ******** setting Function **************
                                         * ****************************************/
        void setBoxLength(Real3 box) { boxLength_ = box;setPBC(true);}
        void setPBC(bool isPbc) {isPeriodic_=isPbc;}

                                        /********************************************
                                         * ******** accessing Function **************
                                         * *****************************************/
        std::vector<vertex>& accessvertices() {return vertices_;}
        std::vector<double3>& accessverticesPos() {return verticesd_;}
        std::vector<triangle>& accesstriangles() {return triangles_;}
        std::vector<Real>& accessTriangleArea() {return triangleArea_;}
        std::vector<Real3>& accessPerVertexDir1() {return PerVertexdir1_;}
        std::vector<Real3>& accessPerVertexDir2() {return PerVertexdir2_;}

                                        /********************************************
                                         * **********  getter Function **************
                                         * *****************************************/
        const std::vector<vertex>& getvertices() const {return vertices_;}
        std::vector<Real3> getVertexPositions();
        std::vector<INT3>  getFaces();
        const std::vector<triangle>& gettriangles() const {return triangles_;}
        const std::vector<Real>& getTriangleArea() const {return triangleArea_;}
        const std::vector<Real3>& getPerVertexDir1() const {return PerVertexdir1_;}
        const std::vector<Real3>& getPerVertexDir2() const {return PerVertexdir2_;}
        const std::vector<Real3>& getFaceNormals() const {return facetNormals_;}
        const std::vector<std::vector<int>>& getMapVertexToFace() const {return MapVertexIndicesToFaceIndices_;}
        const std::vector<Real3>& getCornerAreas() const {return cornerArea_;}
        int getNumVertices() const {return vertices_.size();}
        int getNumTriangles() const {return triangles_.size();}
        Real3 getBoxLength() const {return boxLength_;}
        const std::map<INT2, std::vector<int>> getMapEdgeToFace() const {return MapEdgeToFace_;}

        // Scale all the vertices by a single number 
        void scaleVertices(Real num);

        // Calculate the corner area according to trimesh code C++
        void CalculateCornerArea();

        // calculate volume
        Real calculateVolume();

        // function that calculates the normals of each of the vertex --> weighted by area 
        void CalcVertexNormalsAreaWeighted();

        // function that calculates the normals of each of the vertex --> not weighted by anything and updates the normals in each vertices
        void CalcVertexNormals();

        // Find the vertex direction 
        void CalcPerVertexDir();

        // update the triangles/vertex/edges if they were to be changed
        void update();

        // shift the com with respect to
        void ShiftCOMWithRespectTo(Real3& COM);

        // calculate COM of an object
        Real3 CalculateCOM();

        // set vertices and triangles
        void SetVerticesAndTriangles(const std::vector<vertex>& vertices, const std::vector<triangle>& triangles);
                                        ///////////////////////////////
                                        /////       PBC stuff   ///////
                                        ///////////////////////////////

        // Check whether or not the mesh is periodic 
        bool isPeriodic() const {return isPeriodic_;}
        // find PBC distance between 2 vertices 
        void getVertexDistance(const vertex& v1, const vertex& v2, Real3& distVec, Real& dist) const;
        void getVertexDistance(const Real3& v1, const Real3& v2, Real3& distVec, Real& dist) const;
        void CalculateShift(const Real3& v1, const Real3& v2, Real3& shiftVec);
        // move a vertex into pbc box
        void MoveVertexIntoBox(const Real3& OldVertPos, Real3& NewVertexPos);
        // find shifted vertex position
        Real3 getShiftedVertexPosition(const vertex& v1, const vertex& v2) const;
        Real3 getShiftedVertexPosition(const Real3& v1, const Real3& v2) const;
        // find the shit that takes a position into the box 
        Real3 getShiftIntoBox(const Real3& v1);
    
    private:
        // vertices and triangles in the mesh 
        std::vector<vertex> vertices_;
        std::vector<triangle> triangles_;

        // store vertex positions in double 
        std::vector<double3> verticesd_;

        // the Areas and normals
        std::vector<Real> triangleArea_;
        std::vector<Real3> cornerArea_;
        std::vector<Real3> facetNormals_;
        std::vector<Real3> vertexNormals_;
        std::vector<Real3> PerVertexdir1_;
        std::vector<Real3> PerVertexdir2_;

        // map the edge to face 
        std::map<INT2, std::vector<int>> MapEdgeToFace_;

        // a map from vertex indices to the face indices  
        std::vector<std::vector<int>> MapVertexIndicesToFaceIndices_;

        // map vertex to edges 
        std::vector<std::vector<INT2>> MapVertexToEdges_;

        // the factor to which all the vertices on the mesh will be multiplied by
        Real factor_=1.0;

        // whether or not the mesh is periodic
        bool isPeriodic_=false;
        Real3 boxLength_;
};


namespace MeshTools
{
    using Real3 = CommonTypes::Real3;
    using Real  = CommonTypes::Real;
    using INT3  = CommonTypes::index3;
    using INT2  = CommonTypes::index2;
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Surface_mesh<K::Point_3>                      M;
    typedef K::FT                                               FT;
    typedef std::array<std::size_t,3>                           CGAL_Polygon;
    typedef std::array<FT, 3>                                   Custom_point;
    typedef M::Vertex_index                                     vertex_descriptor;
    typedef M::Face_index                                       face_descriptor;
    using refineptr= std::unique_ptr<MeshRefineStrategy>;

    namespace PMP = CGAL::Polygon_mesh_processing;

    struct Array_traits
    {
        struct Equal_3
        {
            bool operator()(const Custom_point& p, const Custom_point& q) const {
            return (p == q);
            }
        };
        struct Less_xyz_3
        {
            bool operator()(const Custom_point& p, const Custom_point& q) const {
            return std::lexicographical_compare(p.begin(), p.end(), q.begin(), q.end());
            }
        };
        Equal_3 equal_3_object() const { return Equal_3(); }
        Less_xyz_3 less_xyz_3_object() const { return Less_xyz_3(); }
    };

    bool readPLY(std::string& filename, Mesh& mesh);
    bool readPLYlibr(std::string& filename, Mesh& mesh);

    // write PLY file
    void writePLY(std::string filename, const std::vector<Real3>& Vertices, const std::vector<INT3>& faces, Real factor=1.0);
    void writePLY(std::string filename, const std::vector<Real3>& Vertices, const std::vector<INT3>& faces, const std::vector<Real3>& normals);
    void writePLYRGB(std::string filename, const std::vector<Real3>& Vertices, const std::vector<INT3>& faces, const std::vector<Real3>& RGB);
    void writePLY(std::string filename, const Mesh& mesh);

    // write non pbc mesh 
    void writeNonPBCMesh(std::string filename, Mesh& mesh);
    void writeNonPeriodicTriangleIndices(std::string name, Mesh& mesh);

    // auxiliary functions for writing mesh areas
    void writeMeshArea(std::string filename, Mesh& mesh);
    void writeCuttedMesh(std::string filename, Mesh& mesh, Real3& volume);

    // write STL file
    void writeSTL(std::string filename, Mesh& mesh);
    void writeSTL(std::string filename, const std::vector<Real3>& vertices, const std::vector<INT3>& faces);

    // calculate PBC distance
    Real3 calculateShift(const Real3& vec1, const Real3& vec2, const Real3& boxLength);

    // calculate pbc distance 
    void calculateDistance(const Real3& vec1, const Real3& vec2, const Real3& boxlength, Real3& distance, Real& distsq);

    // find if an edge is periodic 
    bool isPeriodicEdge(const Real3& vec1, const Real3& vec2, Real3& newarr, const Real3& boxLength);

    // map vertices to faces
    void MapVerticesToFaces(const Mesh& mesh, std::vector<std::vector<int>>& map);

    // calculate triangle areas and Face normals 
    void CalculateTriangleAreasAndFaceNormals(Mesh& mesh, std::vector<Real>& Areas, std::vector<Real3>& Normals);

    // find vertex neighbors 
    void CalculateVertexNeighbors(const Mesh& mesh, std::vector<std::vector<int>>& neighborIndices);

    // map Edge to faces
    void MapEdgeToFace(const Mesh& mesh, std::map<INT2,std::vector<int>>& magEdgeToFace, std::vector<std::vector<INT2>>& MapVertexToEdge, bool assert=true);

    void MapEdgeToFace(const Mesh& mesh, std::map<INT2,std::vector<int>>& magEdgeToFace, bool assert=true);

    // map edge to opposing vertices --> needs boundary indicators  
    void MapEdgeToOpposingVertices(Mesh& mesh, std::map<INT2, std::vector<int>>& mapEdgeToFace,std::map<INT2, std::vector<int>>& MapEdgeToOppoVertices);

    // find boundary vertices 
    void CalculateBoundaryVertices(const Mesh& mesh, const std::map<INT2, std::vector<int>>& mapEdgeToFace, std::vector<bool>& boundaryIndicator);
    void CalculateBoundaryVertices(const Mesh& mesh, std::vector<bool>& boundaryIndicator);
    void CalculateBoundaryVerticesIndex(const Mesh& mesh, std::vector<int>& boundaryIndices, bool order=false, Real3 center={0,0,0});

    // check if point is on boundary
    bool IsBoundary(int Index, const std::vector<bool>& boundaryIndicator);

    // convert a PBC mesh to non PBC mesh --> usually for visualization purposes
    void ConvertToNonPBCMesh(Mesh& mesh, std::vector<Real3>& vertices, std::vector<INT3>& faces, bool AddNewTriangles=false);
    void ConvertToNonPBCMesh(Mesh& mesh, bool AddNewTriangles=false);

    // check if a particular triangle is periodic
    bool IsPeriodicTriangle(std::vector<vertex>& Vertices,INT3& face, Real3 BoxLength);
    bool IsPeriodicTriangle(std::vector<vertex>& Vertices,INT3& face, Real3 BoxLength, std::map<INT2,bool>& mapEdge);
    bool IsPeriodicTriangle(const Mesh& mesh, int faceindex);
    bool IsPeriodicTriangle(const Mesh& mesh, int faceindex, std::map<INT2,bool>& mapEdge);

    bool FindPeriodicShiftIndex(std::map<INT2,bool>& mapEdge, int& index);

    // shift a triangle 
    void ShiftPeriodicTriangle(const std::vector<vertex>& Vertices, const INT3& faces, Real3 BoxLength, Real3& A, Real3& B, Real3& C);

    // shift a triangle based on a certain point
    void ShiftPeriodicTriangle(const std::vector<vertex>& Vertices, const INT3& faces, Real3 BoxLength, const Real3& point, Real3& A, Real3& B, Real3& C);

    // make an Edge, edge is simply the 2 indices of {{minIndex, maxIndex}}
    INT2 makeEdge(int i, int j);

    // calculate the corner area 
    void CalculateCornerArea(Mesh& mesh, std::vector<Real3>& CornerArea, std::vector<Real>& VertexArea);

    // cut mesh and give a new mesh 
    void CutMesh(Mesh& mesh, std::vector<INT3>& face, std::vector<Real3>& vertices, Real3 volume, bool below=false);

    // cut mesh and give new mesh
    void CutMesh(Mesh& mesh, Real3& volume, bool below=false);

    // Moller Trumbore Ray-Triangle intersection method 
    bool MTRayTriangleIntersection(Real3& A, Real3& B, Real3& C, Real3& O, Real3& D, Real& t, Real& u, Real& v);

    // check degenerate triangles 
    void CheckDegenerateTriangle(Mesh& mesh, \
                                             std::vector<int>& MergeFaces,\
                                             std::vector<INT2>& MergeVertices);

    // decimate degenerate triangles
    bool decimateDegenerateTriangle(Mesh& mesh);

    // calculate barycenter of all the boundary vertices 
    void CalculateBoundaryBarycenter(Mesh& mesh, std::vector<bool>& boundaryIndicator, Real3& barycenter);

    // correct  Mesh --> Many times, we need to cut certain triangles out for certain reasons
    // this function takes care of rearranging the mesh vertices and triangles
    void CorrectMesh(Mesh& mesh, std::vector<int>& FaceIndices);

    // check if a triangle is isolated --> all of its edges are boundaries
    bool IsIsolatedFace(Mesh& mesh, int faceIndex, const std::map<INT2, std::vector<int>>& mapEdgeToFace);

    // check if a triangle is teeth like --> 2 of its edges are boundaries
    bool IsTeethlikeFace(Mesh& mesh, int faceIndex, const std::map<INT2, std::vector<int>>& mapEdgeToFace);

    // regenerate mesh after some faces are cut
    void ReconstructMeshAfterFaceCut(Mesh& mesh);

    // remove vertices that does not meet the minimum criteria of neighbors
    void RemoveMinimumNeighbors(Mesh& mesh, int num_search, int min_num_neighbors);

    // get a sense of the triangular angles 
    void FindTriangleAngles(const Mesh& mesh, std::vector<Real>& angles);

    // return distribution of side lengths
    void FindSideLengths(const Mesh& mesh, std::vector<Real>& SideLength);

    // remove isolated vertices
    std::vector<int> RemoveIsolatedVertices(Mesh& mesh);

    // remove isolated faces 
    void RemoveIsolatedFaces(Mesh& mesh);

    // remove duplicate triangles 
    void RemoveDuplicatedFaces(Mesh& mesh);

    // change winding order of a mesh
    void ChangeWindingOrder(Mesh& m);
    void ChangeWindingOrder(Mesh& m, int num);

    // mesh plane clipping 
    void MeshPlaneClipping(Mesh& m, Real3& point, Real3& normal);

    // mesh plane 
    void TriangleCases(std::vector<INT3>& signs, std::vector<bool>& basic, std::vector<bool>& one_vertex, std::vector<bool>& one_edge);
    void MeshPlaneIntersection(Mesh& m, Real3& point, Real3& normal);

    // // calculate the volume enclosed underneath an interface
    Real CalculateVolumeUnderneath(Mesh& m, int projected_plane=2);
    
    void CalculateCotangentWeights(Mesh& m, const std::vector<std::vector<int>>& neighborIndices, const std::map<INT2, std::vector<int>>& MapEdgeToFace, const std::map<INT2, std::vector<int>>& MapEdgeToOpposingVerts, std::vector<Real3>& dAdpi);
    void CalculateAreaDerivatives(Mesh& m, std::vector<Real3>& dAdpi);

    void CalculateVolumeDerivatives(Mesh& m, const std::vector<std::vector<int>>& MapVertexToFace, std::vector<Real3>& VolumeDerivatives, Real3 shift={0,0,0});
    void CalculateVolumeDerivatives(Mesh& m, std::vector<Real3>& VolumeDerivatives, Real3 shift={0,0,0});

    Real CalculateVolumeDivergenceTheorem(Mesh& m, const std::vector<Real>& vecArea, const std::vector<Real3>& Normal);

    Real CalculateArea(Mesh& m, std::vector<Real>& vecArea, std::vector<Real3>& Normal);

    void CalculateAV(Mesh&, Real& A, Real& V);

    // optimize the mesh with centroidal voronoi calculation
    void CVT_optimize_Mesh(Mesh& m, Real volume_weight=0.0, int nb_iterations=100);

    // CGAL optimize mesh
    void CGAL_optimize_Mesh(Mesh& m, int nb_iterations, Real degree, bool use_restriction=false);

    void MakePBCMesh(Mesh& m);
    
    void ShiftPBCMesh(Mesh& m, Real3& shift);

    // convert from my mesh to cgal mesh
    void ConvertToCGALMesh(Mesh& m, M& cgal_m);

    // find the u,v coordinate of the boundary vertices
    void FindBoundaryUV(Mesh& m, std::vector<Real>& ulist, std::vector<Real>& vlist, const std::vector<bool>& BoundaryIndicator, AFP_shape* shape);
    void FindBoundaryUV(Mesh& m, std::vector<Real>& ulist, std::vector<Real>& vlist, std::vector<int>& BoundaryIndices, AFP_shape* shape, bool order=false);

    // calculate drduv
    void CalculatedrdUV(Mesh& m, AFP_shape* shape, std::vector<int>& BoundaryIndices, std::vector<Real>& ulist, std::vector<Real>& vlist, std::vector<Real3>& drdu, std::vector<Real3>& drdv, bool useNumerical=true);
    void CalculatedrdUV(Mesh& m, AFP_shape* shape, std::vector<Real3>& drdu, std::vector<Real3>& drdv, bool useNumerical=true);

    // calculate dAnbsduv
    void CalculatedAVnbsdUV(Mesh& m,AFP_shape* shape, std::vector<Real2>& dAnbsduv, std::vector<Real2>& dVnbsduv, bool useNumerical=true, Real3 Vshift={0,0,0});
    void CalculatedAVnbsdUV(Mesh& m,AFP_shape* shape, std::vector<int>& BoundaryIndices, std::vector<Real2>& dAnbsduv, std::vector<Real2>& dVnbsduv, bool useNumerical=true, Real3 Vshift={0,0,0});
    void CalculatedAVnbsdUV(Mesh& m,AFP_shape* shape, std::vector<int>& BoundaryIndices, std::vector<Real3>& drdu, std::vector<Real3>& drdv, std::vector<Real2>& dAnbsduv, std::vector<Real2>& dVnbsduv, bool useNumerical=true, Real3 Vshift={0,0,0});
    void CalculatedAVnbsdUV(Mesh& m,AFP_shape* shape, std::vector<int>& BoundaryIndices, std::vector<Real>& ulist, std::vector<Real>& vlist, std::vector<Real3>& drdu, std::vector<Real3>& drdv, std::vector<Real2>& dAnbsduv, std::vector<Real2>& dVnbsduv, bool useNumerical=true, Real3 Vshift={0,0,0});

    // helper function which helps read the inputs for shape
    std::unique_ptr<AFP_shape> ReadAFPShape(CommandLineArguments& cmd);

    refineptr ReadInterfacialMin(CommandLineArguments& cmd);
    refineptr ReadInterfacialMinBoundary(CommandLineArguments& cmd);

    void CalculateContactAngle(Mesh& m, AFP_shape* s, std::vector<Real>& ca);
    void CalculateContactAngleDerivative(Mesh& m, AFP_shape* s, std::vector<Real>& ca, Real k0, Real3 Volume_shift={0,0,0}, bool use_Numerical=true);

    void CalculateAVnbs(Mesh& m, AFP_shape* s, std::vector<int>& BoundaryIndices, std::vector<Real>& ulist,\
                        std::vector<Real>& vlist, Real& A, Real& V, int v_num=1000, bool useNumerical=true, Real3 Vshift={0,0,0});
    void CalculateAVnbs(Mesh& m, AFP_shape* s, Real& A, Real& V, int v_num=1000, bool useNumerical=true, Real3 Vshift={0,0,0});

    Real CalculateVnbsUnderneath(Mesh& m, AFP_shape* s, std::vector<int>& BoundaryIndices, std::vector<Real>& ulist, \
                                  std::vector<Real>& vlist, int projected_plane=2, int v_num=1000, bool useNumerical=true);
    Real CalculateVnbsUnderneath(Mesh& m, AFP_shape* s, int projected_plane=2, int v_num=1000, bool useNumerical=true);

    Real CalculateBoundaryAverageHeight(Mesh& m);

    Real CalculateMaxCurvature(Mesh& m, const std::vector<int>& BoundaryIndices);

    Real CalculateEta(Mesh& m, Real& Pbar, Real3 Boundary_center={0,0,0});

    Real CalculatePbar(Mesh& m, const std::vector<int>& BoundaryIndices, std::vector<Real3>& tdr);
};
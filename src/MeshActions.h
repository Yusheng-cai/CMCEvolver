#pragma once

#include <vector>
#include <string>
#include <array>
#include <functional>
#include <memory>

#ifdef IGL_ENABLED
#include "igl/opengl/glfw/Viewer.h"
#include "igl/read_triangle_mesh.h"
#endif

#include "Eigen/Core"

#include "tools/Assert.h"
#include "MeshGen2d.h"
#include "MeshPlaneIntersection.hpp"
#include "tools/CommandLineArguments.h"
#include "Mesh.h"
#include "tools/CommonTypes.h"
#include "Curvature.h"
#include "tools/InputParser.h"
#include "MeshRefineStrategy.h"
#include "CurvatureCurveFit.h"
#include "CurvatureJetFit.h"
#include "tools/InputParser.h"
#include "tools/CommonOperations.h"
#include "CurvatureEvolution.h"
#include "delaunator.hpp"
#include "tools/Algorithm.h"
#include "LinAlgTools.h"
#include "Bin.h"
#include "Graph.h"
#include "ShortEdgeRemoval.h"
#include "ObtuseTriangleRemoval.h"
#include "LongEdgeRemoval.h"
#include "AFP_shapes.h"
#include "cmc_surface/fast_rdt.h"
#include "InterfacialFE_minimization.h"

#define CGAL_PMP_USE_CERES_SOLVER
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>


namespace MeshActions
{
    using INT3 = std::array<int,3>;
    using INT2 = std::array<int,2>;
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;
    using Real2= CommonTypes::Real2;
    using curveptr = std::unique_ptr<Curvature>;
    using refineptr= std::unique_ptr<MeshRefineStrategy>;
    using double3  = CommonTypes::double3;

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Surface_mesh<K::Point_3>                      M;
    namespace PMP = CGAL::Polygon_mesh_processing;

    // Calculate mesh derivatives
    void MeshDerivatives(CommandLineArguments& cmd);

    // BoundaryRefinement
    void RefineBoundary(CommandLineArguments& cmd);

    // Meshify boundary
    void MeshifyShape(CommandLineArguments& cmd);

    // translate a mesh
    void TranslateMesh(CommandLineArguments& cmd);

    // calculate the mesh of a mesh using the method of curve fitting (second fundamental form)
    void CurveFit(CommandLineArguments& cmd);

    // calculate the curvature of a mesh using derivative
    void CurvatureDerivative(CommandLineArguments& cmd);

    // quadratic curve fit
    void QuadraticCurveFit(CommandLineArguments& cmd);

    // calculate the curvature of a mesh using the method of jet fitting
    void JetFit(CommandLineArguments& cmd);

    // calculate the curvature of a mesh using the method of tensor
    void TensorFit(CommandLineArguments& cmd);

    // calculate the curvature of a mesh using the method of FDM
    void FDMFit(CommandLineArguments& cmd);

    // smooth a mesh using curvature flow (implicit fairing)
    void CurvatureFlow(CommandLineArguments& cmd);

    // function that finds the indices of non periodic triangles in a periodic mesh 
    void FindNonPBCTriangles(CommandLineArguments& cmd);

    void OptimizeMesh(CommandLineArguments& cmd);

    // cut the mesh 
    void CutMesh(CommandLineArguments& cmd);

    // clip the mesh
    void ClipMesh(CommandLineArguments& cmd);

    // convert To STL file
    void ConvertToStL(CommandLineArguments& cmd);

    // convert a pbc mesh to non pbc mesh
    void ConvertToNonPBCMesh(CommandLineArguments& cmd);

    // convert a non pbc mesh to pbc mesh
    void ConvertToPBCMesh(CommandLineArguments& cmd);
    
    // scale mesh by some number 
    void ScaleMesh(CommandLineArguments& cmd);

    // put color onto vertex based on curvature value
    void ColorVertex(CommandLineArguments& cmd);

    // project 3d curvature onto 2d using Moller Trumbore Ray-Triangle intersection algorithm
    void Project3dCurvature(CommandLineArguments& cmd);

    // function that constructs a new mesh such that none of the vertices is within a certain distance of the reference mesh 
    void MeshDistanceCutoff(CommandLineArguments& cmd);

    // flatten a certain dimension
    void FlattenMeshDimension(CommandLineArguments& cmd);

    // Find the minimum distance of a mesh against another reference mesh
    void MinimumMeshDistance(CommandLineArguments& cmd);

    // find mesh plane intersection
    void MeshPlaneIntersection(CommandLineArguments& cmd);

    // find vertex neighbors 
    void FindVertexNeighbors(CommandLineArguments& cmd);

    // find boundary vertices 
    void FindBoundaryVertices(CommandLineArguments& cmd);

    // find the 2d projection of a 3d mesh viewed from a certain angle 
    void Project3dMesh(CommandLineArguments& cmd);

    // cut overlapped regions of 2 proteins based on Ray triangle intersection
    void CutOverlappedRegion(CommandLineArguments& cmd);

    // function that decimates degenerate triangles 
    void DecimateDegenerateTriangles(CommandLineArguments& cmd);

    // curvature evolution --> only use curvefit
    void CurvatureEvolution1(CommandLineArguments& cmd);

    // curvature evolution but using interfacial Fe
    void InterfacialFE_min(CommandLineArguments& cmd);

    void InterfacialFE_min_boundary(CommandLineArguments& cmd);

    void InterfacialFE_min_boundary_L1_constraint(CommandLineArguments& cmd);

    // function that calculates the contact angle
    void CalculateContactAngle(CommandLineArguments& cmd);

    // find isolated triangle 
    void FindIsolatedFace(CommandLineArguments& cmd);

    // cut out teeth like faces
    void CutTeethlikeFace(CommandLineArguments& cmd);

    // replicate periodic mesh --> make it
    void ReplicatePeriodicMesh(CommandLineArguments& cmd);

    // find the distribution of the triangle angles
    void TriangleAngleDistribution(CommandLineArguments& cmd);

    // find the distribution of the side lengths
    void SideLengthsDistribution(CommandLineArguments& cmd);

    // clean up mesh
    void MeshCleanup(CommandLineArguments& cmd);

    // flatten a mesh --> convert one of the indices to be all 0
    void FlattenMesh(CommandLineArguments& cmd);

    // generate a mesh --> conforming triangulations
    // https://doc.cgal.org/latest/Mesh_2/index.html
    void ConformingTriangulations(CommandLineArguments& cmd);

    // split long edges
    void SplitLongEdges(CommandLineArguments& cmd);

    // perform iterative closest point to shift a mesh with respect to another
    void ShiftMeshWithRef(CommandLineArguments& cmd);

    void ShiftPBCMesh(CommandLineArguments& cmd);

    void ViewMeshWithData(CommandLineArguments& cmd);

    // view mesh 
    void ViewMesh(CommandLineArguments& cmd);

    // Calculate Distance between 2 meshes using Ray-Triangle intersection algorithm
    void DistanceBetweenMeshesMT(CommandLineArguments& cmd);

    // change winding order
    void ChangeMeshWindingOrder(CommandLineArguments& cmd);

    // find face normals.size()
    void FindFaceNormals(CommandLineArguments& cmd);

    // calculate the surface properties at which the interface is pinned
    void calculateSurfaceProperties(CommandLineArguments& cmd);

    // calculate the surface area
    void calculateSurfaceArea(CommandLineArguments& cmd);

    // calculate volume
    void calculateInterfaceVolume(CommandLineArguments& cmd);

    // calculate volume underneath
    void calculateInterfaceVolumeUnderneath(CommandLineArguments& cmd);

    // CVT optimization
    void CVT_Mesh_optimization(CommandLineArguments& cmd);

    // Calculate Anbs and Vnbs
    void Mesh_AVnbs(CommandLineArguments& cmd);

    // Calculate Eta
    void Mesh_Eta(CommandLineArguments& cmd);
};
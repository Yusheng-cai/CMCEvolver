add_library(CEsrc SHARED)

set(CEsrc_CPU_source 
    "Bin.cpp"
    "LinAlgTools.cpp"
    "Mesh.cpp"
    "Curvature.cpp"
    "CurvatureFDM.cpp"
    "CurvatureTensor.cpp"
    "Graph.cpp"
    "CurvatureCurveFit.cpp"
    "MeshRefineStrategy.cpp"
    "CurvatureJetFit.cpp"
    "UmbrellaSmoothing.cpp"
    "MeshCurvatureflow.cpp"
    "CurvatureEvolution.cpp"
    "Registry.cpp"
    "MeshGen2d.cpp"
    "MeshActions.cpp"
    "ShortEdgeRemoval.cpp"
    "ObtuseTriangleRemoval.cpp"
    "LongEdgeRemoval.cpp"
    "CurvatureQuadraticFit.cpp"
    "AFP_shapes.cpp"
    "fsolve.cpp"
    "fsolve_wrapper.cpp"
    "InterfacialFE_minimization.cpp"
)


target_sources(CEsrc PRIVATE ${CEsrc_CPU_source})

set(CEsrc_HEADER
    "Bin.h"
    "LinAlgTools.h"
    "Mesh.h"
    "Curvature.h"
    "CurvatureFDM.h"
    "CurvatureTensor.h"
    "Graph.h"
    "CurvatureCurveFit.h"
    "MeshRefineStrategy.h"
    "CurvatureJetFit.h"
    "UmbrellaSmoothing.h"
    "MeshCurvatureflow.h"
    "CurvatureEvolution.h"
    "Registry.h"
    "MeshGen2d.h"
    "MeshActions.h"
    "MeshPlaneIntersection.hpp"
    "ShortEdgeRemoval.h"
    "ObtuseTriangleRemoval.h"
    "LongEdgeRemoval.h"
    "CurvatureQuadraticFit.hpp"
    "AFP_shapes.h"
    "fsolve.hpp"
    "fsolve_wrapper.hpp"
    "InterfacialFE_minimization.h"
)

target_include_directories(CEsrc
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

if (CGAL_FOUND)
   target_include_directories(CEsrc PUBLIC ${CGAL_INCLUDE_DIRS})
endif()

if (Ceres_FOUND)
    target_include_directories(CEsrc PUBLIC ${CERES_INCLUDE_DIRS})
endif()

if (ENABLE_IGL)
    target_link_libraries(CEsrc PUBLIC CEtools CEparallel CGAL::CGAL igl::glfw ${CERES_LIBRARIES})
else()
    target_link_libraries(CEsrc PUBLIC CEtools CEparallel CGAL::CGAL ${CERES_LIBRARIES})
endif()


include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(CEsrc 
    INCLUDE_SUBDIR CEsrc
    INCLUDES ${CEsrc_HEADER})

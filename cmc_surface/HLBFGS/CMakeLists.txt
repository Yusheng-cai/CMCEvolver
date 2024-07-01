add_library(CEHLBFGS SHARED)

target_sources(CEHLBFGS PRIVATE
    "HLBFGS_BLAS.cpp"
    "HLBFGS.cpp"
    "ICFS.cpp"
    "LineSearch.cpp"
    "Lite_Sparse_Matrix.cpp"
)

set(CEHLBFGS_HEADER
    "HLBFGS_BLAS.h"
    "HLBFGS.h"
    "ICFS.h"
    "LineSearch.h"
    "Lite_Sparse_Matrix.h"
    "Sparse_Entry.h"
)

target_include_directories(CEHLBFGS
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(CEHLBFGS
    INCLUDE_SUBDIR CEHLBFGS
    INCLUDES ${CEHLBFGS_HEADER})
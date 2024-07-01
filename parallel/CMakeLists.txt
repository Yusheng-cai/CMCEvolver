add_library(CEparallel SHARED)

target_sources(CEparallel PRIVATE
    "OpenMP.cpp"
)

set(CEPARALLEL_HEADER
    "MPI_environment.h"
    "OpenMP_buffer.h"
    "OpenMP.h")

target_include_directories(CEparallel PUBLIC 
                         "${CMAKE_SOURCE_DIR}")

# target_compile_options(parallel PRIVATE -fPIC)

message("Required package of OpenMP is ${OpenMP_CXX_FOUND}")

# resource: https://gitlab.kitware.com/cmake/cmake/-/issues/21818
if(OpenMP_CXX_FOUND)
    target_link_libraries(CEparallel PUBLIC OpenMP::OpenMP_CXX)
endif()

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(CEparallel 
    INCLUDE_SUBDIR CEparallel 
    INCLUDES ${CEPARALLEL_HEADER})
    
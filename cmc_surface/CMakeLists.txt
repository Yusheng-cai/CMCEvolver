add_library(CEcmc_surface SHARED)

target_sources(CEcmc_surface PRIVATE
    "fast_rdt.cpp"
)

set(CECMC_SURFACE_HEADER
    "fast_rdt.h"
    "geex_utils.h"
    "edgemesh.h")

target_include_directories(CEcmc_surface PUBLIC 
                         "${CMAKE_SOURCE_DIR}")

add_subdirectory(HLBFGS)
target_link_libraries(CEcmc_surface PUBLIC CEHLBFGS)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(CEcmc_surface 
    INCLUDE_SUBDIR CEcmc_surface
    INCLUDES ${CECMC_SURFACE_HEADER})
    
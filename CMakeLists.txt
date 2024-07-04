cmake_minimum_required(VERSION 3.18.1)

# Include the cmake scripts directory
LIST( PREPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/
)

# set an option for using CUDA
option(ENABLE_IGL, "Using igl library" OFF)

# Require C++11
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DMY_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math")

project(CMCEvolver LANGUAGES C CXX)

# libigl
if (ENABLE_IGL)
  include(libigl)
  igl_include(glfw)
  add_compile_options(-DIGL_ENABLED)
endif()

# Requires OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(CGAL NO_DEFAULT_PATH)
if (CGAL_FOUND)
  message("CGAL FOUND.")
  message("CGAL INCLUDE DIR = ${CGAL_INCLUDE_DIRS}")
  message("CGAL LIBRARIES = ${CGAL_LIBRARY}")
  message("CGAL VERSION=${CGAL_VERSION}")
endif()

find_package(Ceres NO_DEFAULT_PATH)
if (Ceres_FOUND)
  message("Ceres FOUND")
  message("Ceres INCLUDE DIR = ${CERES_INCLUDE_DIRS}")
  message("Ceres LIBRARIES = ${CERES_LIBRARIES}")
endif()


# install the libraries as well as the executables in lib & bin folder in build/
# This has to be after project is defined
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")

enable_testing()

add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(parallel)
add_subdirectory(test)
add_subdirectory(cmc_surface)

add_executable(CMCEvolver mesh_op.cpp)

target_include_directories(CMCEvolver PUBLIC ${CMAKE_SOURCE_DIR})

target_link_libraries(CMCEvolver PUBLIC CEsrc CEtools CEparallel CEcmc_surface)

install(TARGETS CMCEvolver DESTINATION bin)

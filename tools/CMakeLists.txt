add_library(CEtools SHARED)

target_sources(CEtools PRIVATE
    "RandomNumber.cpp"
    "InputParser.cpp"
    "CommandLineArguments.cpp"
    "FileSystem.cpp"
    "OutputFunction.cpp"
    "Algorithm.cpp"
)

set(CETOOLS_HEADER
    "Assert.h"
    "CommonTypes.h"
    "GenericFactory.h"
    "InputParser.h"
    "RandomNumber.h"
    "CommandLineArguments.h"
    "Registry.h"
    "FileSystem.h"
    "OutputFunction.h"
    "Algorithm.h"
)

target_include_directories(CEtools
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)
# target_compile_options(tools PRIVATE -fPIC)


include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(CEtools 
    INCLUDE_SUBDIR CEtools
    INCLUDES ${CETOOLS_HEADER})


# CMake version
cmake_minimum_required(VERSION 3.21)
#set the project name and version
project(reads2graph CXX)
find_package(seqan3 REQUIRED)
find_package(sharg REQUIRED)
find_package(boost REQUIRED)
# Add your source files
set(SOURCES
    LoggingLevels.cpp
    Utils.cpp
    reads2graph.cpp
    ReadWrite.cpp
    MinimizerGenerator.cpp
    GraphConstructor.cpp
    OMH.cpp
)
#include_directories
include_directories(include $ENV{CONDA_PREFIX}/include)
link_directories($ENV{CONDA_PREFIX}/lib)

set(CMAKE_CXX_FLAGS "${CAMKE_CXX_FLAGS} -std=c++20 -fopenmp -static-libgcc -static-libstdc++ -isystem $ENV{CONDA_PREFIX}/include/seqan3/submodules/sdsl-lite/include -isystem $ENV{CONDA_PREFIX}/include/seqan3/submodules/cereal/include -fconcepts -g -Wall -lyaml-cpp -Wno-interference-size") #-fsanitize=address

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
aux_source_directory(${PROJECT_SOURCE_DIR}/src c_files)
add_executable(${PROJECT_NAME} ${c_files})
# Link against the required libraries
target_link_libraries(${PROJECT_NAME} -Wl,-Bstatic xxhash -Wl,-Bdynamic)


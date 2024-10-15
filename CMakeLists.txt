# CMake version
cmake_minimum_required(VERSION 3.21)
#set the project name and version
project(reads2graph)
find_package(seqan3 REQUIRED)
find_package(sharg REQUIRED)

# Find and set up XXHash
find_library(XXHASH_LIB NAMES xxhash PATHS $ENV{CONDA_PREFIX}/lib REQUIRED)
find_path(XXHASH_INCLUDE NAMES xxhash.h PATHS $ENV{CONDA_PREFIX}/include REQUIRED)

include_directories(${XXHASH_INCLUDE})

# Add your source files
set(SOURCES
    LoggingLevels.cpp
    Utils.cpp
    reads2graph.cpp
    ReadWrite.cpp
    MinimizerGenerator.cpp
    GraphConstructor.cpp
    gOMH.cpp
    omh.cpp
)
#include_directories
include_directories(include)

set(CMAKE_CXX_FLAGS "${CAMKE_CXX_FLAGS} -std=c++20 -fopenmp -static-libgcc -static-libstdc++ -isystem $ENV{CONDA_PREFIX}/include -isystem $ENV{CONDA_PREFIX}/include/seqan3/submodules/sdsl-lite/include -isystem $ENV{CONDA_PREFIX}/include/seqan3/submodules/cereal/include -fconcepts -g -Wall -lyaml-cpp -Wno-interference-size -Wno-maybe-uninitialized") #-fsanitize=address

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
aux_source_directory(${PROJECT_SOURCE_DIR}/src c_files)
add_executable(${PROJECT_NAME} ${c_files})
target_link_libraries(${PROJECT_NAME} PRIVATE ${XXHASH_LIB})
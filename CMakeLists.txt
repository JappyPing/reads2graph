# CMake version
cmake_minimum_required(VERSION 3.20)

# set(CMAKE_CXX_COMPILER "/home/pping/miniconda3/envs/gcc_env/bin/g++")

#set the project name and version
project(reads2graph CXX)

# # add seqan3 to search path
# list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan3/build_system")
# # add the Sharg Parser to search path
# list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../sharg-parser/build_system")
# Find SeqAn3 (adjust the path accordingly)
find_package(seqan3 REQUIRED)
find_package(sharg REQUIRED)
# find_package(boost REQUIRED)

# Add your source files
set(SOURCES
    LoggingLevels.cpp
    Utils.cpp
    reads2graph.cpp
    ReadWrite.cpp
    GraphManager.cpp
    MinimizerGenerator.cpp
    EdgeConstructor.cpp
    PairWiseEditDis.cpp
    OMH.cpp
)

#include_directories
include_directories(include $ENV{CONDA_PREFIX}/include)
link_directories($ENV{CONDA_PREFIX}/lib)

# # Set include directories for the project
# include_directories(
#     ${SEQAN3_INCLUDE_DIRS}  # Include directories for SeqAn3
#     ${sharg_INCLUDE_DIRS}   # Include directories for sharg
#     include                 # Your project's include directory
#     $ENV{CONDA_PREFIX}/include/lemon  # Explicitly include Lemon's include directory
# )

# # Set link directories for the project
# link_directories(
#     ${SEQAN3_LIB_DIRS}      # Library directories for SeqAn3
#     ${sharg_LIB_DIRS}       # Library directories for sharg
#     $ENV{CONDA_PREFIX}/lib  # Additional library directory if needed
# )



# include_directories(/home/pping/Data/Repo/reads2graph/include)
# include_directories(/home/pping/miniconda3/envs/gcc_env/include)
# link_directories(/home/pping/miniconda3/envs/gcc_env/lib)

# set(ENVS_INCLUDE $ENV{CONDA_PREFIX}/include)
# set(ENVS_LIB $ENV{CONDA_PREFIX}/lib)
# set(ENVS_EXECUTABLE $ENV{CONDA_PREFIX}/bin)

#Set compilling parameters
# set(CMAKE_CXX_STANDARD 20)
# set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_CXX_EXTENSIONS ON)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CAMKE_CXX_FLAGS} -std=c++20 -fopenmp -isystem /home/pping/miniconda3/envs/gcc_env/include/seqan3/submodules/sdsl-lite/include -isystem /home/pping/miniconda3/envs/gcc_env/include/seqan3/submodules/cereal/include -fconcepts -g -Wall -lyaml-cpp -lxxhash -Wno-interference-size") #-fsanitize=address

# set(CMAKE_CXX_FLAGS "${CAMKE_CXX_FLAGS} -std=c++20 -isystem ${CONDA_PREFIX}/envs/gcc_env/include/seqan3/submodules/sdsl-lite/include -isystem ${CONDA_PREFIX}/envs/gcc_env/include/seqan3/submodules/cereal/include -pthread -fconcepts -g -Wall -lyaml-cpp ") #-fsanitize=address

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

aux_source_directory(${PROJECT_SOURCE_DIR}/src c_files)
# target_link_libraries(${PROJECT_NAME} ${CONDA_PREFIX}/lib/libsnap.so)

# add the executable
add_executable(${PROJECT_NAME} ${c_files})
# target_link_libraries(${PROJECT_NAME} $ENV{CONDA_PREFIX}/lib/libigraph.so)
# Link against the required libraries
# target_link_libraries(${PROJECT_NAME}
#     $ENV{CONDA_PREFIX}/lib/libigraph.so
#     seqan3::seqan3        # Link against SeqAn3
#     sharg::sharg          # Link against sharg
#     lemon                  # Link against Lemon
#     # Add other libraries if needed
# )
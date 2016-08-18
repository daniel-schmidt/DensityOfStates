cmake_minimum_required( VERSION 2.8 )
project ( DensityOfStates )
include(CTest)

set( DensityOfStates_VERSION_MAJOR 0 )
set( DensityOfStates_VERSION_MINOR 1 )

set(CMAKE_CXX_FLAGS "-O3 -march=native -std=c++11 -lm -Wall -g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -march=native -std=c++11 -lm -Wall -O0 -g")

include_directories($ENV{EIGEN3_INCLUDE_DIR})
include_directories("FermiOwnCommon/src/Math")
include_directories("FermiOwnCommon/src/Data")

add_subdirectory("FermiOwnCommon")
#add_subdirectory("src")
add_subdirectory("apps")
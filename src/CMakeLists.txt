file(GLOB SOURCES "*.cpp")
add_library( DensityOfStates ${SOURCES} )
target_link_libraries( DensityOfStates AlgorithmsCommon )
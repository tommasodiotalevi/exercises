find_package(ROOT REQUIRED COMPONENTS Hist Graf Graf3d Gpad)

add_subdirectory(src)

add_executable(boltzmannDistribution.out boltzmannDistribution.cpp)
target_link_libraries(boltzmannDistribution.out PRIVATE ROOT::Hist ROOT::Graf ROOT::Graf3d ROOT::Gpad)
target_link_libraries(boltzmannDistribution.out PRIVATE boltzmann_lib)
target_include_directories(boltzmannDistribution.out PRIVATE ${DOCTEST_INCLUDE_DIR})
add_test(NAME boltzmannDistribution.out COMMAND boltzmannDistribution.out 1000 1000,10000,100000)
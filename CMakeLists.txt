cmake_minimum_required(VERSION 3.16)
project(SlopeStabilitySim VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(WIN32)
    add_definitions(-DWIN32_LEAN_AND_MEAN -DNOMINMAX)
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost 1.71 REQUIRED COMPONENTS program_options filesystem)

add_library(slope_lib STATIC
    src/mesh/MeshGenerator.cpp
    src/soil/SoilWaterCurve.cpp
    src/coupling/CouplingAssembler.cpp
    src/solver/NewtonRaphsonSolver.cpp
    src/export/VTKExporter.cpp
    src/stability/SlopeStability.cpp
    src/cli/CLIParser.cpp
)

target_include_directories(slope_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh
    ${CMAKE_CURRENT_SOURCE_DIR}/src/soil
    ${CMAKE_CURRENT_SOURCE_DIR}/src/coupling
    ${CMAKE_CURRENT_SOURCE_DIR}/src/solver
    ${CMAKE_CURRENT_SOURCE_DIR}/src/export
    ${CMAKE_CURRENT_SOURCE_DIR}/src/stability
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cli
)

target_link_libraries(slope_lib PUBLIC
    Eigen3::Eigen
    Boost::program_options
    Boost::filesystem
)

add_executable(slope_sim src/main.cpp)
target_link_libraries(slope_sim PRIVATE slope_lib)

if(BUILD_TESTING)
    enable_testing()
    find_package(GTest REQUIRED)
    add_subdirectory(tests)
endif()

install(TARGETS slope_sim RUNTIME DESTINATION bin)
install(TARGETS slope_lib ARCHIVE DESTINATION lib)
install(DIRECTORY data/ DESTINATION share/slope_sim/data)

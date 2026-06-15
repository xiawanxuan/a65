find_package(GTest REQUIRED)

add_executable(slope_tests
    test_mesh.cpp
    test_soil_water_curve.cpp
    test_coupling.cpp
    test_solver.cpp
    test_stability.cpp
    test_vtk_export.cpp
)

target_link_libraries(slope_tests PRIVATE
    slope_lib
    GTest::GTest
    GTest::Main
)

include(GoogleTest)
gtest_discover_tests(slope_tests)

#include <gtest/gtest.h>
#include "MeshGenerator.h"

class MeshGeneratorTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(MeshGeneratorTest, GenerateSimpleSlopeReturnsValidMesh) {
    auto mesh = MeshGenerator::generateSimpleSlope(10.0, 45.0, 20.0, 10.0, 10, 5);
    ASSERT_NE(mesh, nullptr);
    EXPECT_GT(mesh->numNodes(), 0);
    EXPECT_GT(mesh->numElements(), 0);
}

TEST_F(MeshGeneratorTest, SimpleSlopeHasCorrectNodeCount) {
    int nx = 10, ny = 5;
    auto mesh = MeshGenerator::generateSimpleSlope(10.0, 45.0, 20.0, 10.0, nx, ny);
    int expectedNodes = (nx + 1) * (ny + 1);
    EXPECT_EQ(mesh->numNodes(), expectedNodes);
}

TEST_F(MeshGeneratorTest, SimpleSlopeHasCorrectElementCount) {
    int nx = 10, ny = 5;
    auto mesh = MeshGenerator::generateSimpleSlope(10.0, 45.0, 20.0, 10.0, nx, ny);
    int expectedElements = nx * ny * 2;
    EXPECT_EQ(mesh->numElements(), expectedElements);
}

TEST_F(MeshGeneratorTest, SimpleSlopeHasZones) {
    auto mesh = MeshGenerator::generateSimpleSlope(10.0, 45.0, 20.0, 10.0, 10, 5);
    EXPECT_GE(mesh->getZones().size(), 2u);
}

TEST_F(MeshGeneratorTest, ZonePropertiesAreValid) {
    auto mesh = MeshGenerator::generateSimpleSlope(10.0, 45.0, 20.0, 10.0, 10, 5);
    for (const auto& [id, z] : mesh->getZones()) {
        EXPECT_GT(z.kx, 0.0);
        EXPECT_GT(z.ky, 0.0);
        EXPECT_GT(z.cohesion, 0.0);
        EXPECT_GT(z.frictionAngle, 0.0);
        EXPECT_GT(z.unitWeight, 0.0);
        EXPECT_GT(z.porosity, 0.0);
        EXPECT_LT(z.porosity, 1.0);
    }
}

TEST_F(MeshGeneratorTest, NodeCoordinatesInRange) {
    double height = 10.0, baseLen = 20.0, topLen = 10.0;
    auto mesh = MeshGenerator::generateSimpleSlope(height, 45.0, baseLen, topLen, 10, 5);
    for (const auto& n : mesh->getNodes()) {
        EXPECT_GE(n.x, 0.0);
        EXPECT_LE(n.x, baseLen + height + topLen + 1.0);
        EXPECT_GE(n.y, 0.0);
        EXPECT_LE(n.y, height + 0.1);
    }
}

TEST_F(MeshGeneratorTest, ElementNodeIdsValid) {
    auto mesh = MeshGenerator::generateSimpleSlope(10.0, 45.0, 20.0, 10.0, 10, 5);
    int maxNodeId = mesh->numNodes() - 1;
    for (const auto& e : mesh->getElements()) {
        for (int i = 0; i < 3; ++i) {
            EXPECT_GE(e.nodeIds[i], 0);
            EXPECT_LE(e.nodeIds[i], maxNodeId);
        }
    }
}

TEST_F(MeshGeneratorTest, LoadFromNonexistentFileThrows) {
    EXPECT_THROW(MeshGenerator::loadFromJSON("nonexistent_file.dat"), std::runtime_error);
}

TEST_F(MeshGeneratorTest, GetNodeCoordValid) {
    auto mesh = MeshGenerator::generateSimpleSlope(10.0, 45.0, 20.0, 10.0, 10, 5);
    Eigen::Vector2d c = mesh->getNodeCoord(0);
    EXPECT_DOUBLE_EQ(c(0), mesh->getNodes()[0].x);
    EXPECT_DOUBLE_EQ(c(1), mesh->getNodes()[0].y);
}

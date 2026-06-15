#include <gtest/gtest.h>
#include "CouplingAssembler.h"

class CouplingTest : public ::testing::Test {
protected:
    Eigen::Vector2d p0, p1, p2;

    void SetUp() override {
        p0 = Eigen::Vector2d(0.0, 0.0);
        p1 = Eigen::Vector2d(1.0, 0.0);
        p2 = Eigen::Vector2d(0.0, 1.0);
    }
};

TEST_F(CouplingTest, TriangleAreaCorrect) {
    double A = CouplingAssembler::triangleArea(p0, p1, p2);
    EXPECT_NEAR(A, 0.5, 1.0e-12);
}

TEST_F(CouplingTest, TriangleAreaDegenerate) {
    Eigen::Vector2d pd0(0.0, 0.0);
    Eigen::Vector2d pd1(1.0, 0.0);
    Eigen::Vector2d pd2(2.0, 0.0);
    double A = CouplingAssembler::triangleArea(pd0, pd1, pd2);
    EXPECT_NEAR(A, 0.0, 1.0e-12);
}

TEST_F(CouplingTest, BMatrixSize) {
    auto B = CouplingAssembler::computeBMatrix(p0, p1, p2);
    EXPECT_EQ(B.rows(), 3);
    EXPECT_EQ(B.cols(), 6);
}

TEST_F(CouplingTest, BMatrixConstantStrain) {
    auto B = CouplingAssembler::computeBMatrix(p0, p1, p2);
    Eigen::VectorXd u(6);
    u << 0.01, 0.0, 0.01, 0.0, 0.01, 0.0;
    Eigen::VectorXd eps = B * u;
    EXPECT_NEAR(eps(0), 0.0, 1.0e-10);
}

TEST_F(CouplingTest, ElementMatricesSizes) {
    auto em = CouplingAssembler::computeElementMatrices(
        p0, p1, p2, 50.0e6, 0.3, 1.0e-6, 1.0e-6, 0.35, 3600.0, 0.8, 1.0e-8);
    EXPECT_EQ(em.Kuu.rows(), 6);
    EXPECT_EQ(em.Kuu.cols(), 6);
    EXPECT_EQ(em.Kup.rows(), 6);
    EXPECT_EQ(em.Kup.cols(), 3);
    EXPECT_EQ(em.Kpu.rows(), 3);
    EXPECT_EQ(em.Kpu.cols(), 6);
    EXPECT_EQ(em.Kpp.rows(), 3);
    EXPECT_EQ(em.Kpp.cols(), 3);
}

TEST_F(CouplingTest, KuuSymmetric) {
    auto em = CouplingAssembler::computeElementMatrices(
        p0, p1, p2, 50.0e6, 0.3, 1.0e-6, 1.0e-6, 0.35, 3600.0, 0.8, 1.0e-8);
    Eigen::Matrix<double, 6, 6> diff = em.Kuu - em.Kuu.transpose();
    EXPECT_LT(diff.norm(), 1.0e-10);
}

TEST_F(CouplingTest, KppSymmetric) {
    auto em = CouplingAssembler::computeElementMatrices(
        p0, p1, p2, 50.0e6, 0.3, 1.0e-6, 1.0e-6, 0.35, 3600.0, 0.8, 1.0e-8);
    Eigen::Matrix<double, 3, 3> diff = em.Kpp - em.Kpp.transpose();
    EXPECT_LT(diff.norm(), 1.0e-10);
}

TEST_F(CouplingTest, KupKpuTranspose) {
    auto em = CouplingAssembler::computeElementMatrices(
        p0, p1, p2, 50.0e6, 0.3, 1.0e-6, 1.0e-6, 0.35, 3600.0, 0.8, 1.0e-8);
    Eigen::Matrix<double, 3, 6> diff = em.Kpu - em.Kup.transpose();
    EXPECT_LT(diff.norm(), 1.0e-10);
}

TEST_F(CouplingTest, AssemblerConstruction) {
    CouplingAssembler assembler(4);
    EXPECT_NO_THROW(assembler.buildGlobalSystem());
}

TEST_F(CouplingTest, AssembleSingleElement) {
    CouplingAssembler assembler(3);
    auto em = CouplingAssembler::computeElementMatrices(
        p0, p1, p2, 50.0e6, 0.3, 1.0e-6, 1.0e-6, 0.35, 3600.0, 0.8, 1.0e-8);
    int nodes[3] = {0, 1, 2};
    std::vector<int> fixed;
    assembler.assembleElement(em, nodes, fixed);
    EXPECT_NO_THROW(assembler.buildGlobalSystem());
}

TEST_F(CouplingTest, DisplacementsAndPressures) {
    CouplingAssembler assembler(3);
    Eigen::VectorXd disp = assembler.getDisplacements();
    Eigen::VectorXd pres = assembler.getPorePressures();
    EXPECT_EQ(disp.size(), 6);
    EXPECT_EQ(pres.size(), 3);
}

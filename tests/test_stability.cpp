#include <gtest/gtest.h>
#include "SlopeStability.h"
#include <cmath>

class StabilityTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(StabilityTest, PrincipalStressesHydrostatic) {
    double s1, s3;
    SlopeStability::computePrincipalStresses(100.0, 100.0, 0.0, s1, s3);
    EXPECT_NEAR(s1, 100.0, 1.0e-10);
    EXPECT_NEAR(s3, 100.0, 1.0e-10);
}

TEST_F(StabilityTest, PrincipalStressesPureShear) {
    double s1, s3;
    SlopeStability::computePrincipalStresses(0.0, 0.0, 50.0, s1, s3);
    EXPECT_NEAR(s1, 50.0, 1.0e-10);
    EXPECT_NEAR(s3, -50.0, 1.0e-10);
}

TEST_F(StabilityTest, PrincipalStressesGeneral) {
    double s1, s3;
    SlopeStability::computePrincipalStresses(80.0, 20.0, 30.0, s1, s3);
    EXPECT_GT(s1, s3);
    EXPECT_NEAR(s1, 100.0, 1.0e-10);
    EXPECT_NEAR(s3, 0.0, 1.0e-10);
}

TEST_F(StabilityTest, MohrCoulombFOS_Hydrostatic) {
    double fos = SlopeStability::mohrCoulombFOS(100.0, 100.0, 25.0e3, 28.0 * M_PI / 180.0);
    EXPECT_GT(fos, 1.0e6);
}

TEST_F(StabilityTest, MohrCoulumbFOS_ShearStress) {
    double phi = 28.0 * M_PI / 180.0;
    double c = 25.0e3;
    double fos = SlopeStability::mohrCoulombFOS(200.0e3, 50.0e3, c, phi);
    EXPECT_GT(fos, 0.0);
    EXPECT_LT(fos, 100.0);
}

TEST_F(StabilityTest, ComputeStabilityBasic) {
    int numElem = 4;
    Eigen::VectorXd sigmaXX(numElem), sigmaYY(numElem), sigmaXY(numElem), pp(numElem);
    sigmaXX << 100.0e3, 120.0e3, 80.0e3, 150.0e3;
    sigmaYY << 50.0e3, 60.0e3, 40.0e3, 70.0e3;
    sigmaXY << 10.0e3, 15.0e3, 20.0e3, 25.0e3;
    pp << 20.0e3, 25.0e3, 30.0e3, 35.0e3;

    std::vector<double> cohesion(numElem, 25.0e3);
    std::vector<double> friction(numElem, 28.0);
    std::vector<double> unitWeight(numElem, 18.5e3);

    std::vector<double> nodeX = {0, 1, 0, 1, 0, 1};
    std::vector<double> nodeY = {0, 0, 1, 1, 2, 2};
    std::vector<int> elemNodes = {0, 1, 2, 1, 3, 2, 2, 3, 4, 3, 5, 4};

    auto result = SlopeStability::computeStability(
        sigmaXX, sigmaYY, sigmaXY, pp,
        cohesion, friction, unitWeight,
        nodeX, nodeY, elemNodes, 3);

    EXPECT_GT(result.globalFOS, 0.0);
    EXPECT_GT(result.minFOS, 0.0);
    EXPECT_EQ(result.elementFOS.size(), static_cast<size_t>(numElem));
}

TEST_F(StabilityTest, IdentifySlipSurface) {
    std::vector<double> fos = {2.0, 1.5, 0.8, 1.2, 3.0, 0.9};
    std::vector<double> nx = {0, 1, 2, 3, 4, 5};
    std::vector<double> ny = {0, 0, 1, 1, 2, 2};
    std::vector<int> en = {0, 1, 2, 3, 4, 5};

    auto critical = SlopeStability::identifySlipSurface(fos, nx, ny, en, 1, 1.5);
    EXPECT_GE(critical.size(), 2u);
}

TEST_F(StabilityTest, GlobalFOSWeighted) {
    std::vector<double> fos = {2.0, 3.0};
    std::vector<double> areas = {1.0, 1.0};
    double gf = SlopeStability::globalFOS_SRM(fos, areas);
    EXPECT_NEAR(gf, 2.5, 1.0e-10);
}

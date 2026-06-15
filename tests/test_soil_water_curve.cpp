#include <gtest/gtest.h>
#include "SoilWaterCurve.h"
#include <cmath>

class SoilWaterCurveTest : public ::testing::Test {
protected:
    SoilWaterParams params_;

    void SetUp() override {
        params_ = SoilWaterCurve::fromZoneProperties(0.5, 1.3, 0.1, 0.35, 1.0e-6);
    }
};

TEST_F(SoilWaterCurveTest, FromZonePropertiesSetsM) {
    EXPECT_DOUBLE_EQ(params_.m, 1.0 - 1.0 / 1.3);
}

TEST_F(SoilWaterCurveTest, SaturatedDegreeOfSaturation) {
    double Se = SoilWaterCurve::degreeOfSaturation(params_, 0.0);
    EXPECT_NEAR(Se, 1.0, 1.0e-10);
}

TEST_F(SoilWaterCurveTest, UnsaturatedDegreeOfSaturation) {
    double Se = SoilWaterCurve::degreeOfSaturation(params_, 50.0e3);
    EXPECT_LT(Se, 1.0);
    EXPECT_GT(Se, params_.residualSaturation);
}

TEST_F(SoilWaterCurveTest, HighSuctionApproachesResidual) {
    double Se = SoilWaterCurve::degreeOfSaturation(params_, 1.0e8);
    EXPECT_NEAR(Se, params_.residualSaturation, 0.05);
}

TEST_F(SoilWaterCurveTest, RelativePermeabilitySaturated) {
    double kr = SoilWaterCurve::relativePermeability(params_, 0.0);
    EXPECT_NEAR(kr, 1.0, 1.0e-6);
}

TEST_F(SoilWaterCurveTest, RelativePermeabilityDecreasesWithSuction) {
    double kr0 = SoilWaterCurve::relativePermeability(params_, 10.0e3);
    double kr1 = SoilWaterCurve::relativePermeability(params_, 100.0e3);
    EXPECT_LT(kr1, kr0);
}

TEST_F(SoilWaterCurveTest, RelativePermeabilityZeroAtHighSuction) {
    double kr = SoilWaterCurve::relativePermeability(params_, 1.0e8);
    EXPECT_NEAR(kr, 0.0, 0.01);
}

TEST_F(SoilWaterCurveTest, MoistureContentConsistent) {
    double suction = 30.0e3;
    double Se = SoilWaterCurve::degreeOfSaturation(params_, suction);
    double theta = SoilWaterCurve::moistureContent(params_, suction);
    EXPECT_NEAR(theta, params_.porosity * Se, 1.0e-10);
}

TEST_F(SoilWaterCurveTest, dSdPNegativeOrZeroForPositiveSuction) {
    double dSd = SoilWaterCurve::dSdP(params_, 50.0e3);
    EXPECT_GE(dSd, 0.0);
}

TEST_F(SoilWaterCurveTest, dSdPZeroAtSaturation) {
    double dSd = SoilWaterCurve::dSdP(params_, 0.0);
    EXPECT_DOUBLE_EQ(dSd, 0.0);
}

TEST_F(SoilWaterCurveTest, ComputeElementSaturationSize) {
    int n = 5;
    Eigen::VectorXd pressure = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd elevation = Eigen::VectorXd::Constant(n, 5.0);
    Eigen::VectorXd Se = SoilWaterCurve::computeElementSaturation(params_, pressure, elevation);
    EXPECT_EQ(Se.size(), n);
}

TEST_F(SoilWaterCurveTest, ComputeElementPermeabilitySize) {
    int n = 5;
    Eigen::VectorXd pressure = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd elevation = Eigen::VectorXd::Constant(n, 5.0);
    Eigen::VectorXd Kr = SoilWaterCurve::computeElementPermeability(params_, pressure, elevation);
    EXPECT_EQ(Kr.size(), n);
}

#include <gtest/gtest.h>
#include "SeismicLoad.h"
#include "CouplingAssembler.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

static const double PI = 3.14159265358979323846;

static bool fileExists(const std::string& path) {
    std::ifstream f(path.c_str());
    return f.good();
}

static std::vector<double> simpleLinspace(double start, double end, int n) {
    std::vector<double> v(n);
    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; ++i) v[i] = start + step * i;
    return v;
}

TEST(SeismicLoadTest, HarmonicGeneration) {
    SeismicLoad s;
    double freq = 2.0;
    double amplitude = 5.0;
    int cycles = 10;
    double dt = 0.005;
    double duration = cycles / freq;
    int expectedPts = static_cast<int>(duration / dt) + 1;

    s.generateHarmonic(duration, dt, freq, amplitude, cycles, 0);

    EXPECT_EQ(s.numPoints(), expectedPts);
    EXPECT_NEAR(s.getDuration(), duration, 1e-9);
    EXPECT_GT(s.getPGA(), 0.0);
    EXPECT_LE(s.getPGA(), amplitude * 1.0 + 1e-9);

    Eigen::Vector2d accAtZero = s.getAcceleration(0.0);
    EXPECT_NEAR(accAtZero(0), 0.0, 1e-6);

    double tQuarter = 1.0 / (4.0 * freq);
    Eigen::Vector2d accPeak = s.getAcceleration(tQuarter);
    EXPECT_NEAR(std::abs(accPeak(0)), amplitude, amplitude * 0.05);
}

TEST(SeismicLoadTest, SineSweepBounds) {
    SeismicLoad s;
    double fStart = 0.5;
    double fEnd = 10.0;
    double ampl = 2.5;
    double duration = 10.0;
    double dt = 0.01;

    s.generateSineSweep(duration, dt, fStart, fEnd, ampl, 0);

    EXPECT_GT(s.numPoints(), 0);
    EXPECT_NEAR(s.getDuration(), duration, dt * 2.0);

    double measuredPGA = s.getPGA();
    EXPECT_LE(measuredPGA, ampl * 1.1 + 1e-9);
    EXPECT_GT(measuredPGA, ampl * 0.3);

    int n = s.numPoints();
    std::vector<double> times;
    for (int i = 0; i < n; ++i) times.push_back(i * dt);
    EXPECT_EQ(s.numPoints(), times.size());
}

TEST(SeismicLoadTest, ArtificialMotionRealistic) {
    SeismicLoad s;
    double magnitude = 7.0;
    double distance = 30.0;
    double duration = 20.0;
    double dt = 0.02;

    s.generateArtificialMotion(duration, dt, magnitude, distance, 0);

    EXPECT_GT(s.numPoints(), 0);
    EXPECT_NEAR(s.getDuration(), duration, dt * 2.0);

    double pga = s.getPGA();
    double predictedPGA = 10.0 * std::pow(10.0, 0.25 * magnitude - 1.3) *
                         std::exp(-0.02 * distance);
    double ratio = predictedPGA > 1e-9 ? pga / predictedPGA : 1.0;
    EXPECT_GT(ratio, 0.1);
    EXPECT_LT(ratio, 10.0);

    Eigen::Vector2d accAtStart = s.getAcceleration(0.0);
    EXPECT_NEAR(accAtStart(0), 0.0, 0.5);
    EXPECT_NEAR(accAtStart(1), 0.0, 0.5);

    Eigen::Vector2d accAtEnd = s.getAcceleration(duration - 1e-6);
    EXPECT_LT(std::abs(accAtEnd(0)), pga * 0.5 + 1e-6);
}

TEST(SeismicLoadTest, LoadTimeHistoryFile) {
    std::string testPath = "D:/trae3/a65/data/el_centro.dat";
    if (!fileExists(testPath)) {
        GTEST_SKIP() << "el_centro.dat not present, skipping file load test";
    }

    SeismicLoad s;
    bool ok = s.loadTimeHistory(testPath);
    EXPECT_TRUE(ok);
    EXPECT_GT(s.numPoints(), 50);
    EXPECT_GT(s.getDuration(), 15.0);
    EXPECT_GT(s.getPGA(), 0.5);

    Eigen::Vector2d accStart = s.getAcceleration(0.0);
    EXPECT_NEAR(accStart(0), 0.0, 0.05);

    Eigen::Vector2d accMid = s.getAcceleration(4.2);
    EXPECT_GT(std::abs(accMid(0)), 2.0);
}

TEST(SeismicLoadTest, AccelerationInterpolation) {
    SeismicLoad s;
    double freq = 5.0;
    double ampl = 10.0;
    double dt = 0.002;
    double duration = 2.0;

    s.generateHarmonic(duration, dt, freq, ampl, static_cast<int>(freq * duration), 0);

    for (int i = 0; i < 50; ++i) {
        double t = 0.01 * i + 0.001;
        double expected = ampl * std::sin(2.0 * PI * freq * t);
        Eigen::Vector2d got = s.getAcceleration(t);
        EXPECT_NEAR(got(0), expected, ampl * 0.05 + 1e-9) << " at t=" << t;
    }

    Eigen::Vector2d before = s.getAcceleration(-0.1);
    EXPECT_NEAR(before(0), 0.0, 1e-9);
    EXPECT_NEAR(before(1), 0.0, 1e-9);

    Eigen::Vector2d after = s.getAcceleration(duration + 1.0);
    EXPECT_NEAR(after(0), 0.0, 1e-9);
    EXPECT_NEAR(after(1), 0.0, 1e-9);
}

TEST(SeismicLoadTest, VelocityAndDisplacementIntegration) {
    SeismicLoad s;
    double freq = 1.0;
    double ampl = 4.0 * PI * PI;
    double dt = 0.001;
    double duration = 5.0;

    s.generateHarmonic(duration, dt, freq, ampl, static_cast<int>(freq * duration), 0);

    for (int i = 0; i < 20; ++i) {
        double t = 0.1 * i + 0.05;
        Eigen::Vector2d v = s.getVelocity(t);
        double expectedVel = -2.0 * PI * std::cos(2.0 * PI * freq * t);
        if (std::abs(expectedVel) > 1e-6) {
            double ratio = std::abs(v(0)) / std::abs(expectedVel);
            EXPECT_GT(ratio, 0.5) << " at t=" << t;
            EXPECT_LT(ratio, 1.5) << " at t=" << t;
        }
    }

    for (int i = 0; i < 10; ++i) {
        double t = 0.5 * i + 0.25;
        Eigen::Vector2d d = s.getDisplacement(t);
        double expectedDisp = std::sin(2.0 * PI * freq * t);
        if (std::abs(expectedDisp) > 1e-6) {
            double ratio = std::abs(d(0)) / std::abs(expectedDisp);
            EXPECT_GT(ratio, 0.3) << " at t=" << t;
            EXPECT_LT(ratio, 2.5) << " at t=" << t;
        }
    }
}

TEST(SeismicLoadTest, NodeLumpedMassPositive) {
    std::vector<double> x = {0.0, 1.0, 0.0, 1.0, 0.5, 0.5};
    std::vector<double> y = {0.0, 0.0, 1.0, 1.0, 0.5, 0.5};
    int nodesPerElem = 3;
    std::vector<int> conn = {
        0, 1, 4,
        1, 3, 5,
        3, 2, 5,
        2, 0, 4,
        4, 5, 2,
        4, 5, 1
    };
    std::vector<double> unitWeight = {
        20000.0, 20000.0, 20000.0, 20000.0, 20000.0, 20000.0
    };

    Eigen::VectorXd mass = SeismicLoad::computeNodeLumpedMass(x, y, conn, nodesPerElem, unitWeight, 1.0);

    EXPECT_EQ(mass.size(), 6);
    for (int i = 0; i < 6; ++i) {
        EXPECT_GT(mass(i), 0.0) << "mass(" << i << ") should be positive";
    }

    double totalArea = 0.0;
    for (int e = 0; e < 6; ++e) {
        int i0 = conn[3*e], i1 = conn[3*e+1], i2 = conn[3*e+2];
        double area = 0.5 * std::abs(
            (x[i1]-x[i0])*(y[i2]-y[i0]) - (x[i2]-x[i0])*(y[i1]-y[i0])
        );
        totalArea += area;
    }
    double expectedTotalMass = 20000.0 / 9.81 * totalArea;
    double computedTotalMass = mass.sum();
    EXPECT_NEAR(computedTotalMass, expectedTotalMass, expectedTotalMass * 0.01);
}

TEST(SeismicLoadTest, ElementMassMatrixSPD) {
    Eigen::Vector2d p0(0.0, 0.0);
    Eigen::Vector2d p1(1.0, 0.0);
    Eigen::Vector2d p2(0.0, 1.0);
    double unitWeight = 20000.0;
    double rho = unitWeight / 9.81;

    Eigen::Matrix<double, 6, 6> lumpedM = SeismicLoad::computeElementMassMatrix(p0, p1, p2, unitWeight, true);
    for (int i = 0; i < 6; ++i) {
        EXPECT_GT(lumpedM(i, i), 0.0);
        for (int j = 0; j < 6; ++j) {
            if (i != j) EXPECT_DOUBLE_EQ(lumpedM(i, j), 0.0);
        }
    }
    double totalM = 0.0;
    for (int i = 0; i < 6; ++i) totalM += lumpedM(i, i);
    double expectedM = rho * 0.5;
    EXPECT_NEAR(totalM, expectedM, expectedM * 1e-9);

    Eigen::Matrix<double, 6, 6> consistentM = SeismicLoad::computeElementMassMatrix(p0, p1, p2, unitWeight, false);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> es(consistentM);
    Eigen::VectorXd evals = es.eigenvalues();
    for (int i = 0; i < 6; ++i) {
        EXPECT_GT(evals(i), 0.0) << "Consistent mass eigenvalue " << i << " not positive";
    }
    double totalConsistent = 0.0;
    for (int i = 0; i < 6; ++i) totalConsistent += consistentM(i, i);
    EXPECT_NEAR(totalConsistent, expectedM, expectedM * 0.01);
}

TEST(SeismicLoadTest, InertialForceEquilibrium) {
    std::vector<double> x = {0.0, 2.0, 1.0};
    std::vector<double> y = {0.0, 0.0, 2.0};
    int nodesPerElem = 3;
    std::vector<int> conn = {0, 1, 2};
    std::vector<double> unitWeight = {19620.0};
    double area = 0.5 * std::abs((2.0)*(2.0) - (1.0)*(0.0));
    std::vector<double> areas = {area};

    SeismicLoad s;
    s.generateHarmonic(1.0, 0.01, 1.0, 9.81, 1, 0);

    double testTime = 0.25;
    Eigen::VectorXd fx, fy;
    s.assembleInertialForceToNodes(x, y, conn, nodesPerElem, unitWeight, areas, testTime, fx, fy, 1.0);

    EXPECT_EQ(fx.size(), 3);
    EXPECT_EQ(fy.size(), 3);

    double totalFx = fx(0) + fx(1) + fx(2);
    double totalFy = fy(0) + fy(1) + fy(2);

    double rho = 19620.0 / 9.81;
    double totalMass = rho * area;
    Eigen::Vector2d inputAcc = s.getAcceleration(testTime);
    double expectedFx = -totalMass * inputAcc(0);
    double expectedFy = -totalMass * inputAcc(1);

    EXPECT_NEAR(totalFx, expectedFx, std::abs(expectedFx) * 0.02 + 1e-9);
    EXPECT_NEAR(totalFy, expectedFy, std::abs(expectedFy) * 0.02 + 1e-9);
}

TEST(SeismicLoadTest, CouplingAssemblerSeismicInterface) {
    CouplingAssembler assembler(3);

    Eigen::VectorXd inertiaFx(3);
    Eigen::VectorXd inertiaFy(3);
    inertiaFx << 10.0, -5.0, -5.0;
    inertiaFy << 3.0, -1.5, -1.5;

    assembler.applySeismicInertia(inertiaFx, inertiaFy, 0.25, 0.5);

    std::vector<int> fixedDofs;
    Eigen::Vector2d p0(0.0, 0.0);
    Eigen::Vector2d p1(1.0, 0.0);
    Eigen::Vector2d p2(0.0, 1.0);
    ElementMatrices em = CouplingAssembler::computeDynamicElementMatrices(
        p0, p1, p2, 1e7, 0.25, 1e-6, 1e-6, 0.35, 20000.0, 1.0,
        0.5, 0.0, 1.0, true
    );
    std::vector<int> nodeIds = {0, 1, 2};
    assembler.assembleMassMatrix(em.Muu, nodeIds);
    assembler.assembleElement(em, nodeIds, fixedDofs);
    assembler.applyBoundaryConditions(fixedDofs, {});
    assembler.buildGlobalSystem();

    GlobalSystem sys = assembler.getSystem();
    EXPECT_EQ(sys.K.rows(), 9);
    EXPECT_EQ(sys.F.size(), 9);

    for (int i = 0; i < 3; ++i) {
        EXPECT_GE(sys.F(2 * i) - inertiaFx(i), -1e-9);
        EXPECT_GE(sys.F(2 * i + 1) - inertiaFy(i), -1e-9);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

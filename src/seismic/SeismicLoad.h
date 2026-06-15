#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>

struct AccelerationPoint {
    double time;
    double ax;
    double ay;
};

struct SeismicConfig {
    std::string timeHistoryFile;
    double startTime = 0.0;
    double duration = 30.0;
    double amplitudeFactor = 1.0;
    double gravity = 9.81;
    int integrationOrder = 2;
    bool applyX = true;
    bool applyY = false;
};

class SeismicLoad {
public:
    SeismicLoad() = default;

    bool loadTimeHistory(const std::string& filename);

    void generateSineSweep(
        double duration, double dt,
        double freqStart, double freqEnd,
        double amplitude, double directionDeg = 0.0);

    void generateHarmonic(
        double duration, double dt,
        double frequency, double amplitude,
        int numCycles = -1, double directionDeg = 0.0);

    void generateArtificialMotion(
        double duration, double dt,
        double magnitude, double epicenterDistance,
        double directionDeg = 0.0);

    Eigen::Vector2d getAcceleration(double time) const;

    Eigen::Vector2d getVelocity(double time) const;

    Eigen::Vector2d getDisplacement(double time) const;

    Eigen::VectorXd computeInertialForce(
        const Eigen::VectorXd& nodeMass,
        double time) const;

    void assembleInertialForceToNodes(
        const std::vector<double>& nodeX,
        const std::vector<double>& nodeY,
        const std::vector<int>& elementNodeIds,
        int nodesPerElement,
        const std::vector<double>& elementUnitWeight,
        const std::vector<double>& elementArea,
        double currentTime,
        Eigen::VectorXd& fx,
        Eigen::VectorXd& fy,
        double amplitudeScale = 1.0) const;

    Eigen::Matrix<double, 6, 6> computeElementMassMatrix(
        const Eigen::Vector2d& p0,
        const Eigen::Vector2d& p1,
        const Eigen::Vector2d& p2,
        double unitWeight,
        double lumped = true) const;

    const std::vector<AccelerationPoint>& getTimeHistory() const { return timeHistory_; }
    int numPoints() const { return static_cast<int>(timeHistory_.size()); }
    double getDuration() const;
    double getPGA() const;

    static Eigen::VectorXd computeNodeLumpedMass(
        const std::vector<double>& nodeX,
        const std::vector<double>& nodeY,
        const std::vector<int>& elementNodeIds,
        int nodesPerElement,
        const std::vector<double>& elementUnitWeight,
        double thickness = 1.0);

private:
    std::vector<AccelerationPoint> timeHistory_;
    mutable std::vector<AccelerationPoint> velocityHistory_;
    mutable std::vector<AccelerationPoint> displacementHistory_;
    mutable bool historyIntegrated_ = false;

    void integrateHistory() const;
    int findTimeIndex(double time) const;
    static double interpolateLinear(double t0, double t1, double v0, double v1, double t);
};

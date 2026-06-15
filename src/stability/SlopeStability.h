#pragma once

#include <Eigen/Dense>
#include <vector>

struct StressState {
    double sigmaXX;
    double sigmaYY;
    double sigmaXY;
    double porePressure;
};

struct StabilityResult {
    double globalFOS;
    double minFOS;
    int minFOSElement;
    std::vector<double> elementFOS;
    std::vector<int> criticalElements;
    std::vector<double> slipSurfaceY;
};

class SlopeStability {
public:
    static StabilityResult computeStability(
        const Eigen::VectorXd& sigmaXX,
        const Eigen::VectorXd& sigmaYY,
        const Eigen::VectorXd& sigmaXY,
        const Eigen::VectorXd& porePressure,
        const std::vector<double>& cohesion,
        const std::vector<double>& frictionAngle,
        const std::vector<double>& unitWeight,
        const std::vector<double>& nodeX,
        const std::vector<double>& nodeY,
        const std::vector<int>& elementNodeIds,
        int nodesPerElement);

    static double mohrCoulombFOS(
        double sigma1, double sigma3,
        double cohesion, double frictionAngleRad);

    static StressState computeElementStress(
        const Eigen::VectorXd& dispX,
        const Eigen::VectorXd& dispY,
        const Eigen::VectorXd& porePressure,
        const int nodeIds[3],
        double E, double nu_poisson);

    static double globalFOS_SRM(
        const std::vector<double>& elementFOS,
        const std::vector<double>& elementArea);

    static std::vector<int> identifySlipSurface(
        const std::vector<double>& elementFOS,
        const std::vector<double>& nodeX,
        const std::vector<double>& nodeY,
        const std::vector<int>& elementNodeIds,
        int nodesPerElement,
        double fosThreshold);

    static void computePrincipalStresses(
        double sigmaXX, double sigmaYY, double sigmaXY,
        double& sigma1, double& sigma3);

    static Eigen::VectorXd computeAllElementStress(
        const Eigen::VectorXd& dispX,
        const Eigen::VectorXd& dispY,
        const Eigen::VectorXd& porePressure,
        const std::vector<int>& elementNodeIds,
        int nodesPerElement,
        const std::vector<double>& E,
        const std::vector<double>& nu_poisson,
        const std::vector<double>& nodeX,
        const std::vector<double>& nodeY,
        Eigen::VectorXd& outSigmaXX,
        Eigen::VectorXd& outSigmaYY,
        Eigen::VectorXd& outSigmaXY);
};

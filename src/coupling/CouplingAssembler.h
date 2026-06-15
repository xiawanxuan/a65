#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <memory>

struct ElementMatrices {
    Eigen::Matrix<double, 6, 6> Kuu;
    Eigen::Matrix<double, 6, 3> Kup;
    Eigen::Matrix<double, 3, 6> Kpu;
    Eigen::Matrix<double, 3, 3> Kpp;
    Eigen::Matrix<double, 3, 3> Mpp;
    Eigen::Matrix<double, 6, 6> Muu;
    Eigen::VectorXd Fu;
    Eigen::VectorXd Fp;
};

struct GlobalSystem {
    Eigen::SparseMatrix<double> K;
    Eigen::VectorXd F;
    Eigen::VectorXd U;
    int ndof;
    int nu;
    int np;
};

class CouplingAssembler {
public:
    CouplingAssembler(int numNodes);

    void assembleElement(
        const ElementMatrices& em,
        const int nodeIds[3],
        const std::vector<int>& fixedDofs);

    void applyGravity(const std::vector<double>& nodeY,
                      const std::vector<double>& elementUnitWeight,
                      const std::vector<int[3]>& elementNodeIds);

    void applyBoundaryConditions(
        const std::vector<int>& fixedDofs,
        const std::vector<double>& fixedValues);

    void applyRainfallBC(
        const std::vector<int>& surfaceNodes,
        double rainfallRate,
        double dt,
        double currentTime);

    void applySeismicInertia(
        const Eigen::VectorXd& fx,
        const Eigen::VectorXd& fy,
        double beta = 1.0,
        double gamma = 0.5);

    void assembleMassMatrix(
        const Eigen::Matrix<double, 6, 6>& Muu,
        const int nodeIds[3]);

    void applyNewmarkAcceleration(
        const Eigen::VectorXd& acceleration,
        double beta, double gamma, double dt);

    void setPreviousState(
        const Eigen::VectorXd& prevU,
        const Eigen::VectorXd& prevV,
        const Eigen::VectorXd& prevA);

    void getNewmarkPredictor(
        Eigen::VectorXd& predictor,
        double beta, double gamma, double dt);

    void updateNewmarkState(
        const Eigen::VectorXd& deltaU,
        Eigen::VectorXd& newU,
        Eigen::VectorXd& newV,
        Eigen::VectorXd& newA,
        double beta, double gamma, double dt);

    static ElementMatrices computeDynamicElementMatrices(
        const Eigen::Vector2d& p0,
        const Eigen::Vector2d& p1,
        const Eigen::Vector2d& p2,
        double E, double nu_poisson,
        double kx, double ky,
        double porosity, double unitWeight,
        double dt,
        double Se, double dSdP_val,
        double biotAlpha = 1.0,
        double massLumping = true);

    void buildGlobalSystem();

    GlobalSystem getSystem() const { return system_; }

    Eigen::VectorXd getDisplacements() const;
    Eigen::VectorXd getPorePressures() const;

    void updateSolution(const Eigen::VectorXd& deltaU);

    static ElementMatrices computeElementMatrices(
        const Eigen::Vector2d& p0,
        const Eigen::Vector2d& p1,
        const Eigen::Vector2d& p2,
        double E, double nu_poisson,
        double kx, double ky,
        double porosity, double dt,
        double Se, double dSdP_val,
        double biotAlpha = 1.0);

    static Eigen::Matrix<double, 3, 6> computeBMatrix(
        const Eigen::Vector2d& p0,
        const Eigen::Vector2d& p1,
        const Eigen::Vector2d& p2);

    static double triangleArea(
        const Eigen::Vector2d& p0,
        const Eigen::Vector2d& p1,
        const Eigen::Vector2d& p2);

private:
    int numNodes_;
    std::vector<Eigen::Triplet<double>> triplets_;
    Eigen::VectorXd rhs_;
    Eigen::VectorXd solution_;
    GlobalSystem system_;

    int displacementDofIndex(int nodeId, int comp) const { return 2 * nodeId + comp; }
    int pressureDofIndex(int nodeId) const { return 2 * numNodes_ + nodeId; }
    int totalDofs() const { return 3 * numNodes_; }
};

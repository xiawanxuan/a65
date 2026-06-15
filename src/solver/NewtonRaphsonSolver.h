#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include <vector>
#include <string>

struct SolverConfig {
    int maxIterations = 50;
    double tolerance = 1.0e-8;
    double lineSearchMaxSteps = 5;
    bool useLineSearch = true;
    int blockSize = 0;
};

struct IterationRecord {
    int iteration;
    double residualNorm;
    double correctionNorm;
    bool converged;
    std::string message;
};

struct SolverResult {
    Eigen::VectorXd solution;
    std::vector<IterationRecord> history;
    int totalIterations;
    bool converged;
    double finalResidual;
};

class NewtonRaphsonSolver {
public:
    explicit NewtonRaphsonSolver(const SolverConfig& config = SolverConfig{});

    SolverResult solve(
        const Eigen::SparseMatrix<double>& K,
        const Eigen::VectorXd& residual);

    SolverResult solveWithCallback(
        const Eigen::SparseMatrix<double>& K,
        const Eigen::VectorXd& residual,
        std::function<void(int, const Eigen::VectorXd&)> callback);

    void setInitialGuess(const Eigen::VectorXd& guess);
    void setResidualFunction(
        std::function<Eigen::VectorXd(const Eigen::VectorXd&)> R);
    void setJacobianFunction(
        std::function<Eigen::SparseMatrix<double>(const Eigen::VectorXd&)> J);

    SolverResult solveNonlinear();

    static Eigen::VectorXd solveLinear(
        const Eigen::SparseMatrix<double>& K,
        const Eigen::VectorXd& rhs);

private:
    SolverConfig config_;
    Eigen::VectorXd currentSolution_;
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)> residualFunc_;
    std::function<Eigen::SparseMatrix<double>(const Eigen::VectorXd&)> jacobianFunc_;

    double computeResidualNorm(const Eigen::VectorXd& R) const;
    double performLineSearch(
        const Eigen::VectorXd& x,
        const Eigen::VectorXd& dx,
        const Eigen::VectorXd& R);
};

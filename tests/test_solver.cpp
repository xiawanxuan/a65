#include <gtest/gtest.h>
#include "NewtonRaphsonSolver.h"
#include <Eigen/Dense>
#include <cmath>

class SolverTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(SolverTest, LinearSolveIdentity) {
    int n = 10;
    Eigen::SparseMatrix<double> K(n, n);
    K.setIdentity();
    Eigen::VectorXd rhs = Eigen::VectorXd::Ones(n);
    Eigen::VectorXd x = NewtonRaphsonSolver::solveLinear(K, rhs);
    EXPECT_NEAR(x.norm(), std::sqrt(static_cast<double>(n)), 1.0e-8);
}

TEST_F(SolverTest, LinearSolveDiagonal) {
    int n = 5;
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < n; ++i) {
        triplets.push_back(Eigen::Triplet<double>(i, i, static_cast<double>(i + 1)));
    }
    Eigen::SparseMatrix<double> K(n, n);
    K.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::VectorXd rhs(n);
    rhs << 1.0, 2.0, 3.0, 4.0, 5.0;
    Eigen::VectorXd x = NewtonRaphsonSolver::solveLinear(K, rhs);
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(x(i), 1.0, 1.0e-8);
    }
}

TEST_F(SolverTest, SolveLinearSystem) {
    SolverConfig cfg;
    NewtonRaphsonSolver solver(cfg);

    int n = 4;
    Eigen::SparseMatrix<double> K(n, n);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.push_back(Eigen::Triplet<double>(0, 0, 4.0));
    triplets.push_back(Eigen::Triplet<double>(0, 1, 1.0));
    triplets.push_back(Eigen::Triplet<double>(1, 0, 1.0));
    triplets.push_back(Eigen::Triplet<double>(1, 1, 3.0));
    triplets.push_back(Eigen::Triplet<double>(2, 2, 5.0));
    triplets.push_back(Eigen::Triplet<double>(2, 3, 2.0));
    triplets.push_back(Eigen::Triplet<double>(3, 2, 2.0));
    triplets.push_back(Eigen::Triplet<double>(3, 3, 6.0));
    K.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd rhs(n);
    rhs << 1.0, 2.0, 3.0, 4.0;

    auto result = solver.solve(K, rhs);
    EXPECT_TRUE(result.converged);
}

TEST_F(SolverTest, NonlinearSolveQuadratic) {
    SolverConfig cfg;
    cfg.maxIterations = 100;
    cfg.tolerance = 1.0e-10;
    NewtonRaphsonSolver solver(cfg);

    Eigen::VectorXd x0 = Eigen::VectorXd::Constant(1, 3.0);
    solver.setInitialGuess(x0);

    solver.setResidualFunction([](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        Eigen::VectorXd R(1);
        R(0) = x(0) * x(0) - 4.0;
        return R;
    });

    solver.setJacobianFunction([](const Eigen::VectorXd& x) -> Eigen::SparseMatrix<double> {
        Eigen::SparseMatrix<double> J(1, 1);
        J.insert(0, 0) = 2.0 * x(0);
        J.makeCompressed();
        return J;
    });

    auto result = solver.solveNonlinear();
    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(std::abs(result.solution(0)), 2.0, 1.0e-6);
}

TEST_F(SolverTest, NonlinearSolveLinearSystem) {
    SolverConfig cfg;
    cfg.maxIterations = 20;
    cfg.tolerance = 1.0e-10;
    NewtonRaphsonSolver solver(cfg);

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(2);
    solver.setInitialGuess(x0);

    solver.setResidualFunction([](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        Eigen::VectorXd R(2);
        R(0) = 2.0 * x(0) + x(1) - 5.0;
        R(1) = x(0) + 3.0 * x(1) - 7.0;
        return R;
    });

    solver.setJacobianFunction([](const Eigen::VectorXd&) -> Eigen::SparseMatrix<double> {
        Eigen::SparseMatrix<double> J(2, 2);
        J.insert(0, 0) = 2.0;
        J.insert(0, 1) = 1.0;
        J.insert(1, 0) = 1.0;
        J.insert(1, 1) = 3.0;
        J.makeCompressed();
        return J;
    });

    auto result = solver.solveNonlinear();
    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(result.solution(0), 1.6, 1.0e-6);
    EXPECT_NEAR(result.solution(1), 1.8, 1.0e-6);
}

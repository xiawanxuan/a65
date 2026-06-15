#include "NewtonRaphsonSolver.h"
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include <iostream>
#include <algorithm>

NewtonRaphsonSolver::NewtonRaphsonSolver(const SolverConfig& config)
    : config_(config)
{
}

double NewtonRaphsonSolver::computeResidualNorm(const Eigen::VectorXd& R) const
{
    return R.norm();
}

Eigen::VectorXd NewtonRaphsonSolver::solveLinear(
    const Eigen::SparseMatrix<double>& K,
    const Eigen::VectorXd& rhs)
{
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(K);
    solver.factorize(K);

    if (solver.info() != Eigen::Success) {
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> iterSolver;
        iterSolver.setMaxIterations(1000);
        iterSolver.setTolerance(1.0e-12);
        iterSolver.compute(K);
        Eigen::VectorXd x = iterSolver.solve(rhs);
        return x;
    }

    return solver.solve(rhs);
}

SolverResult NewtonRaphsonSolver::solve(
    const Eigen::SparseMatrix<double>& K,
    const Eigen::VectorXd& residual)
{
    SolverResult result;
    result.solution = solveLinear(K, residual);
    result.converged = true;
    result.totalIterations = 1;
    result.finalResidual = residual.norm();

    IterationRecord rec;
    rec.iteration = 1;
    rec.residualNorm = result.finalResidual;
    rec.correctionNorm = result.solution.norm();
    rec.converged = true;
    rec.message = "Linear solve completed";
    result.history.push_back(rec);

    return result;
}

SolverResult NewtonRaphsonSolver::solveWithCallback(
    const Eigen::SparseMatrix<double>& K,
    const Eigen::VectorXd& residual,
    std::function<void(int, const Eigen::VectorXd&)> callback)
{
    SolverResult result = solve(K, residual);
    if (callback) callback(1, result.solution);
    return result;
}

void NewtonRaphsonSolver::setInitialGuess(const Eigen::VectorXd& guess)
{
    currentSolution_ = guess;
}

void NewtonRaphsonSolver::setResidualFunction(
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)> R)
{
    residualFunc_ = std::move(R);
}

void NewtonRaphsonSolver::setJacobianFunction(
    std::function<Eigen::SparseMatrix<double>(const Eigen::VectorXd&)> J)
{
    jacobianFunc_ = std::move(J);
}

SolverResult NewtonRaphsonSolver::solveNonlinear()
{
    SolverResult result;
    result.converged = false;
    result.totalIterations = 0;

    Eigen::VectorXd x = currentSolution_;

    for (int iter = 0; iter < config_.maxIterations; ++iter) {
        Eigen::VectorXd R = residualFunc_(x);
        double rNorm = computeResidualNorm(R);

        IterationRecord rec;
        rec.iteration = iter + 1;
        rec.residualNorm = rNorm;
        rec.converged = false;

        if (rNorm < config_.tolerance) {
            rec.converged = true;
            rec.message = "Converged: residual below tolerance";
            result.history.push_back(rec);
            result.converged = true;
            result.totalIterations = iter + 1;
            result.finalResidual = rNorm;
            result.solution = x;
            return result;
        }

        Eigen::SparseMatrix<double> J = jacobianFunc_(x);
        Eigen::VectorXd dx = solveLinear(J, -R);
        double dxNorm = dx.norm();
        rec.correctionNorm = dxNorm;

        if (config_.useLineSearch) {
            double alpha = performLineSearch(x, dx, R);
            dx *= alpha;
            rec.correctionNorm = dx.norm();
        }

        x += dx;
        rec.message = "Iteration " + std::to_string(iter + 1);
        result.history.push_back(rec);

        if (dxNorm < config_.tolerance * 1.0e-2) {
            result.converged = true;
            result.totalIterations = iter + 1;
            result.finalResidual = rNorm;
            result.solution = x;

            IterationRecord convRec;
            convRec.iteration = iter + 1;
            convRec.residualNorm = rNorm;
            convRec.correctionNorm = dxNorm;
            convRec.converged = true;
            convRec.message = "Converged: correction below tolerance";
            result.history.push_back(convRec);
            return result;
        }
    }

    result.solution = x;
    result.totalIterations = config_.maxIterations;
    result.finalResidual = computeResidualNorm(residualFunc_(x));

    IterationRecord failRec;
    failRec.iteration = config_.maxIterations;
    failRec.residualNorm = result.finalResidual;
    failRec.converged = false;
    failRec.message = "Failed to converge within max iterations";
    result.history.push_back(failRec);

    return result;
}

double NewtonRaphsonSolver::performLineSearch(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& dx,
    const Eigen::VectorXd& R)
{
    double alpha = 1.0;
    double r0 = 0.5 * R.dot(R);

    for (int k = 0; k < config_.lineSearchMaxSteps; ++k) {
        Eigen::VectorXd xNew = x + alpha * dx;
        Eigen::VectorXd Rnew = residualFunc_(xNew);
        double rNew = 0.5 * Rnew.dot(Rnew);

        if (rNew < r0 * (1.0 - 1.0e-4 * alpha)) {
            return alpha;
        }
        alpha *= 0.5;
    }

    return alpha;
}

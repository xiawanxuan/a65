#include "CouplingAssembler.h"
#include <cmath>
#include <stdexcept>

CouplingAssembler::CouplingAssembler(int numNodes)
    : numNodes_(numNodes)
{
    int ndof = totalDofs();
    rhs_.resize(ndof);
    rhs_.setZero();
    solution_.resize(ndof);
    solution_.setZero();
    system_.ndof = ndof;
    system_.nu = 2 * numNodes;
    system_.np = numNodes;
}

double CouplingAssembler::triangleArea(
    const Eigen::Vector2d& p0,
    const Eigen::Vector2d& p1,
    const Eigen::Vector2d& p2)
{
    return 0.5 * std::abs((p1(0) - p0(0)) * (p2(1) - p0(1)) - (p2(0) - p0(0)) * (p1(1) - p0(1)));
}

Eigen::Matrix<double, 3, 6> CouplingAssembler::computeBMatrix(
    const Eigen::Vector2d& p0,
    const Eigen::Vector2d& p1,
    const Eigen::Vector2d& p2)
{
    double x1 = p0(0), y1 = p0(1);
    double x2 = p1(0), y2 = p1(1);
    double x3 = p2(0), y3 = p2(1);

    double A = 2.0 * triangleArea(p0, p1, p2);
    if (A < 1.0e-30) A = 1.0e-30;

    double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
    double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;

    Eigen::Matrix<double, 3, 6> B = Eigen::Matrix<double, 3, 6>::Zero();
    B(0, 0) = b1 / A; B(0, 2) = b2 / A; B(0, 4) = b3 / A;
    B(1, 1) = c1 / A; B(1, 3) = c2 / A; B(1, 5) = c3 / A;
    B(2, 0) = c1 / A; B(2, 1) = b1 / A; B(2, 2) = c2 / A;
    B(2, 3) = b2 / A; B(2, 4) = c3 / A; B(2, 5) = b3 / A;

    return B;
}

ElementMatrices CouplingAssembler::computeElementMatrices(
    const Eigen::Vector2d& p0,
    const Eigen::Vector2d& p1,
    const Eigen::Vector2d& p2,
    double E, double nu_poisson,
    double kx, double ky,
    double porosity, double dt,
    double Se, double dSdP_val,
    double biotAlpha)
{
    ElementMatrices em;
    double A = triangleArea(p0, p1, p2);
    if (A < 1.0e-30) A = 1.0e-30;

    auto B = computeBMatrix(p0, p1, p2);

    double factor = E * (1.0 - nu_poisson) / ((1.0 + nu_poisson) * (1.0 - 2.0 * nu_poisson));
    Eigen::Matrix<double, 3, 3> D = Eigen::Matrix<double, 3, 3>::Zero();
    D(0, 0) = factor;
    D(1, 1) = factor;
    D(0, 1) = factor * nu_poisson / (1.0 - nu_poisson);
    D(1, 0) = D(0, 1);
    D(2, 2) = factor * (1.0 - 2.0 * nu_poisson) / (2.0 * (1.0 - nu_poisson));

    em.Kuu = A * B.transpose() * D * B;

    Eigen::Matrix<double, 6, 3> Np = Eigen::Matrix<double, 6, 3>::Zero();
    Np(0, 0) = 1.0 / 3.0; Np(2, 1) = 1.0 / 3.0; Np(4, 2) = 1.0 / 3.0;
    Np(1, 0) = 1.0 / 3.0; Np(3, 1) = 1.0 / 3.0; Np(5, 2) = 1.0 / 3.0;

    Eigen::Matrix<double, 3, 3> Np_block = (1.0 / 3.0) * Eigen::Matrix<double, 3, 3>::Ones();

    em.Kup = -biotAlpha * A * B.transpose() * Eigen::Vector3d(1, 1, 0).asDiagonal() * Np_block;
    em.Kpu = em.Kup.transpose();

    double x1 = p0(0), y1 = p0(1);
    double x2 = p1(0), y2 = p1(1);
    double x3 = p2(0), y3 = p2(1);
    double A2 = 2.0 * A;
    double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
    double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;

    Eigen::Matrix<double, 3, 3> gradOp = Eigen::Matrix<double, 3, 3>::Zero();
    gradOp(0, 0) = b1; gradOp(0, 1) = b2; gradOp(0, 2) = b3;
    gradOp(1, 0) = c1; gradOp(1, 1) = c2; gradOp(1, 2) = c3;
    gradOp /= A2;

    Eigen::Matrix<double, 3, 3> permTensor = Eigen::Matrix<double, 3, 3>::Zero();
    permTensor(0, 0) = kx;
    permTensor(1, 1) = ky;

    em.Kpp = -dt * A * gradOp.transpose() * permTensor * gradOp;

    double storage = porosity * dSdP_val / (9.81e3);
    em.Mpp = storage * A * Np_block;
    em.Kpp += em.Mpp;

    em.Fu = Eigen::VectorXd::Zero(6);
    em.Fu(1) = 0.0;
    em.Fu(3) = 0.0;
    em.Fu(5) = 0.0;

    em.Fp = Eigen::VectorXd::Zero(3);

    return em;
}

void CouplingAssembler::assembleElement(
    const ElementMatrices& em,
    const int nodeIds[3],
    const std::vector<int>& fixedDofs)
{
    std::vector<int> dofIndices(9);
    for (int i = 0; i < 3; ++i) {
        dofIndices[2 * i] = displacementDofIndex(nodeIds[i], 0);
        dofIndices[2 * i + 1] = displacementDofIndex(nodeIds[i], 1);
    }
    for (int i = 0; i < 3; ++i) {
        dofIndices[6 + i] = pressureDofIndex(nodeIds[i]);
    }

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            triplets_.push_back(Eigen::Triplet<double>(dofIndices[i], dofIndices[j], em.Kuu(i, j)));
        }
        rhs_(dofIndices[i]) += em.Fu(i);
    }

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 3; ++j) {
            triplets_.push_back(Eigen::Triplet<double>(dofIndices[i], dofIndices[6 + j], em.Kup(i, j)));
            triplets_.push_back(Eigen::Triplet<double>(dofIndices[6 + j], dofIndices[i], em.Kpu(j, i)));
        }
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            triplets_.push_back(Eigen::Triplet<double>(dofIndices[6 + i], dofIndices[6 + j], em.Kpp(i, j)));
        }
        rhs_(dofIndices[6 + i]) += em.Fp(i);
    }
}

void CouplingAssembler::applyGravity(
    const std::vector<double>& nodeY,
    const std::vector<double>& elementUnitWeight,
    const std::vector<int[3]>& elementNodeIds)
{
    for (size_t e = 0; e < elementUnitWeight.size(); ++e) {
        double gamma = elementUnitWeight[e];
        double bodyForceY = -gamma / 3.0;
        for (int i = 0; i < 3; ++i) {
            int nid = elementNodeIds[e][i];
            int ydof = displacementDofIndex(nid, 1);
            rhs_(ydof) += bodyForceY;
        }
    }
}

void CouplingAssembler::applyBoundaryConditions(
    const std::vector<int>& fixedDofs,
    const std::vector<double>& fixedValues)
{
    for (size_t i = 0; i < fixedDofs.size(); ++i) {
        int dof = fixedDofs[i];
        double val = fixedValues[i];

        triplets_.push_back(Eigen::Triplet<double>(dof, dof, 1.0e30));
        rhs_(dof) = val * 1.0e30;
    }
}

void CouplingAssembler::applyRainfallBC(
    const std::vector<int>& surfaceNodes,
    double rainfallRate,
    double dt,
    double currentTime)
{
    double fluxIntensity = rainfallRate * 9.81e3;
    for (int nid : surfaceNodes) {
        int pdof = pressureDofIndex(nid);
        rhs_(pdof) += fluxIntensity * dt;
    }
}

void CouplingAssembler::buildGlobalSystem()
{
    int ndof = totalDofs();
    system_.K.resize(ndof, ndof);
    system_.K.setFromTriplets(triplets_.begin(), triplets_.end());
    system_.F = rhs_;
    system_.U = solution_;
    system_.ndof = ndof;
}

Eigen::VectorXd CouplingAssembler::getDisplacements() const
{
    Eigen::VectorXd u(2 * numNodes_);
    for (int i = 0; i < 2 * numNodes_; ++i) {
        u(i) = solution_(i);
    }
    return u;
}

Eigen::VectorXd CouplingAssembler::getPorePressures() const
{
    Eigen::VectorXd p(numNodes_);
    for (int i = 0; i < numNodes_; ++i) {
        p(i) = solution_(2 * numNodes_ + i);
    }
    return p;
}

void CouplingAssembler::updateSolution(const Eigen::VectorXd& deltaU)
{
    solution_ += deltaU;
}

void CouplingAssembler::applySeismicInertia(
    const Eigen::VectorXd& fx,
    const Eigen::VectorXd& fy,
    double beta, double gamma)
{
    for (int i = 0; i < numNodes_; ++i) {
        if (i < fx.size()) {
            int xdof = displacementDofIndex(i, 0);
            rhs_(xdof) += fx(i);
        }
        if (i < fy.size()) {
            int ydof = displacementDofIndex(i, 1);
            rhs_(ydof) += fy(i);
        }
    }
}

void CouplingAssembler::assembleMassMatrix(
    const Eigen::Matrix<double, 6, 6>& Muu,
    const int nodeIds[3])
{
    std::vector<int> dofIndices(6);
    for (int i = 0; i < 3; ++i) {
        dofIndices[2 * i] = displacementDofIndex(nodeIds[i], 0);
        dofIndices[2 * i + 1] = displacementDofIndex(nodeIds[i], 1);
    }

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (std::abs(Muu(i, j)) > 1.0e-30) {
                triplets_.push_back(Eigen::Triplet<double>(dofIndices[i], dofIndices[j], Muu(i, j)));
            }
        }
    }
}

void CouplingAssembler::applyNewmarkAcceleration(
    const Eigen::VectorXd& acceleration,
    double beta, double gamma, double dt)
{
    for (int i = 0; i < numNodes_; ++i) {
        if (2 * i + 1 >= acceleration.size()) continue;
        double ax = acceleration(2 * i);
        double ay = acceleration(2 * i + 1);

        double massCoeff = 1.0 / (beta * dt * dt);
        double xdof = displacementDofIndex(i, 0);
        double ydof = displacementDofIndex(i, 1);

        triplets_.push_back(Eigen::Triplet<double>(xdof, xdof, massCoeff));
        triplets_.push_back(Eigen::Triplet<double>(ydof, ydof, massCoeff));
    }
}

void CouplingAssembler::setPreviousState(
    const Eigen::VectorXd& prevU,
    const Eigen::VectorXd& prevV,
    const Eigen::VectorXd& prevA)
{
    if (prevU.size() >= 2 * numNodes_) {
        for (int i = 0; i < 2 * numNodes_; ++i) {
            solution_(i) = prevU(i);
        }
    }
}

void CouplingAssembler::getNewmarkPredictor(
    Eigen::VectorXd& predictor,
    double beta, double gamma, double dt)
{
    predictor.resize(2 * numNodes_);
    predictor.setZero();

    Eigen::VectorXd uPrev = solution_.head(2 * numNodes_);
    predictor = uPrev;
}

void CouplingAssembler::updateNewmarkState(
    const Eigen::VectorXd& deltaU,
    Eigen::VectorXd& newU,
    Eigen::VectorXd& newV,
    Eigen::VectorXd& newA,
    double beta, double gamma, double dt)
{
    int nu = 2 * numNodes_;
    newU = solution_.head(nu) + deltaU.head(nu);

    if (newV.size() != nu) newV.resize(nu);
    if (newA.size() != nu) newA.resize(nu);

    double a1 = 1.0 / (beta * dt * dt);
    double a2 = gamma / (beta * dt);
    double a3 = 1.0 / (beta * dt);
    double a4 = gamma / beta - 1.0;

    for (int i = 0; i < nu; ++i) {
        newA(i) = a1 * deltaU(i) - a3 * newV(i) - a4 * newA(i);
        newV(i) += gamma * dt * newA(i);
    }

    solution_.head(nu) = newU;
}

ElementMatrices CouplingAssembler::computeDynamicElementMatrices(
    const Eigen::Vector2d& p0,
    const Eigen::Vector2d& p1,
    const Eigen::Vector2d& p2,
    double E, double nu_poisson,
    double kx, double ky,
    double porosity, double unitWeight,
    double dt,
    double Se, double dSdP_val,
    double biotAlpha,
    double massLumping)
{
    ElementMatrices em = computeElementMatrices(
        p0, p1, p2, E, nu_poisson, kx, ky, porosity, dt, Se, dSdP_val, biotAlpha);

    double A = triangleArea(p0, p1, p2);
    if (A < 1.0e-30) A = 1.0e-30;
    double rho = unitWeight / 9.81;
    double nodeMass = rho * A / 3.0;

    em.Muu.setZero();
    if (massLumping) {
        for (int i = 0; i < 3; ++i) {
            em.Muu(2 * i, 2 * i) = nodeMass;
            em.Muu(2 * i + 1, 2 * i + 1) = nodeMass;
        }
    } else {
        double c12 = rho * A / 12.0;
        double c11 = rho * A / 6.0;
        for (int i = 0; i < 3; ++i) {
            em.Muu(2 * i, 2 * i) = c11;
            em.Muu(2 * i + 1, 2 * i + 1) = c11;
            for (int j = i + 1; j < 3; ++j) {
                em.Muu(2 * i, 2 * j) = c12;
                em.Muu(2 * j, 2 * i) = c12;
                em.Muu(2 * i + 1, 2 * j + 1) = c12;
                em.Muu(2 * j + 1, 2 * i + 1) = c12;
            }
        }
    }

    return em;
}

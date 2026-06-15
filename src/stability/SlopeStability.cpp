#include "SlopeStability.h"
#include <cmath>
#include <algorithm>
#include <numeric>

void SlopeStability::computePrincipalStresses(
    double sigmaXX, double sigmaYY, double sigmaXY,
    double& sigma1, double& sigma3)
{
    double avg = 0.5 * (sigmaXX + sigmaYY);
    double diff = 0.5 * (sigmaXX - sigmaYY);
    double R = std::sqrt(diff * diff + sigmaXY * sigmaXY);
    sigma1 = avg + R;
    sigma3 = avg - R;
}

double SlopeStability::mohrCoulombFOS(
    double sigma1, double sigma3,
    double cohesion, double frictionAngleRad)
{
    if (sigma1 < sigma3) std::swap(sigma1, sigma3);

    double phi = frictionAngleRad;
    double shearStrength = cohesion * std::cos(phi) + 0.5 * (sigma1 + sigma3) * std::sin(phi);
    double shearStress = 0.5 * (sigma1 - sigma3);

    if (shearStress < 1.0e-10) return 1.0e10;
    return shearStrength / shearStress;
}

StabilityResult SlopeStability::computeStability(
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
    int nodesPerElement)
{
    StabilityResult result;
    int numElements = static_cast<int>(sigmaXX.size());
    result.elementFOS.resize(numElements);
    result.minFOS = 1.0e10;
    result.minFOSElement = 0;

    std::vector<double> areas(numElements, 0.0);

    for (int e = 0; e < numElements; ++e) {
        double s1, s3;
        double effectiveSigmaXX = sigmaXX(e) + porePressure(e);
        double effectiveSigmaYY = sigmaYY(e) + porePressure(e);
        computePrincipalStresses(effectiveSigmaXX, effectiveSigmaYY, sigmaXY(e), s1, s3);

        double phiRad = frictionAngle[e] * M_PI / 180.0;
        double fos = mohrCoulombFOS(s1, s3, cohesion[e], phiRad);
        result.elementFOS[e] = fos;

        if (fos < result.minFOS) {
            result.minFOS = fos;
            result.minFOSElement = e;
        }

        if (nodesPerElement == 3) {
            int n0 = elementNodeIds[e * nodesPerElement];
            int n1 = elementNodeIds[e * nodesPerElement + 1];
            int n2 = elementNodeIds[e * nodesPerElement + 2];
            double x0 = nodeX[n0], y0 = nodeY[n0];
            double x1 = nodeX[n1], y1 = nodeY[n1];
            double x2 = nodeX[n2], y2 = nodeY[n2];
            areas[e] = 0.5 * std::abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
        }
    }

    result.globalFOS = globalFOS_SRM(result.elementFOS, areas);
    result.criticalElements = identifySlipSurface(
        result.elementFOS, nodeX, nodeY, elementNodeIds, nodesPerElement,
        result.globalFOS * 1.2);

    return result;
}

double SlopeStability::globalFOS_SRM(
    const std::vector<double>& elementFOS,
    const std::vector<double>& elementArea)
{
    double weightedSum = 0.0;
    double totalArea = 0.0;

    for (size_t e = 0; e < elementFOS.size(); ++e) {
        double w = elementArea[e];
        double fos = std::min(elementFOS[e], 10.0);
        weightedSum += fos * w;
        totalArea += w;
    }

    if (totalArea < 1.0e-30) return 0.0;
    return weightedSum / totalArea;
}

std::vector<int> SlopeStability::identifySlipSurface(
    const std::vector<double>& elementFOS,
    const std::vector<double>& nodeX,
    const std::vector<double>& nodeY,
    const std::vector<int>& elementNodeIds,
    int nodesPerElement,
    double fosThreshold)
{
    std::vector<int> critical;
    for (size_t e = 0; e < elementFOS.size(); ++e) {
        if (elementFOS[e] < fosThreshold) {
            critical.push_back(static_cast<int>(e));
        }
    }
    return critical;
}

StressState SlopeStability::computeElementStress(
    const Eigen::VectorXd& dispX,
    const Eigen::VectorXd& dispY,
    const Eigen::VectorXd& porePressure,
    const int nodeIds[3],
    double E, double nu_poisson)
{
    StressState state{};
    return state;
}

Eigen::VectorXd SlopeStability::computeAllElementStress(
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
    Eigen::VectorXd& outSigmaXY)
{
    int numElements = static_cast<int>(elementNodeIds.size()) / nodesPerElement;
    outSigmaXX.resize(numElements);
    outSigmaYY.resize(numElements);
    outSigmaXY.resize(numElements);
    Eigen::VectorXd devStress(numElements);

    for (int e = 0; e < numElements; ++e) {
        int n0 = elementNodeIds[e * nodesPerElement];
        int n1 = elementNodeIds[e * nodesPerElement + 1];
        int n2 = elementNodeIds[e * nodesPerElement + 2];

        double x0 = nodeX[n0], y0 = nodeY[n0];
        double x1 = nodeX[n1], y1 = nodeY[n1];
        double x2 = nodeX[n2], y2 = nodeY[n2];

        double A2 = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
        if (std::abs(A2) < 1.0e-30) A2 = 1.0e-30;

        double b1 = y1 - y2, b2 = y2 - y0, b3 = y0 - y1;
        double c1 = x2 - x1, c2 = x0 - x2, c3 = x1 - x0;

        double ux0 = dispX(n0), uy0 = dispY(n0);
        double ux1 = dispX(n1), uy1 = dispY(n1);
        double ux2 = dispX(n2), uy2 = dispY(n2);

        double epsXX = (b1 * ux0 + b2 * ux1 + b3 * ux2) / A2;
        double epsYY = (c1 * uy0 + c2 * uy1 + c3 * uy2) / A2;
        double gammaXY = ((c1 * ux0 + c2 * ux1 + c3 * ux2) + (b1 * uy0 + b2 * uy1 + b3 * uy2)) / A2;

        double nu = nu_poisson[e];
        double eMod = E[e];
        double factor = eMod / ((1.0 + nu) * (1.0 - 2.0 * nu));

        outSigmaXX(e) = factor * ((1.0 - nu) * epsXX + nu * epsYY);
        outSigmaYY(e) = factor * (nu * epsXX + (1.0 - nu) * epsYY);
        outSigmaXY(e) = factor * (1.0 - 2.0 * nu) / 2.0 * gammaXY;

        devStress(e) = std::sqrt(outSigmaXX(e) * outSigmaXX(e) +
                                 outSigmaYY(e) * outSigmaYY(e) +
                                 3.0 * outSigmaXY(e) * outSigmaXY(e));
    }

    return devStress;
}

#include "SoilWaterCurve.h"
#include <cmath>
#include <algorithm>

SoilWaterParams SoilWaterCurve::fromZoneProperties(
    double alpha, double nVG, double sRes, double porosity, double ksat)
{
    SoilWaterParams params;
    params.alpha = alpha;
    params.n = nVG;
    params.m = 1.0 - 1.0 / nVG;
    params.residualSaturation = sRes;
    params.porosity = porosity;
    params.ksat = ksat;
    return params;
}

double SoilWaterCurve::suctionFromPressure(double porePressure, double elevation) {
    double totalHead = porePressure / 9.81e3 + elevation;
    return std::max(0.0, -totalHead * 9.81e3);
}

double SoilWaterCurve::effectiveSaturation(const SoilWaterParams& params, double suction) {
    if (suction <= 0.0) return 1.0;
    double h = suction / (params.alpha * 9.81e3);
    double se = std::pow(1.0 + std::pow(h, params.n), -params.m);
    return std::clamp(se, 0.0, 1.0);
}

double SoilWaterCurve::degreeOfSaturation(const SoilWaterParams& params, double suction) {
    double se = effectiveSaturation(params, suction);
    return params.residualSaturation + (1.0 - params.residualSaturation) * se;
}

double SoilWaterCurve::dSdP(const SoilWaterParams& params, double suction) {
    if (suction <= 0.0) return 0.0;
    double h = suction / (params.alpha * 9.81e3);
    double hn = std::pow(h, params.n);
    double base = 1.0 + hn;
    double dSe = -params.m * params.n * std::pow(base, -params.m - 1.0) * hn / h;
    double dSuction_dP = -1.0 / (params.alpha * 9.81e3);
    return (1.0 - params.residualSaturation) * dSe * dSuction_dP;
}

double SoilWaterCurve::relativePermeability(const SoilWaterParams& params, double suction) {
    double se = effectiveSaturation(params, suction);
    if (se <= 0.0) return 0.0;
    if (se >= 1.0) return 1.0;
    double kr = std::sqrt(se) * std::pow(1.0 - std::pow(1.0 - std::pow(se, 1.0 / params.m), params.m), 2.0);
    return std::clamp(kr, 0.0, 1.0);
}

double SoilWaterCurve::dKrdS(const SoilWaterParams& params, double suction) {
    double eps = 1.0;
    double kr0 = relativePermeability(params, suction);
    double kr1 = relativePermeability(params, suction + eps);
    return (kr1 - kr0) / eps;
}

double SoilWaterCurve::moistureContent(const SoilWaterParams& params, double suction) {
    return params.porosity * degreeOfSaturation(params, suction);
}

Eigen::VectorXd SoilWaterCurve::computeElementSaturation(
    const SoilWaterParams& params,
    const Eigen::VectorXd& porePressure,
    const Eigen::VectorXd& elevation)
{
    int n = static_cast<int>(porePressure.size());
    Eigen::VectorXd Se(n);
    for (int i = 0; i < n; ++i) {
        double suction = suctionFromPressure(porePressure(i), elevation(i));
        Se(i) = degreeOfSaturation(params, suction);
    }
    return Se;
}

Eigen::VectorXd SoilWaterCurve::computeElementPermeability(
    const SoilWaterParams& params,
    const Eigen::VectorXd& porePressure,
    const Eigen::VectorXd& elevation)
{
    int n = static_cast<int>(porePressure.size());
    Eigen::VectorXd Kr(n);
    for (int i = 0; i < n; ++i) {
        double suction = suctionFromPressure(porePressure(i), elevation(i));
        Kr(i) = relativePermeability(params, suction);
    }
    return Kr;
}

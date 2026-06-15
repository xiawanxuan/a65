#pragma once

#include <Eigen/Dense>
#include <vector>

struct SoilWaterParams {
    double alpha;
    double n;
    double m;
    double residualSaturation;
    double porosity;
    double ksat;
};

class SoilWaterCurve {
public:
    static SoilWaterParams fromZoneProperties(
        double alpha, double nVG, double sRes, double porosity, double ksat);

    static double degreeOfSaturation(const SoilWaterParams& params, double suction);

    static double dSdP(const SoilWaterParams& params, double suction);

    static double relativePermeability(const SoilWaterParams& params, double suction);

    static double dKrdS(const SoilWaterParams& params, double suction);

    static double moistureContent(const SoilWaterParams& params, double suction);

    static Eigen::VectorXd computeElementSaturation(
        const SoilWaterParams& params,
        const Eigen::VectorXd& porePressure,
        const Eigen::VectorXd& elevation);

    static Eigen::VectorXd computeElementPermeability(
        const SoilWaterParams& params,
        const Eigen::VectorXd& porePressure,
        const Eigen::VectorXd& elevation);

private:
    static double effectiveSaturation(const SoilWaterParams& params, double suction);
    static double suctionFromPressure(double porePressure, double elevation);
};

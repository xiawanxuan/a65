#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

struct VTKFieldData {
    std::string name;
    int numComponents;
    std::vector<double> data;
};

class VTKExporter {
public:
    static void exportUnstructuredGrid(
        const std::string& filename,
        const std::vector<double>& nodeX,
        const std::vector<double>& nodeY,
        const std::vector<int>& elementNodeIds,
        int nodesPerElement,
        const std::vector<VTKFieldData>& pointData,
        const std::vector<VTKFieldData>& cellData = {});

    static void exportTimeStep(
        const std::string& outputDir,
        int timeStep,
        double currentTime,
        const std::vector<double>& nodeX,
        const std::vector<double>& nodeY,
        const std::vector<int>& elementNodeIds,
        int nodesPerElement,
        const Eigen::VectorXd& displacementX,
        const Eigen::VectorXd& displacementY,
        const Eigen::VectorXd& porePressure,
        const Eigen::VectorXd& saturation,
        const Eigen::VectorXd& effectiveStressXX,
        const Eigen::VectorXd& effectiveStressYY,
        const Eigen::VectorXd& effectiveStressXY,
        const Eigen::VectorXd& safetyFactor);

    static void writePVDFile(
        const std::string& filename,
        const std::vector<std::pair<double, std::string>>& timeStepFiles);

private:
    static std::string formatFloat(double val, int precision = 15);
    static std::string getVTKCellType(int nodesPerElement);
    static void serializeFieldData(
        std::ofstream& ofs,
        const VTKFieldData& field,
        const std::string& indent);
};

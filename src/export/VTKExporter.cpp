#include "VTKExporter.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <cmath>
#include <limits>
#include <algorithm>

namespace fs = std::filesystem;

std::string VTKExporter::formatFloat(double val, int precision) {
    if (std::isnan(val)) {
        return "0.000000000000000e+00";
    }
    if (std::isinf(val)) {
        if (val > 0.0) {
            return "1.000000000000000e+100";
        } else {
            return "-1.000000000000000e+100";
        }
    }

    if (val == 0.0) {
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(precision) << 0.0;
        std::string s = oss.str();
        std::string::size_type ePos = s.find('e');
        if (ePos == std::string::npos) ePos = s.find('E');
        if (ePos != std::string::npos) {
            s.replace(ePos, 1, "e");
        }
        return s;
    }

    std::ostringstream oss;
    oss << std::scientific << std::uppercase << std::setprecision(precision) << val;
    std::string s = oss.str();

    std::string::size_type ePos = s.find('E');
    if (ePos != std::string::npos) {
        s.replace(ePos, 1, "e");

        std::string expStr = s.substr(ePos + 1);
        if (!expStr.empty()) {
            char expSign = '+';
            size_t numStart = 0;
            if (expStr[0] == '+' || expStr[0] == '-') {
                expSign = expStr[0];
                numStart = 1;
            }
            std::string expNum = expStr.substr(numStart);
            if (expNum.size() < 2) {
                expNum = std::string(2 - expNum.size(), '0') + expNum;
            }
            s = s.substr(0, ePos + 1) + expSign + expNum;
        }
    }

    return s;
}

std::string VTKExporter::getVTKCellType(int nodesPerElement) {
    if (nodesPerElement == 3) return "5";
    if (nodesPerElement == 4) return "9";
    return "7";
}

void VTKExporter::serializeFieldData(
    std::ofstream& ofs,
    const VTKFieldData& field,
    const std::string& indent)
{
    int stride = field.numComponents;
    size_t totalSize = field.data.size();
    size_t numTuples = totalSize / stride;

    const int valuesPerLine = 3;
    int lineCount = 0;

    for (size_t i = 0; i < numTuples; ++i) {
        if (lineCount == 0) {
            ofs << indent;
        }

        for (int c = 0; c < stride; ++c) {
            size_t idx = i * stride + c;
            if (idx < totalSize) {
                double val = field.data[idx];
                ofs << formatFloat(val) << " ";
            }
        }

        lineCount++;
        if (lineCount >= valuesPerLine) {
            ofs << "\n";
            lineCount = 0;
        }
    }

    if (lineCount > 0) {
        ofs << "\n";
    }
}

void VTKExporter::exportUnstructuredGrid(
    const std::string& filename,
    const std::vector<double>& nodeX,
    const std::vector<double>& nodeY,
    const std::vector<int>& elementNodeIds,
    int nodesPerElement,
    const std::vector<VTKFieldData>& pointData,
    const std::vector<VTKFieldData>& cellData)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open()) return;

    int numNodes = static_cast<int>(nodeX.size());
    int numElements = static_cast<int>(elementNodeIds.size()) / nodesPerElement;

    ofs << "<?xml version=\"1.0\"?>\n";
    ofs << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    ofs << "  <UnstructuredGrid>\n";
    ofs << "    <Piece NumberOfPoints=\"" << numNodes << "\" NumberOfCells=\"" << numElements << "\">\n";

    ofs << "      <Points>\n";
    ofs << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < numNodes; ++i) {
        ofs << "          " << formatFloat(nodeX[i]) << " "
            << formatFloat(nodeY[i]) << " 0.0\n";
    }
    ofs << "        </DataArray>\n";
    ofs << "      </Points>\n";

    if (!pointData.empty()) {
        ofs << "      <PointData>\n";
        for (const auto& pd : pointData) {
            ofs << "        <DataArray type=\"Float64\" Name=\"" << pd.name
                << "\" NumberOfComponents=\"" << pd.numComponents
                << "\" format=\"ascii\">\n";
            serializeFieldData(ofs, pd, "          ");
            ofs << "        </DataArray>\n";
        }
        ofs << "      </PointData>\n";
    }

    if (!cellData.empty()) {
        ofs << "      <CellData>\n";
        for (const auto& cd : cellData) {
            ofs << "        <DataArray type=\"Float64\" Name=\"" << cd.name
                << "\" NumberOfComponents=\"" << cd.numComponents
                << "\" format=\"ascii\">\n";
            serializeFieldData(ofs, cd, "          ");
            ofs << "        </DataArray>\n";
        }
        ofs << "      </CellData>\n";
    }

    ofs << "      <Cells>\n";
    ofs << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    ofs << "          ";
    for (int e = 0; e < numElements; ++e) {
        for (int n = 0; n < nodesPerElement; ++n) {
            ofs << elementNodeIds[e * nodesPerElement + n] << " ";
        }
    }
    ofs << "\n        </DataArray>\n";

    ofs << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    ofs << "          ";
    for (int e = 1; e <= numElements; ++e) {
        ofs << e * nodesPerElement << " ";
    }
    ofs << "\n        </DataArray>\n";

    ofs << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    ofs << "          ";
    std::string cellType = getVTKCellType(nodesPerElement);
    for (int e = 0; e < numElements; ++e) {
        ofs << cellType << " ";
    }
    ofs << "\n        </DataArray>\n";
    ofs << "      </Cells>\n";

    ofs << "    </Piece>\n";
    ofs << "  </UnstructuredGrid>\n";
    ofs << "</VTKFile>\n";

    ofs.close();
}

void VTKExporter::exportTimeStep(
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
    const Eigen::VectorXd& safetyFactor)
{
    fs::create_directories(outputDir);

    std::ostringstream fname;
    fname << outputDir << "/slope_step_"
          << std::setfill('0') << std::setw(6) << timeStep << ".vtu";

    int numNodes = static_cast<int>(nodeX.size());
    int numElements = static_cast<int>(elementNodeIds.size()) / nodesPerElement;

    VTKFieldData dispField;
    dispField.name = "Displacement";
    dispField.numComponents = 3;
    dispField.data.resize(numNodes * 3, 0.0);
    for (int i = 0; i < numNodes; ++i) {
        dispField.data[i * 3] = displacementX(i);
        dispField.data[i * 3 + 1] = displacementY(i);
    }

    VTKFieldData pressureField;
    pressureField.name = "PorePressure";
    pressureField.numComponents = 1;
    pressureField.data.resize(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        pressureField.data[i] = porePressure(i);
    }

    VTKFieldData satField;
    satField.name = "Saturation";
    satField.numComponents = 1;
    satField.data.resize(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        satField.data[i] = saturation(i);
    }

    VTKFieldData stressField;
    stressField.name = "EffectiveStress";
    stressField.numComponents = 3;
    stressField.data.resize(numElements * 3, 0.0);
    for (int i = 0; i < numElements; ++i) {
        stressField.data[i * 3] = effectiveStressXX(i);
        stressField.data[i * 3 + 1] = effectiveStressYY(i);
        stressField.data[i * 3 + 2] = effectiveStressXY(i);
    }

    VTKFieldData fosField;
    fosField.name = "SafetyFactor";
    fosField.numComponents = 1;
    fosField.data.resize(numElements);
    for (int i = 0; i < numElements; ++i) {
        fosField.data[i] = safetyFactor(i);
    }

    exportUnstructuredGrid(
        fname.str(), nodeX, nodeY, elementNodeIds, nodesPerElement,
        {dispField, pressureField, satField},
        {stressField, fosField});
}

void VTKExporter::writePVDFile(
    const std::string& filename,
    const std::vector<std::pair<double, std::string>>& timeStepFiles)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open()) return;

    ofs << "<?xml version=\"1.0\"?>\n";
    ofs << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
    ofs << "  <Collection>\n";
    for (const auto& [time, file] : timeStepFiles) {
        ofs << "    <DataSet timestep=\"" << formatFloat(time)
            << "\" file=\"" << file << "\"/>\n";
    }
    ofs << "  </Collection>\n";
    ofs << "</VTKFile>\n";
    ofs.close();
}

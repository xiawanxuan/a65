#include <gtest/gtest.h>
#include "VTKExporter.h"
#include <fstream>
#include <filesystem>
#include <sstream>
#include <cmath>

class VTKExportTest : public ::testing::Test {
protected:
    std::string testDir;

    void SetUp() override {
        testDir = "test_output_vtk";
        std::filesystem::create_directories(testDir);
    }

    void TearDown() override {
        std::filesystem::remove_all(testDir);
    }

    static int countOccurrences(const std::string& str, const std::string& sub) {
        int count = 0;
        size_t pos = 0;
        while ((pos = str.find(sub, pos)) != std::string::npos) {
            count++;
            pos += sub.length();
        }
        return count;
    }
};

TEST_F(VTKExportTest, ExportUnstructuredGridCreatesFile) {
    std::string filename = testDir + "/test.vtu";
    std::vector<double> nodeX = {0.0, 1.0, 0.0, 1.0};
    std::vector<double> nodeY = {0.0, 0.0, 1.0, 1.0};
    std::vector<int> elemNodes = {0, 1, 2, 1, 3, 2};

    VTKExporter::exportUnstructuredGrid(filename, nodeX, nodeY, elemNodes, 3, {});

    EXPECT_TRUE(std::filesystem::exists(filename));
}

TEST_F(VTKExportTest, ExportUnstructuredGridWithPointData) {
    std::string filename = testDir + "/test_pd.vtu";
    std::vector<double> nodeX = {0.0, 1.0, 0.0};
    std::vector<double> nodeY = {0.0, 0.0, 1.0};
    std::vector<int> elemNodes = {0, 1, 2};

    VTKFieldData pd;
    pd.name = "Pressure";
    pd.numComponents = 1;
    pd.data = {100.0, 200.0, 150.0};

    VTKExporter::exportUnstructuredGrid(filename, nodeX, nodeY, elemNodes, 3, {pd});

    EXPECT_TRUE(std::filesystem::exists(filename));
    std::ifstream ifs(filename);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("Pressure"), std::string::npos);
}

TEST_F(VTKExportTest, ExportUnstructuredGridWithCellData) {
    std::string filename = testDir + "/test_cd.vtu";
    std::vector<double> nodeX = {0.0, 1.0, 0.0, 1.0};
    std::vector<double> nodeY = {0.0, 0.0, 1.0, 1.0};
    std::vector<int> elemNodes = {0, 1, 2, 1, 3, 2};

    VTKFieldData cd;
    cd.name = "SafetyFactor";
    cd.numComponents = 1;
    cd.data = {1.5, 2.0};

    VTKExporter::exportUnstructuredGrid(filename, nodeX, nodeY, elemNodes, 3, {}, {cd});

    EXPECT_TRUE(std::filesystem::exists(filename));
    std::ifstream ifs(filename);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("SafetyFactor"), std::string::npos);
}

TEST_F(VTKExportTest, ExportTimeStepCreatesVTU) {
    int numNodes = 3;
    int numElem = 1;
    std::vector<double> nodeX = {0.0, 1.0, 0.0};
    std::vector<double> nodeY = {0.0, 0.0, 1.0};
    std::vector<int> elemNodes = {0, 1, 2};

    Eigen::VectorXd dx = Eigen::VectorXd::Zero(numNodes);
    Eigen::VectorXd dy = Eigen::VectorXd::Zero(numNodes);
    Eigen::VectorXd pp = Eigen::VectorXd::Zero(numNodes);
    Eigen::VectorXd sat = Eigen::VectorXd::Ones(numNodes);
    Eigen::VectorXd sxx = Eigen::VectorXd::Zero(numElem);
    Eigen::VectorXd syy = Eigen::VectorXd::Zero(numElem);
    Eigen::VectorXd sxy = Eigen::VectorXd::Zero(numElem);
    Eigen::VectorXd fos = Eigen::VectorXd::Ones(numElem);

    VTKExporter::exportTimeStep(testDir, 0, 0.0, nodeX, nodeY, elemNodes, 3,
                                 dx, dy, pp, sat, sxx, syy, sxy, fos);

    EXPECT_TRUE(std::filesystem::exists(testDir + "/slope_step_000000.vtu"));
}

TEST_F(VTKExportTest, WritePVDFile) {
    std::string pvdFile = testDir + "/test.pvd";
    std::vector<std::pair<double, std::string>> entries = {
        {0.0, "step_0.vtu"},
        {3600.0, "step_1.vtu"}
    };

    VTKExporter::writePVDFile(pvdFile, entries);

    EXPECT_TRUE(std::filesystem::exists(pvdFile));
    std::ifstream ifs(pvdFile);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("Collection"), std::string::npos);
    EXPECT_NE(content.find("step_0.vtu"), std::string::npos);
}

TEST_F(VTKExportTest, LowPorePressureNotTruncated) {
    std::string filename = testDir + "/low_pp.vtu";
    std::vector<double> nodeX = {0.0, 1.0, 0.0};
    std::vector<double> nodeY = {0.0, 0.0, 1.0};
    std::vector<int> elemNodes = {0, 1, 2};

    VTKFieldData pd;
    pd.name = "PorePressure";
    pd.numComponents = 1;
    pd.data = {
        -1.0e-3,
        -1.23456789012345e-7,
        -1.0e-10
    };

    VTKExporter::exportUnstructuredGrid(filename, nodeX, nodeY, elemNodes, 3, {pd});

    std::ifstream ifs(filename);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());

    EXPECT_NE(content.find("-1.0"), std::string::npos);
    EXPECT_NE(content.find("e-03"), std::string::npos);
    EXPECT_NE(content.find("e-07"), std::string::npos);
    EXPECT_NE(content.find("e-10"), std::string::npos);
}

TEST_F(VTKExportTest, NegativePorePressurePreserved) {
    std::string filename = testDir + "/neg_pp.vtu";
    std::vector<double> nodeX = {0.0, 1.0, 0.0};
    std::vector<double> nodeY = {0.0, 0.0, 1.0};
    std::vector<int> elemNodes = {0, 1, 2};

    VTKFieldData pd;
    pd.name = "PorePressure";
    pd.numComponents = 1;
    pd.data = {-100.0e3, -50.0, -0.001};

    VTKExporter::exportUnstructuredGrid(filename, nodeX, nodeY, elemNodes, 3, {pd});

    std::ifstream ifs(filename);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());

    int minusCount = countOccurrences(content, "-");
    EXPECT_GE(minusCount, 3);
    EXPECT_NE(content.find("-1.00"), std::string::npos);
    EXPECT_NE(content.find("-5.00"), std::string::npos);
    EXPECT_NE(content.find("-1.00"), std::string::npos);
}

TEST_F(VTKExportTest, FieldDataAllValuesWritten) {
    std::string filename = testDir + "/all_values.vtu";
    int numNodes = 50;
    std::vector<double> nodeX(numNodes);
    std::vector<double> nodeY(numNodes);
    std::vector<int> elemNodes;

    for (int i = 0; i < numNodes; ++i) {
        nodeX[i] = static_cast<double>(i);
        nodeY[i] = 0.0;
    }
    for (int i = 0; i < numNodes - 1; ++i) {
        elemNodes.push_back(i);
        elemNodes.push_back(i + 1);
        elemNodes.push_back(i);
    }

    VTKFieldData pd;
    pd.name = "PorePressure";
    pd.numComponents = 1;
    pd.data.resize(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        pd.data[i] = -std::pow(0.5, i) * 100.0;
    }

    VTKExporter::exportUnstructuredGrid(filename, nodeX, nodeY, elemNodes, 3, {pd});

    std::ifstream ifs(filename);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());

    int eCount = countOccurrences(content, "e-");
    EXPECT_GE(eCount, 10);
    EXPECT_NE(content.find("PorePressure"), std::string::npos);
}

TEST_F(VTKExportTest, MultiComponentFieldComplete) {
    std::string filename = testDir + "/multi_comp.vtu";
    std::vector<double> nodeX = {0.0, 1.0, 0.0, 1.0};
    std::vector<double> nodeY = {0.0, 0.0, 1.0, 1.0};
    std::vector<int> elemNodes = {0, 1, 2, 1, 3, 2};

    VTKFieldData pd;
    pd.name = "Displacement";
    pd.numComponents = 3;
    pd.data = {
        0.001, -0.002, 0.0,
        0.003, -0.004, 0.0,
        0.005, -0.006, 0.0,
        0.007, -0.008, 0.0
    };

    VTKExporter::exportUnstructuredGrid(filename, nodeX, nodeY, elemNodes, 3, {pd});

    std::ifstream ifs(filename);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());

    int negativeCount = countOccurrences(content, "-");
    EXPECT_GE(negativeCount, 4);
}

TEST_F(VTKExportTest, VerySmallValuesSerialized) {
    std::string filename = testDir + "/tiny_values.vtu";
    std::vector<double> nodeX = {0.0, 1.0, 0.0};
    std::vector<double> nodeY = {0.0, 0.0, 1.0};
    std::vector<int> elemNodes = {0, 1, 2};

    VTKFieldData pd;
    pd.name = "TinyValues";
    pd.numComponents = 1;
    pd.data = {1.234567890123456e-15, 0.0, -1.234567890123456e-15};

    VTKExporter::exportUnstructuredGrid(filename, nodeX, nodeY, elemNodes, 3, {pd});

    std::ifstream ifs(filename);
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());

    EXPECT_NE(content.find("e-15"), std::string::npos);
    EXPECT_NE(content.find("1.23456"), std::string::npos);
}

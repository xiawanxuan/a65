#include <gtest/gtest.h>
#include "VTKExporter.h"
#include <fstream>
#include <filesystem>

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

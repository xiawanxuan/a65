#pragma once

#include <string>
#include <vector>

struct SimConfig {
    std::string meshFile;
    std::string outputDir = "output";
    double totalTime = 86400.0;
    double timeStep = 3600.0;
    double rainfallRate = 1.0e-6;
    double rainfallStartTime = 0.0;
    double rainfallEndTime = 43200.0;
    int maxNewtonIterations = 50;
    double newtonTolerance = 1.0e-8;
    bool useLineSearch = true;
    int blockSize = 0;
    bool generateMesh = false;
    double slopeHeight = 10.0;
    double slopeAngle = 45.0;
    double baseLength = 20.0;
    double topLength = 10.0;
    int meshNX = 20;
    int meshNY = 10;
    bool verbose = false;
    bool exportIntermediate = false;
};

class CLIParser {
public:
    static SimConfig parse(int argc, char* argv[]);
    static void printHelp();
    static bool validate(const SimConfig& config);
};

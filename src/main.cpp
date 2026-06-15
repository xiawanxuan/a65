#include "MeshGenerator.h"
#include "SoilWaterCurve.h"
#include "CouplingAssembler.h"
#include "NewtonRaphsonSolver.h"
#include "VTKExporter.h"
#include "SlopeStability.h"
#include "CLIParser.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <filesystem>

namespace fs = std::filesystem;

static std::vector<int> findSurfaceNodes(const Mesh& mesh) {
    std::vector<int> surfaceNodes;
    double maxY = -1.0e30;
    for (const auto& n : mesh.getNodes()) {
        if (n.y > maxY) maxY = n.y;
    }

    double tol = 1.0e-6;
    for (const auto& n : mesh.getNodes()) {
        if (std::abs(n.y - maxY) < tol) {
            surfaceNodes.push_back(n.id);
        }
    }

    double minY = 1.0e30;
    for (const auto& n : mesh.getNodes()) {
        if (n.y < minY) minY = n.y;
    }
    double midY = 0.5 * (minY + maxY);
    for (const auto& n : mesh.getNodes()) {
        if (std::abs(n.y - midY) < 0.05 * (maxY - minY)) {
            bool already = false;
            for (int sid : surfaceNodes) {
                if (sid == n.id) { already = true; break; }
            }
            if (!already) surfaceNodes.push_back(n.id);
        }
    }

    return surfaceNodes;
}

static std::vector<int> findFixedDofs(const Mesh& mesh) {
    std::vector<int> fixedDofs;
    double minY = 1.0e30;
    double minX = 1.0e30;
    double maxX = -1.0e30;

    for (const auto& n : mesh.getNodes()) {
        if (n.y < minY) minY = n.y;
        if (n.x < minX) minX = n.x;
        if (n.x > maxX) maxX = n.x;
    }

    double tol = 1.0e-6;
    int numNodes = mesh.numNodes();

    for (const auto& n : mesh.getNodes()) {
        if (std::abs(n.y - minY) < tol) {
            fixedDofs.push_back(2 * n.id);
            fixedDofs.push_back(2 * n.id + 1);
        }
        if (std::abs(n.x - minX) < tol) {
            fixedDofs.push_back(2 * n.id);
        }
        if (std::abs(n.x - maxX) < tol) {
            fixedDofs.push_back(2 * n.id);
        }
    }

    return fixedDofs;
}

int main(int argc, char* argv[])
{
    SimConfig config = CLIParser::parse(argc, argv);
    if (!CLIParser::validate(config)) return 1;

    std::cout << "========================================\n";
    std::cout << "  Slope Stability Simulation\n";
    std::cout << "  Rainfall-Induced Seepage-Stress Coupling\n";
    std::cout << "========================================\n\n";

    std::unique_ptr<Mesh> mesh;
    if (config.generateMesh) {
        std::cout << "[Mesh] Generating simple slope mesh...\n";
        mesh = MeshGenerator::generateSimpleSlope(
            config.slopeHeight, config.slopeAngle,
            config.baseLength, config.topLength,
            config.meshNX, config.meshNY);
    } else {
        std::cout << "[Mesh] Loading mesh from: " << config.meshFile << "\n";
        mesh = MeshGenerator::loadFromJSON(config.meshFile);
    }

    int numNodes = mesh->numNodes();
    int numElements = mesh->numElements();
    std::cout << "[Mesh] Nodes: " << numNodes
              << ", Elements: " << numElements << "\n\n";

    std::vector<double> nodeX(numNodes), nodeY(numNodes);
    Eigen::VectorXd elevation(numNodes);
    for (const auto& n : mesh->getNodes()) {
        nodeX[n.id] = n.x;
        nodeY[n.id] = n.y;
        elevation(n.id) = n.y;
    }

    std::vector<int> flatElementNodes;
    for (const auto& e : mesh->getElements()) {
        flatElementNodes.push_back(e.nodeIds[0]);
        flatElementNodes.push_back(e.nodeIds[1]);
        flatElementNodes.push_back(e.nodeIds[2]);
    }

    std::vector<double> elemE(numElements), elemNu(numElements);
    std::vector<double> elemCohesion(numElements), elemFriction(numElements);
    std::vector<double> elemUnitWeight(numElements), elemKx(numElements), elemKy(numElements);
    std::vector<double> elemPorosity(numElements), elemAlpha(numElements), elemN(numElements), elemSres(numElements);

    for (const auto& e : mesh->getElements()) {
        const auto& z = mesh->getZoneProperties(e.zoneId);
        elemE[e.id] = 50.0e6;
        elemNu[e.id] = 0.3;
        elemCohesion[e.id] = z.cohesion;
        elemFriction[e.id] = z.frictionAngle;
        elemUnitWeight[e.id] = z.unitWeight;
        elemKx[e.id] = z.kx;
        elemKy[e.id] = z.ky;
        elemPorosity[e.id] = z.porosity;
        elemAlpha[e.id] = z.vanGenuchtenAlpha;
        elemN[e.id] = z.vanGenuchtenN;
        elemSres[e.id] = z.residualSaturation;
    }

    std::vector<int> surfaceNodes = findSurfaceNodes(*mesh);
    std::vector<int> fixedDofs = findFixedDofs(*mesh);
    std::vector<double> fixedValues(fixedDofs.size(), 0.0);

    std::vector<int[3]> elemNodeArray(numElements);
    std::vector<double> elemGamma(numElements);
    for (const auto& e : mesh->getElements()) {
        elemNodeArray[e.id][0] = e.nodeIds[0];
        elemNodeArray[e.id][1] = e.nodeIds[1];
        elemNodeArray[e.id][2] = e.nodeIds[2];
        elemGamma[e.id] = elemUnitWeight[e.id];
    }

    fs::create_directories(config.outputDir);

    Eigen::VectorXd porePressure = Eigen::VectorXd::Zero(numNodes);
    for (const auto& n : mesh->getNodes()) {
        double depth = config.slopeHeight - n.y;
        if (depth > 0) {
            porePressure(n.id) = depth * 9.81e3 * 0.5;
        }
    }

    SolverConfig solverCfg;
    solverCfg.maxIterations = config.maxNewtonIterations;
    solverCfg.tolerance = config.newtonTolerance;
    solverCfg.useLineSearch = config.useLineSearch;
    solverCfg.blockSize = config.blockSize;
    NewtonRaphsonSolver solver(solverCfg);

    std::vector<std::pair<double, std::string>> pvdEntries;
    int numTimeSteps = static_cast<int>(std::ceil(config.totalTime / config.timeStep));

    std::cout << "[Sim] Time steps: " << numTimeSteps
              << ", dt = " << config.timeStep << " s\n";
    std::cout << "[Sim] Rainfall: " << config.rainfallRate << " m/s from t="
              << config.rainfallStartTime << " to t=" << config.rainfallEndTime << "\n\n";

    for (int step = 0; step < numTimeSteps; ++step) {
        double currentTime = step * config.timeStep;
        double nextTime = (step + 1) * config.timeStep;
        double dt = config.timeStep;

        std::cout << "--- Time step " << step + 1 << "/" << numTimeSteps
                  << " (t = " << std::fixed << std::setprecision(1) << currentTime
                  << " -> " << nextTime << " s) ---\n";

        CouplingAssembler assembler(numNodes);

        for (const auto& e : mesh->getElements()) {
            int n0 = e.nodeIds[0], n1 = e.nodeIds[1], n2 = e.nodeIds[2];
            Eigen::Vector2d p0(nodeX[n0], nodeY[n0]);
            Eigen::Vector2d p1(nodeX[n1], nodeY[n1]);
            Eigen::Vector2d p2(nodeX[n2], nodeY[n2]);

            const auto& z = mesh->getZoneProperties(e.zoneId);
            auto swParams = SoilWaterCurve::fromZoneProperties(
                z.vanGenuchtenAlpha, z.vanGenuchtenN,
                z.residualSaturation, z.porosity, z.kx);

            double avgPressure = (porePressure(n0) + porePressure(n1) + porePressure(n2)) / 3.0;
            double avgElev = (nodeY[n0] + nodeY[n1] + nodeY[n2]) / 3.0;
            double suction = SoilWaterCurve::suctionFromPressure(avgPressure, avgElev);
            double Se = SoilWaterCurve::degreeOfSaturation(swParams, suction);
            double dSdP_val = SoilWaterCurve::dSdP(swParams, suction);

            auto em = CouplingAssembler::computeElementMatrices(
                p0, p1, p2,
                elemE[e.id], elemNu[e.id],
                z.kx * SoilWaterCurve::relativePermeability(swParams, suction),
                z.ky * SoilWaterCurve::relativePermeability(swParams, suction),
                z.porosity, dt, Se, dSdP_val);

            assembler.assembleElement(em, e.nodeIds, fixedDofs);
        }

        assembler.applyGravity(nodeY, elemGamma, elemNodeArray);
        assembler.applyBoundaryConditions(fixedDofs, fixedValues);

        bool isRaining = (currentTime >= config.rainfallStartTime &&
                          currentTime < config.rainfallEndTime);
        if (isRaining) {
            assembler.applyRainfallBC(surfaceNodes, config.rainfallRate, dt, currentTime);
            if (config.verbose) {
                std::cout << "  [Rain] Applying rainfall: " << config.rainfallRate << " m/s\n";
            }
        }

        assembler.buildGlobalSystem();

        GlobalSystem sys = assembler.getSystem();
        Eigen::VectorXd residual = sys.F - sys.K * sys.U;

        auto result = solver.solve(sys.K, residual);
        assembler.updateSolution(result.solution);

        porePressure = assembler.getPorePressures();

        Eigen::VectorXd dispX = assembler.getDisplacements().head(numNodes);
        Eigen::VectorXd dispYVec(numNodes);
        for (int i = 0; i < numNodes; ++i) {
            dispYVec(i) = assembler.getDisplacements()(2 * i + 1);
        }
        Eigen::VectorXd dispXVec(numNodes);
        for (int i = 0; i < numNodes; ++i) {
            dispXVec(i) = assembler.getDisplacements()(2 * i);
        }

        Eigen::VectorXd sigmaXX, sigmaYY, sigmaXY;
        SlopeStability::computeAllElementStress(
            dispXVec, dispYVec, porePressure,
            flatElementNodes, 3, elemE, elemNu,
            nodeX, nodeY, sigmaXX, sigmaYY, sigmaXY);

        auto stability = SlopeStability::computeStability(
            sigmaXX, sigmaYY, sigmaXY, porePressure,
            elemCohesion, elemFriction, elemUnitWeight,
            nodeX, nodeY, flatElementNodes, 3);

        std::cout << "  Global FOS: " << std::setprecision(4) << stability.globalFOS
                  << ", Min FOS: " << stability.minFOS
                  << " (elem " << stability.minFOSElement << ")\n";

        if (!stability.criticalElements.empty()) {
            std::cout << "  Critical elements: " << stability.criticalElements.size() << "\n";
        }

        Eigen::VectorXd saturation(numNodes);
        for (int i = 0; i < numNodes; ++i) {
            double suction_i = SoilWaterCurve::suctionFromPressure(porePressure(i), elevation(i));
            auto dummyParams = SoilWaterCurve::fromZoneProperties(0.5, 1.3, 0.1, 0.35, 1.0e-6);
            saturation(i) = SoilWaterCurve::degreeOfSaturation(dummyParams, suction_i);
        }

        Eigen::VectorXd elemFOS = Eigen::Map<Eigen::VectorXd>(
            stability.elementFOS.data(), stability.elementFOS.size());

        VTKExporter::exportTimeStep(
            config.outputDir, step, currentTime,
            nodeX, nodeY, flatElementNodes, 3,
            dispXVec, dispYVec, porePressure, saturation,
            sigmaXX, sigmaYY, sigmaXY, elemFOS);

        std::ostringstream vtuName;
        vtuName << "slope_step_" << std::setfill('0') << std::setw(6) << step << ".vtu";
        pvdEntries.emplace_back(currentTime, vtuName.str());

        if (config.exportIntermediate && result.history.size() > 1) {
            for (const auto& rec : result.history) {
                std::cout << "    Iter " << rec.iteration
                          << ": |R|=" << std::scientific << rec.residualNorm
                          << ", |du|=" << rec.correctionNorm;
                if (rec.converged) std::cout << " [CONVERGED]";
                std::cout << "\n";
            }
        }
    }

    std::string pvdPath = config.outputDir + "/slope_simulation.pvd";
    VTKExporter::writePVDFile(pvdPath, pvdEntries);

    std::cout << "\n========================================\n";
    std::cout << "  Simulation Complete\n";
    std::cout << "  Output: " << config.outputDir << "\n";
    std::cout << "  PVD file: " << pvdPath << "\n";
    std::cout << "========================================\n";

    return 0;
}

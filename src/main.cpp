#include "MeshGenerator.h"
#include "SoilWaterCurve.h"
#include "CouplingAssembler.h"
#include "NewtonRaphsonSolver.h"
#include "VTKExporter.h"
#include "SlopeStability.h"
#include "CLIParser.h"
#include "SeismicLoad.h"

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
    std::cout << "  Rainfall-Seismic Seepage-Stress Coupling\n";
    if (config.enableSeismic) std::cout << "  [Seismic Coupling ENABLED]\n";
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
    std::vector<double> elemArea(numElements);

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

        int n0 = e.nodeIds[0], n1 = e.nodeIds[1], n2 = e.nodeIds[2];
        double x0 = nodeX[n0], y0 = nodeY[n0];
        double x1 = nodeX[n1], y1 = nodeY[n1];
        double x2 = nodeX[n2], y2 = nodeY[n2];
        elemArea[e.id] = 0.5 * std::abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
    }

    SeismicLoad seismic;
    Eigen::VectorXd nodeLumpedMass(numNodes);
    nodeLumpedMass.setZero();

    if (config.enableSeismic) {
        std::cout << "[Seismic] Initializing seismic loading module...\n";
        double accelScale = config.seismicAmplitude * 9.81;
        double seismicDt = std::min(config.timeStep, 0.02);

        switch (config.seismicMode) {
        case 0:
            std::cout << "[Seismic] Mode: Harmonic excitation, f="
                      << config.seismicHarmonicFreq << " Hz\n";
            seismic.generateHarmonic(config.seismicDuration, seismicDt,
                config.seismicHarmonicFreq, accelScale,
                config.seismicHarmonicCycles, config.seismicDirection);
            break;
        case 1:
            std::cout << "[Seismic] Mode: Sine sweep, " << config.seismicFreqStart
                      << "-" << config.seismicFreqEnd << " Hz\n";
            seismic.generateSineSweep(config.seismicDuration, seismicDt,
                config.seismicFreqStart, config.seismicFreqEnd,
                accelScale, config.seismicDirection);
            break;
        case 2:
            std::cout << "[Seismic] Mode: Artificial motion, M="
                      << config.seismicMagnitude << ", R=" << config.seismicDistance << " km\n";
            seismic.generateArtificialMotion(config.seismicDuration, seismicDt,
                config.seismicMagnitude, config.seismicDistance, config.seismicDirection);
            break;
        case 3:
            std::cout << "[Seismic] Mode: Loading external file: " << config.seismicFile << "\n";
            if (!seismic.loadTimeHistory(config.seismicFile)) {
                std::cerr << "[Seismic] Failed to load seismic file!\n";
                return 2;
            }
            break;
        default:
            std::cout << "[Seismic] Mode: Harmonic (default)\n";
            seismic.generateHarmonic(config.seismicDuration, seismicDt,
                config.seismicHarmonicFreq, accelScale,
                config.seismicHarmonicCycles, config.seismicDirection);
        }

        std::cout << "[Seismic] PGA: " << seismic.getPGA()
                  << " m/s^2, Duration: " << seismic.getDuration() << " s, Points: "
                  << seismic.numPoints() << "\n";

        nodeLumpedMass = SeismicLoad::computeNodeLumpedMass(
            nodeX, nodeY, flatElementNodes, 3, elemUnitWeight, 1.0);

        double totalMass = nodeLumpedMass.sum();
        std::cout << "[Seismic] Total lumped mass: " << totalMass << " kg\n\n";
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

    Eigen::VectorXd velocityX = Eigen::VectorXd::Zero(numNodes);
    Eigen::VectorXd velocityY = Eigen::VectorXd::Zero(numNodes);
    Eigen::VectorXd accelerationX = Eigen::VectorXd::Zero(numNodes);
    Eigen::VectorXd accelerationY = Eigen::VectorXd::Zero(numNodes);
    Eigen::VectorXd dispXPrev_ = Eigen::VectorXd::Zero(numNodes);
    Eigen::VectorXd dispYPrev_ = Eigen::VectorXd::Zero(numNodes);

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
              << config.rainfallStartTime << " to t=" << config.rainfallEndTime << " s\n";
    if (config.enableSeismic) {
        std::cout << "[Sim] Seismic: from t=" << config.seismicStartTime
                  << " to t=" << config.seismicStartTime + config.seismicDuration
                  << " s (mode=" << config.seismicMode << ")\n";
    }
    std::cout << "\n";

    for (int step = 0; step < numTimeSteps; ++step) {
        double currentTime = step * config.timeStep;
        double nextTime = (step + 1) * config.timeStep;
        double dt = config.timeStep;
        double midTime = 0.5 * (currentTime + nextTime);

        std::cout << "--- Time step " << step + 1 << "/" << numTimeSteps
                  << " (t = " << std::fixed << std::setprecision(2) << currentTime
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

            ElementMatrices em;
            if (config.enableSeismic) {
                em = CouplingAssembler::computeDynamicElementMatrices(
                    p0, p1, p2,
                    elemE[e.id], elemNu[e.id],
                    z.kx * SoilWaterCurve::relativePermeability(swParams, suction),
                    z.ky * SoilWaterCurve::relativePermeability(swParams, suction),
                    z.porosity, z.unitWeight,
                    dt, Se, dSdP_val, 1.0, true);
            } else {
                em = CouplingAssembler::computeElementMatrices(
                    p0, p1, p2,
                    elemE[e.id], elemNu[e.id],
                    z.kx * SoilWaterCurve::relativePermeability(swParams, suction),
                    z.ky * SoilWaterCurve::relativePermeability(swParams, suction),
                    z.porosity, dt, Se, dSdP_val);
            }

            assembler.assembleElement(em, e.nodeIds, fixedDofs);
        }

        assembler.applyGravity(nodeY, elemGamma, elemNodeArray);

        bool isSeismicActive = false;
        Eigen::VectorXd inertiaFx = Eigen::VectorXd::Zero(numNodes);
        Eigen::VectorXd inertiaFy = Eigen::VectorXd::Zero(numNodes);
        Eigen::VectorXd inputAx = Eigen::VectorXd::Zero(numNodes);
        Eigen::VectorXd inputAy = Eigen::VectorXd::Zero(numNodes);

        if (config.enableSeismic &&
            midTime >= config.seismicStartTime &&
            midTime < config.seismicStartTime + config.seismicDuration)
        {
            isSeismicActive = true;
            double relTime = midTime - config.seismicStartTime;
            Eigen::Vector2d inputAcc = seismic.getAcceleration(relTime);

            if (!config.seismicApplyX) inputAcc(0) = 0.0;
            if (!config.seismicApplyY) inputAcc(1) = 0.0;

            for (int i = 0; i < numNodes; ++i) {
                inputAx(i) = inputAcc(0);
                inputAy(i) = inputAcc(1);
            }

            seismic.assembleInertialForceToNodes(
                nodeX, nodeY, flatElementNodes, 3,
                elemUnitWeight, elemArea, relTime,
                inertiaFx, inertiaFy, config.seismicAmplitude);

            assembler.applySeismicInertia(inertiaFx, inertiaFy,
                config.newmarkBeta, config.newmarkGamma);

            if (config.verbose) {
                std::cout << "  [Seismic] Input acc: (" << inputAcc(0) << ", "
                          << inputAcc(1) << ") m/s^2\n";
            }
        }

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

        Eigen::VectorXd dispXVec = assembler.getDisplacements().head(numNodes);
        Eigen::VectorXd dispYVec(numNodes);
        for (int i = 0; i < numNodes; ++i) {
            dispXVec(i) = assembler.getDisplacements()(2 * i);
            dispYVec(i) = assembler.getDisplacements()(2 * i + 1);
        }

        if (config.enableSeismic) {
            double safeDt = std::max(dt, 1.0e-12);
            Eigen::VectorXd newVelX(numNodes), newVelY(numNodes);
            Eigen::VectorXd newAccX(numNodes), newAccY(numNodes);
            if (step > 0) {
                for (int i = 0; i < numNodes; ++i) {
                    newVelX(i) = (dispXVec(i) - dispXPrev_(i)) / safeDt;
                    newVelY(i) = (dispYVec(i) - dispYPrev_(i)) / safeDt;
                    newAccX(i) = (newVelX(i) - velocityX(i)) / safeDt;
                    newAccY(i) = (newVelY(i) - velocityY(i)) / safeDt;
                }
            } else {
                for (int i = 0; i < numNodes; ++i) {
                    newVelX(i) = dispXVec(i) / safeDt;
                    newVelY(i) = dispYVec(i) / safeDt;
                    newAccX(i) = newVelX(i) / safeDt;
                    newAccY(i) = newVelY(i) / safeDt;
                }
            }
            velocityX = newVelX;
            velocityY = newVelY;
            accelerationX = newAccX;
            accelerationY = newAccY;
            dispXPrev_ = dispXVec;
            dispYPrev_ = dispYVec;
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

        if (isSeismicActive) {
            std::cout << "  Seismic: Input="
                      << std::fixed << std::setprecision(3) << inputAx(0) << " m/s^2"
                      << " | FOS drop factor: " << std::setprecision(3)
                      << (stability.minFOS < 1.5 ? "WARNING-LOW" : "SAFE") << "\n";
        }

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

        if (config.enableSeismic) {
            VTKExporter::exportTimeStepSeismic(
                config.outputDir, step, currentTime,
                nodeX, nodeY, flatElementNodes, 3,
                dispXVec, dispYVec, porePressure, saturation,
                sigmaXX, sigmaYY, sigmaXY, elemFOS,
                accelerationX, accelerationY,
                velocityX, velocityY,
                inertiaFx, inertiaFy,
                nodeLumpedMass,
                inputAx, inputAy,
                true);
        } else {
            VTKExporter::exportTimeStep(
                config.outputDir, step, currentTime,
                nodeX, nodeY, flatElementNodes, 3,
                dispXVec, dispYVec, porePressure, saturation,
                sigmaXX, sigmaYY, sigmaXY, elemFOS);
        }

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
    if (config.enableSeismic) {
        std::cout << "  Seismic fields: Acceleration, Velocity, InertiaForce, NodeMass, InputAcceleration\n";
    }
    std::cout << "========================================\n";

    return 0;
}

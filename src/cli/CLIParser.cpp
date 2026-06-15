#include "CLIParser.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>

namespace po = boost::program_options;

SimConfig CLIParser::parse(int argc, char* argv[])
{
    SimConfig config;

    po::options_description general("General Options");
    general.add_options()
        ("help,h", "Show help message")
        ("mesh,m", po::value<std::string>(&config.meshFile), "Input mesh file path")
        ("output,o", po::value<std::string>(&config.outputDir)->default_value("output"),
         "Output directory for VTK files")
        ("verbose,v", po::bool_switch(&config.verbose), "Enable verbose output");

    po::options_description time("Time Options");
    time.add_options()
        ("total-time,T", po::value<double>(&config.totalTime)->default_value(86400.0),
         "Total simulation time in seconds")
        ("time-step,dt", po::value<double>(&config.timeStep)->default_value(3600.0),
         "Time step size in seconds")
        ("rainfall-rate", po::value<double>(&config.rainfallRate)->default_value(1.0e-6),
         "Rainfall intensity in m/s")
        ("rainfall-start", po::value<double>(&config.rainfallStartTime)->default_value(0.0),
         "Rainfall start time in seconds")
        ("rainfall-end", po::value<double>(&config.rainfallEndTime)->default_value(43200.0),
         "Rainfall end time in seconds");

    po::options_description solver("Solver Options");
    solver.add_options()
        ("max-iter", po::value<int>(&config.maxNewtonIterations)->default_value(50),
         "Maximum Newton-Raphson iterations per step")
        ("tolerance", po::value<double>(&config.newtonTolerance)->default_value(1.0e-8),
         "Newton-Raphson convergence tolerance")
        ("line-search", po::bool_switch(&config.useLineSearch)->default_value(true),
         "Enable line search in Newton iterations")
        ("block-size", po::value<int>(&config.blockSize)->default_value(0),
         "Block size for chunked loading (0=load all)");

    po::options_description meshgen("Mesh Generation Options");
    meshgen.add_options()
        ("generate-mesh", po::bool_switch(&config.generateMesh),
         "Generate simple slope mesh instead of loading")
        ("slope-height", po::value<double>(&config.slopeHeight)->default_value(10.0),
         "Slope height for mesh generation")
        ("slope-angle", po::value<double>(&config.slopeAngle)->default_value(45.0),
         "Slope angle in degrees")
        ("base-length", po::value<double>(&config.baseLength)->default_value(20.0),
         "Base length for mesh generation")
        ("top-length", po::value<double>(&config.topLength)->default_value(10.0),
         "Top platform length")
        ("mesh-nx", po::value<int>(&config.meshNX)->default_value(20),
         "Number of elements in x direction")
        ("mesh-ny", po::value<int>(&config.meshNY)->default_value(10),
         "Number of elements in y direction");

    po::options_description exportOpt("Export Options");
    exportOpt.add_options()
        ("export-intermediate", po::bool_switch(&config.exportIntermediate),
         "Export intermediate iteration results");

    po::options_description seismic("Seismic Loading Options");
    seismic.add_options()
        ("enable-seismic", po::bool_switch(&config.enableSeismic),
         "Enable seismic loading coupling")
        ("seismic-file", po::value<std::string>(&config.seismicFile),
         "Seismic acceleration time history file")
        ("seismic-start", po::value<double>(&config.seismicStartTime)->default_value(0.0),
         "Seismic loading start time in seconds")
        ("seismic-duration", po::value<double>(&config.seismicDuration)->default_value(30.0),
         "Seismic loading duration for synthetic motion in seconds")
        ("seismic-amplitude", po::value<double>(&config.seismicAmplitude)->default_value(1.0),
         "Seismic acceleration amplitude scaling factor (g)")
        ("seismic-direction", po::value<double>(&config.seismicDirection)->default_value(0.0),
         "Seismic input direction in degrees (0=horizontal only, 90=vertical only)")
        ("seismic-mode", po::value<int>(&config.seismicMode)->default_value(0),
         "Seismic motion generation: 0=harmonic, 1=sine-sweep, 2=artificial motion, 3=external file")
        ("seismic-freq-start", po::value<double>(&config.seismicFreqStart)->default_value(1.0),
         "Sine sweep start frequency in Hz")
        ("seismic-freq-end", po::value<double>(&config.seismicFreqEnd)->default_value(20.0),
         "Sine sweep end frequency in Hz")
        ("seismic-magnitude", po::value<double>(&config.seismicMagnitude)->default_value(6.5),
         "Earthquake magnitude for artificial motion")
        ("seismic-distance", po::value<double>(&config.seismicDistance)->default_value(20.0),
         "Epicentral distance for artificial motion in km")
        ("seismic-harmonic-freq", po::value<double>(&config.seismicHarmonicFreq)->default_value(2.5),
         "Harmonic excitation frequency in Hz")
        ("seismic-harmonic-cycles", po::value<int>(&config.seismicHarmonicCycles)->default_value(10),
         "Number of harmonic cycles (-1 for full duration)")
        ("seismic-apply-x", po::bool_switch(&config.seismicApplyX)->default_value(true),
         "Apply seismic excitation in X direction")
        ("seismic-apply-y", po::bool_switch(&config.seismicApplyY)->default_value(false),
         "Apply seismic excitation in Y direction")
        ("newmark-beta", po::value<double>(&config.newmarkBeta)->default_value(0.25),
         "Newmark-beta integration parameter")
        ("newmark-gamma", po::value<double>(&config.newmarkGamma)->default_value(0.5),
         "Newmark-gamma integration parameter")
        ("rayleigh-alpha", po::value<double>(&config.rayleighAlpha)->default_value(0.05),
         "Rayleigh damping mass coefficient")
        ("rayleigh-beta", po::value<double>(&config.rayleighBeta)->default_value(0.02),
         "Rayleigh damping stiffness coefficient");

    po::options_description all;
    all.add(general).add(time).add(solver).add(meshgen).add(exportOpt).add(seismic);

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, all), vm);
        po::notify(vm);
    } catch (const po::error& e) {
        std::cerr << "Error parsing command line: " << e.what() << std::endl;
        std::cerr << "Use --help for usage information." << std::endl;
        std::exit(1);
    }

    if (vm.count("help")) {
        std::cout << "Slope Stability Simulation - Rainfall-Induced Seepage-Stress Coupling\n\n";
        std::cout << all << std::endl;
        std::exit(0);
    }

    return config;
}

void CLIParser::printHelp()
{
    std::cout << "Usage: slope_sim [options]\n";
    std::cout << "Run 'slope_sim --help' for full options.\n";
}

bool CLIParser::validate(const SimConfig& config)
{
    if (!config.generateMesh && config.meshFile.empty()) {
        std::cerr << "Error: Either --mesh or --generate-mesh must be specified.\n";
        return false;
    }
    if (config.totalTime <= 0.0) {
        std::cerr << "Error: Total time must be positive.\n";
        return false;
    }
    if (config.timeStep <= 0.0) {
        std::cerr << "Error: Time step must be positive.\n";
        return false;
    }
    if (config.timeStep > config.totalTime) {
        std::cerr << "Error: Time step cannot exceed total time.\n";
        return false;
    }
    if (config.maxNewtonIterations < 1) {
        std::cerr << "Error: Max iterations must be >= 1.\n";
        return false;
    }
    if (config.newtonTolerance <= 0.0) {
        std::cerr << "Error: Tolerance must be positive.\n";
        return false;
    }
    if (config.enableSeismic) {
        if (config.seismicMode == 3 && config.seismicFile.empty()) {
            std::cerr << "Error: Seismic mode set to external file but --seismic-file not provided.\n";
            return false;
        }
        if (config.seismicAmplitude < 0.0) {
            std::cerr << "Error: Seismic amplitude factor must be non-negative.\n";
            return false;
        }
        if (config.newmarkBeta <= 0.0 || config.newmarkGamma <= 0.0) {
            std::cerr << "Error: Newmark integration parameters must be positive.\n";
            return false;
        }
    }
    return true;
}

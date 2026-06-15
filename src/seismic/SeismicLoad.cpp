#include "SeismicLoad.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <random>
#include <numbers>

double SeismicLoad::interpolateLinear(double t0, double t1, double v0, double v1, double t) {
    if (std::abs(t1 - t0) < 1.0e-30) return v0;
    double alpha = (t - t0) / (t1 - t0);
    return v0 + alpha * (v1 - v0);
}

bool SeismicLoad::loadTimeHistory(const std::string& filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) return false;

    timeHistory_.clear();
    std::string line;

    while (std::getline(ifs, line)) {
        size_t commentPos = line.find('#');
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }

        std::string trimmed = line;
        trimmed.erase(0, trimmed.find_first_not_of(" \t\r\n"));
        trimmed.erase(trimmed.find_last_not_of(" \t\r\n") + 1);
        if (trimmed.empty()) continue;

        std::istringstream iss(trimmed);
        AccelerationPoint pt{};
        iss >> pt.time >> pt.ax;
        if (iss.good()) {
            iss >> pt.ay;
        } else {
            pt.ay = 0.0;
        }
        timeHistory_.push_back(pt);
    }

    historyIntegrated_ = false;
    velocityHistory_.clear();
    displacementHistory_.clear();
    return !timeHistory_.empty();
}

void SeismicLoad::generateSineSweep(
    double duration, double dt,
    double freqStart, double freqEnd,
    double amplitude, double directionDeg)
{
    timeHistory_.clear();
    double dirRad = directionDeg * M_PI / 180.0;
    double cosD = std::cos(dirRad);
    double sinD = std::sin(dirRad);

    int nSteps = static_cast<int>(std::ceil(duration / dt)) + 1;
    double k = (freqEnd - freqStart) / duration;

    for (int i = 0; i <= nSteps; ++i) {
        double t = i * dt;
        double omega = 2.0 * M_PI * (freqStart + 0.5 * k * t);
        double phi = 2.0 * M_PI * (freqStart * t + 0.5 * k * t * t);
        double a = amplitude * std::sin(phi);

        AccelerationPoint pt{};
        pt.time = t;
        pt.ax = a * cosD;
        pt.ay = a * sinD;
        timeHistory_.push_back(pt);
    }

    historyIntegrated_ = false;
    velocityHistory_.clear();
    displacementHistory_.clear();
}

void SeismicLoad::generateHarmonic(
    double duration, double dt,
    double frequency, double amplitude,
    int numCycles, double directionDeg)
{
    timeHistory_.clear();
    double dirRad = directionDeg * M_PI / 180.0;
    double cosD = std::cos(dirRad);
    double sinD = std::sin(dirRad);
    double omega = 2.0 * M_PI * frequency;

    double totalDuration = duration;
    if (numCycles > 0) {
        totalDuration = numCycles / frequency;
    }

    int nSteps = static_cast<int>(std::ceil(totalDuration / dt)) + 1;
    for (int i = 0; i <= nSteps; ++i) {
        double t = i * dt;
        double a = amplitude * std::sin(omega * t);

        double ramp = 1.0;
        double tRamp = std::min(0.1 * totalDuration, 1.0);
        if (t < tRamp) ramp = t / tRamp;
        if (t > totalDuration - tRamp) ramp = (totalDuration - t) / tRamp;
        a *= std::max(0.0, ramp);

        AccelerationPoint pt{};
        pt.time = t;
        pt.ax = a * cosD;
        pt.ay = a * sinD;
        timeHistory_.push_back(pt);
    }

    historyIntegrated_ = false;
    velocityHistory_.clear();
    displacementHistory_.clear();
}

void SeismicLoad::generateArtificialMotion(
    double duration, double dt,
    double magnitude, double epicenterDistance,
    double directionDeg)
{
    timeHistory_.clear();
    double dirRad = directionDeg * M_PI / 180.0;
    double cosD = std::cos(dirRad);
    double sinD = std::sin(dirRad);

    double PGA = 0.01 * std::pow(10.0, 0.5 * magnitude - 1.0) *
                 std::exp(-0.02 * epicenterDistance);

    std::mt19937 rng(42);
    std::normal_distribution<double> normal(0.0, 1.0);

    int nSteps = static_cast<int>(std::ceil(duration / dt)) + 1;
    int nFreq = nSteps / 2;
    double freqStep = 1.0 / duration;

    std::vector<std::complex<double>> spectrumX(nFreq, std::complex<double>(0, 0));
    std::vector<std::complex<double>> spectrumY(nFreq, std::complex<double>(0, 0));

    double beta1 = 3.0 - magnitude * 0.1;
    double omegaC = 5.0 + magnitude * 2.0;
    double omegaMax = 100.0;

    for (int i = 0; i < nFreq; ++i) {
        double omega = 2.0 * M_PI * i * freqStep;
        if (omega < 1.0e-6) continue;
        if (omega > omegaMax) continue;

        double psd = std::pow(omega / omegaC, beta1) /
                    std::pow(1.0 + std::pow(omega / omegaC, 2.0), (beta1 + 4.0) / 2.0);
        double amp = std::sqrt(psd * freqStep) * PGA;

        double phase = 2.0 * M_PI * normal(rng);
        spectrumX[i] = std::polar(amp * cosD, phase);
        spectrumY[i] = std::polar(amp * sinD, phase + 2.0 * M_PI * normal(rng) * 0.1);
    }

    for (int i = 0; i <= nSteps; ++i) {
        double t = i * dt;
        double ax = 0.0, ay = 0.0;
        for (int j = 1; j < nFreq; ++j) {
            double omega = 2.0 * M_PI * j * freqStep;
            ax += spectrumX[j].real() * std::cos(omega * t) - spectrumX[j].imag() * std::sin(omega * t);
            ay += spectrumY[j].real() * std::cos(omega * t) - spectrumY[j].imag() * std::sin(omega * t);
        }
        ax *= 2.0;
        ay *= 2.0;

        double ramp = 1.0;
        double tRamp = std::min(2.0, 0.1 * duration);
        if (t < tRamp) ramp = t / tRamp;
        if (t > duration - tRamp) ramp = (duration - t) / tRamp;

        AccelerationPoint pt{};
        pt.time = t;
        pt.ax = ax * ramp;
        pt.ay = ay * ramp;
        timeHistory_.push_back(pt);
    }

    double computedPGA = getPGA();
    if (computedPGA > 1.0e-20) {
        double scale = PGA / computedPGA;
        for (auto& pt : timeHistory_) {
            pt.ax *= scale;
            pt.ay *= scale;
        }
    }

    historyIntegrated_ = false;
    velocityHistory_.clear();
    displacementHistory_.clear();
}

int SeismicLoad::findTimeIndex(double time) const {
    if (timeHistory_.empty()) return -1;
    if (time <= timeHistory_.front().time) return 0;
    if (time >= timeHistory_.back().time) return static_cast<int>(timeHistory_.size()) - 2;

    int lo = 0, hi = static_cast<int>(timeHistory_.size()) - 1;
    while (lo <= hi) {
        int mid = (lo + hi) / 2;
        if (timeHistory_[mid].time <= time &&
            (mid + 1 >= static_cast<int>(timeHistory_.size()) || timeHistory_[mid + 1].time > time)) {
            return mid;
        }
        if (timeHistory_[mid].time < time) lo = mid + 1;
        else hi = mid - 1;
    }
    return lo - 1;
}

Eigen::Vector2d SeismicLoad::getAcceleration(double time) const {
    if (timeHistory_.empty()) return Eigen::Vector2d(0.0, 0.0);

    int idx = findTimeIndex(time);
    if (idx < 0) return Eigen::Vector2d(0.0, 0.0);

    int nextIdx = std::min(idx + 1, static_cast<int>(timeHistory_.size()) - 1);
    const auto& p0 = timeHistory_[idx];
    const auto& p1 = timeHistory_[nextIdx];

    double ax = interpolateLinear(p0.time, p1.time, p0.ax, p1.ax, time);
    double ay = interpolateLinear(p0.time, p1.time, p0.ay, p1.ay, time);
    return Eigen::Vector2d(ax, ay);
}

void SeismicLoad::integrateHistory() const {
    if (timeHistory_.size() < 2 || historyIntegrated_) return;

    velocityHistory_.resize(timeHistory_.size());
    displacementHistory_.resize(timeHistory_.size());

    velocityHistory_[0] = {timeHistory_[0].time, 0.0, 0.0};
    displacementHistory_[0] = {timeHistory_[0].time, 0.0, 0.0};

    for (size_t i = 1; i < timeHistory_.size(); ++i) {
        double dt = timeHistory_[i].time - timeHistory_[i - 1].time;
        if (dt < 0.0) dt = 0.0;

        double vx0 = velocityHistory_[i - 1].ax;
        double vy0 = velocityHistory_[i - 1].ay;
        double ax0 = timeHistory_[i - 1].ax;
        double ay0 = timeHistory_[i - 1].ay;
        double ax1 = timeHistory_[i].ax;
        double ay1 = timeHistory_[i].ay;

        velocityHistory_[i].time = timeHistory_[i].time;
        velocityHistory_[i].ax = vx0 + 0.5 * dt * (ax0 + ax1);
        velocityHistory_[i].ay = vy0 + 0.5 * dt * (ay0 + ay1);

        double dx0 = displacementHistory_[i - 1].ax;
        double dy0 = displacementHistory_[i - 1].ay;

        displacementHistory_[i].time = timeHistory_[i].time;
        displacementHistory_[i].ax = dx0 + dt * vx0 + dt * dt / 6.0 * (2.0 * ax0 + ax1);
        displacementHistory_[i].ay = dy0 + dt * vy0 + dt * dt / 6.0 * (2.0 * ay0 + ay1);
    }

    historyIntegrated_ = true;
}

Eigen::Vector2d SeismicLoad::getVelocity(double time) const {
    integrateHistory();
    if (velocityHistory_.empty()) return Eigen::Vector2d(0.0, 0.0);

    int idx = findTimeIndex(time);
    if (idx < 0) return Eigen::Vector2d(0.0, 0.0);

    int nextIdx = std::min(idx + 1, static_cast<int>(velocityHistory_.size()) - 1);
    const auto& p0 = velocityHistory_[idx];
    const auto& p1 = velocityHistory_[nextIdx];

    double vx = interpolateLinear(p0.time, p1.time, p0.ax, p1.ax, time);
    double vy = interpolateLinear(p0.time, p1.time, p0.ay, p1.ay, time);
    return Eigen::Vector2d(vx, vy);
}

Eigen::Vector2d SeismicLoad::getDisplacement(double time) const {
    integrateHistory();
    if (displacementHistory_.empty()) return Eigen::Vector2d(0.0, 0.0);

    int idx = findTimeIndex(time);
    if (idx < 0) return Eigen::Vector2d(0.0, 0.0);

    int nextIdx = std::min(idx + 1, static_cast<int>(displacementHistory_.size()) - 1);
    const auto& p0 = displacementHistory_[idx];
    const auto& p1 = displacementHistory_[nextIdx];

    double dx = interpolateLinear(p0.time, p1.time, p0.ax, p1.ax, time);
    double dy = interpolateLinear(p0.time, p1.time, p0.ay, p1.ay, time);
    return Eigen::Vector2d(dx, dy);
}

Eigen::VectorXd SeismicLoad::computeInertialForce(
    const Eigen::VectorXd& nodeMass,
    double time) const
{
    Eigen::Vector2d acc = getAcceleration(time);
    int n = static_cast<int>(nodeMass.size());
    Eigen::VectorXd force(2 * n);
    force.setZero();

    for (int i = 0; i < n; ++i) {
        force(2 * i) = -nodeMass(i) * acc(0);
        force(2 * i + 1) = -nodeMass(i) * acc(1);
    }

    return force;
}

void SeismicLoad::assembleInertialForceToNodes(
    const std::vector<double>& nodeX,
    const std::vector<double>& nodeY,
    const std::vector<int>& elementNodeIds,
    int nodesPerElement,
    const std::vector<double>& elementUnitWeight,
    const std::vector<double>& elementArea,
    double currentTime,
    Eigen::VectorXd& fx,
    Eigen::VectorXd& fy,
    double amplitudeScale) const
{
    int numNodes = static_cast<int>(nodeX.size());
    int numElements = static_cast<int>(elementNodeIds.size()) / nodesPerElement;

    fx.setZero(numNodes);
    fy.setZero(numNodes);

    Eigen::Vector2d acc = getAcceleration(currentTime);
    acc *= amplitudeScale;

    for (int e = 0; e < numElements; ++e) {
        double gamma = elementUnitWeight[e];
        double A = elementArea[e];
        double nodeMass = gamma * A / (3.0 * 9.81);

        for (int n = 0; n < nodesPerElement; ++n) {
            int nid = elementNodeIds[e * nodesPerElement + n];
            fx(nid) += -nodeMass * acc(0);
            fy(nid) += -nodeMass * acc(1);
        }
    }
}

Eigen::Matrix<double, 6, 6> SeismicLoad::computeElementMassMatrix(
    const Eigen::Vector2d& p0,
    const Eigen::Vector2d& p1,
    const Eigen::Vector2d& p2,
    double unitWeight,
    bool lumped) const
{
    double A = 0.5 * std::abs((p1(0) - p0(0)) * (p2(1) - p0(1)) - (p2(0) - p0(0)) * (p1(1) - p0(1)));
    double rho = unitWeight / 9.81;
    double massPerNode = rho * A / 3.0;

    Eigen::Matrix<double, 6, 6> M = Eigen::Matrix<double, 6, 6>::Zero();
    if (lumped) {
        for (int i = 0; i < 3; ++i) {
            M(2 * i, 2 * i) = massPerNode;
            M(2 * i + 1, 2 * i + 1) = massPerNode;
        }
    } else {
        double coeff = rho * A / 12.0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double mij = (i == j) ? 2.0 * coeff : coeff;
                M(2 * i, 2 * j) = mij;
                M(2 * i + 1, 2 * j + 1) = mij;
            }
        }
    }
    return M;
}

double SeismicLoad::getDuration() const {
    if (timeHistory_.empty()) return 0.0;
    return timeHistory_.back().time - timeHistory_.front().time;
}

double SeismicLoad::getPGA() const {
    double pga = 0.0;
    for (const auto& pt : timeHistory_) {
        double a = std::sqrt(pt.ax * pt.ax + pt.ay * pt.ay);
        if (a > pga) pga = a;
    }
    return pga;
}

Eigen::VectorXd SeismicLoad::computeNodeLumpedMass(
    const std::vector<double>& nodeX,
    const std::vector<double>& nodeY,
    const std::vector<int>& elementNodeIds,
    int nodesPerElement,
    const std::vector<double>& elementUnitWeight,
    double thickness)
{
    int numNodes = static_cast<int>(nodeX.size());
    int numElements = static_cast<int>(elementNodeIds.size()) / nodesPerElement;

    Eigen::VectorXd mass(numNodes);
    mass.setZero();

    for (int e = 0; e < numElements; ++e) {
        int n0 = elementNodeIds[e * nodesPerElement];
        int n1 = elementNodeIds[e * nodesPerElement + 1];
        int n2 = elementNodeIds[e * nodesPerElement + 2];

        double x0 = nodeX[n0], y0 = nodeY[n0];
        double x1 = nodeX[n1], y1 = nodeY[n1];
        double x2 = nodeX[n2], y2 = nodeY[n2];

        double A = 0.5 * std::abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
        double rho = elementUnitWeight[e] / 9.81;
        double nodeMass = rho * A * thickness / 3.0;

        for (int n = 0; n < nodesPerElement; ++n) {
            int nid = elementNodeIds[e * nodesPerElement + n];
            mass(nid) += nodeMass;
        }
    }

    return mass;
}

#include "MeshGenerator.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>

void Mesh::addNode(const Node& node) {
    nodes_.push_back(node);
}

void Mesh::addElement(const TriangleElement& elem) {
    elements_.push_back(elem);
}

void Mesh::addZone(int zoneId, const ZoneProperties& props) {
    zones_[zoneId] = props;
}

void Mesh::addBoundaryCondition(const BoundaryCondition& bc) {
    boundaryConditions_.push_back(bc);
}

Eigen::Vector2d Mesh::getNodeCoord(int nodeId) const {
    for (const auto& n : nodes_) {
        if (n.id == nodeId) return Eigen::Vector2d(n.x, n.y);
    }
    throw std::runtime_error("Node ID not found: " + std::to_string(nodeId));
}

const ZoneProperties& Mesh::getZoneProperties(int zoneId) const {
    auto it = zones_.find(zoneId);
    if (it == zones_.end())
        throw std::runtime_error("Zone ID not found: " + std::to_string(zoneId));
    return it->second;
}

void MeshGenerator::triangulateQuad(
    int n0, int n1, int n2, int n3,
    int zoneId, int& elemId,
    std::vector<TriangleElement>& elements)
{
    TriangleElement e1, e2;
    e1.id = elemId++;
    e1.nodeIds[0] = n0; e1.nodeIds[1] = n1; e1.nodeIds[2] = n2;
    e1.zoneId = zoneId;
    e2.id = elemId++;
    e2.nodeIds[0] = n0; e2.nodeIds[1] = n2; e2.nodeIds[2] = n3;
    e2.zoneId = zoneId;
    elements.push_back(e1);
    elements.push_back(e2);
}

std::unique_ptr<Mesh> MeshGenerator::generateSimpleSlope(
    double height, double slopeAngleDeg,
    double baseLength, double topLength,
    int nx, int ny)
{
    auto mesh = std::make_unique<Mesh>();
    double slopeRad = slopeAngleDeg * M_PI / 180.0;
    double slopeHRun = height / std::tan(slopeRad);
    double totalWidth = baseLength + slopeHRun + topLength;

    int nodeId = 0;
    for (int j = 0; j <= ny; ++j) {
        double y = height * static_cast<double>(j) / ny;
        for (int i = 0; i <= nx; ++i) {
            double x = totalWidth * static_cast<double>(i) / nx;
            Node n;
            n.id = nodeId++;
            n.x = x;
            n.y = y;
            mesh->addNode(n);
        }
    }

    int cols = nx + 1;
    int elemId = 0;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int n0 = j * cols + i;
            int n1 = j * cols + i + 1;
            int n2 = (j + 1) * cols + i + 1;
            int n3 = (j + 1) * cols + i;

            double cx = 0.0, cy = 0.0;
            for (int nid : {n0, n1, n2, n3}) {
                cx += mesh->getNodes()[nid].x;
                cy += mesh->getNodes()[nid].y;
            }
            cx /= 4.0; cy /= 4.0;

            int zoneId = 1;
            double slopeStartX = baseLength;
            double slopeEndX = baseLength + slopeHRun;
            if (cy < height && cx >= slopeStartX && cx <= slopeEndX) {
                zoneId = 2;
            }

            triangulateQuad(n0, n1, n2, n3, zoneId, elemId,
                            const_cast<std::vector<TriangleElement>&>(mesh->getElements()));
        }
    }

    ZoneProperties zone1{};
    zone1.zoneId = 1;
    zone1.kx = 1.0e-6; zone1.ky = 1.0e-6;
    zone1.cohesion = 25.0e3;
    zone1.frictionAngle = 28.0;
    zone1.dilationAngle = 0.0;
    zone1.unitWeight = 18.5e3;
    zone1.porosity = 0.35;
    zone1.vanGenuchtenAlpha = 0.5;
    zone1.vanGenuchtenN = 1.3;
    zone1.residualSaturation = 0.1;
    mesh->addZone(1, zone1);

    ZoneProperties zone2{};
    zone2.zoneId = 2;
    zone2.kx = 5.0e-6; zone2.ky = 5.0e-6;
    zone2.cohesion = 15.0e3;
    zone2.frictionAngle = 22.0;
    zone2.dilationAngle = 0.0;
    zone2.unitWeight = 19.0e3;
    zone2.porosity = 0.40;
    zone2.vanGenuchtenAlpha = 0.8;
    zone2.vanGenuchtenN = 1.2;
    zone2.residualSaturation = 0.15;
    mesh->addZone(2, zone2);

    return mesh;
}

std::unique_ptr<Mesh> MeshGenerator::loadFromJSON(const std::string& filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open())
        throw std::runtime_error("Cannot open mesh file: " + filename);

    auto mesh = std::make_unique<Mesh>();
    std::string line, section;

    while (std::getline(ifs, line)) {
        std::string trimmed = line;
        trimmed.erase(0, trimmed.find_first_not_of(" \t\r\n"));
        trimmed.erase(trimmed.find_last_not_of(" \t\r\n") + 1);
        if (trimmed.empty() || trimmed[0] == '#') continue;

        if (trimmed == "[NODES]") { section = "NODES"; continue; }
        if (trimmed == "[ELEMENTS]") { section = "ELEMENTS"; continue; }
        if (trimmed == "[ZONES]") { section = "ZONES"; continue; }
        if (trimmed == "[BOUNDARY]") { section = "BOUNDARY"; continue; }

        std::istringstream iss(trimmed);

        if (section == "NODES") {
            Node n;
            iss >> n.id >> n.x >> n.y;
            mesh->addNode(n);
        } else if (section == "ELEMENTS") {
            TriangleElement e;
            iss >> e.id >> e.nodeIds[0] >> e.nodeIds[1] >> e.nodeIds[2] >> e.zoneId;
            mesh->addElement(e);
        } else if (section == "ZONES") {
            ZoneProperties z{};
            iss >> z.zoneId >> z.kx >> z.ky
                >> z.cohesion >> z.frictionAngle >> z.dilationAngle
                >> z.unitWeight >> z.porosity
                >> z.vanGenuchtenAlpha >> z.vanGenuchtenN >> z.residualSaturation;
            mesh->addZone(z.zoneId, z);
        } else if (section == "BOUNDARY") {
            BoundaryCondition bc;
            std::string typeStr;
            iss >> bc.nodeId >> typeStr >> bc.value >> bc.startTime >> bc.endTime;
            bc.type = typeStr;
            mesh->addBoundaryCondition(bc);
        }
    }

    return mesh;
}

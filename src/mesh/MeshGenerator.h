#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <map>
#include <memory>

struct Node {
    int id;
    double x, y;
};

struct TriangleElement {
    int id;
    int nodeIds[3];
    int zoneId;
};

struct ZoneProperties {
    int zoneId;
    double kx;
    double ky;
    double cohesion;
    double frictionAngle;
    double dilationAngle;
    double unitWeight;
    double porosity;
    double vanGenuchtenAlpha;
    double vanGenuchtenN;
    double residualSaturation;
};

struct BoundaryCondition {
    int nodeId;
    std::string type;
    double value;
    double startTime;
    double endTime;
};

class Mesh {
public:
    Mesh() = default;

    void addNode(const Node& node);
    void addElement(const TriangleElement& elem);
    void addZone(int zoneId, const ZoneProperties& props);
    void addBoundaryCondition(const BoundaryCondition& bc);

    const std::vector<Node>& getNodes() const { return nodes_; }
    const std::vector<TriangleElement>& getElements() const { return elements_; }
    const std::map<int, ZoneProperties>& getZones() const { return zones_; }
    const std::vector<BoundaryCondition>& getBoundaryConditions() const { return boundaryConditions_; }

    Eigen::Vector2d getNodeCoord(int nodeId) const;
    const ZoneProperties& getZoneProperties(int zoneId) const;
    int numNodes() const { return static_cast<int>(nodes_.size()); }
    int numElements() const { return static_cast<int>(elements_.size()); }

private:
    std::vector<Node> nodes_;
    std::vector<TriangleElement> elements_;
    std::map<int, ZoneProperties> zones_;
    std::vector<BoundaryCondition> boundaryConditions_;
};

class MeshGenerator {
public:
    static std::unique_ptr<Mesh> loadFromJSON(const std::string& filename);
    static std::unique_ptr<Mesh> generateSimpleSlope(
        double height, double slopeAngleDeg,
        double baseLength, double topLength,
        int nx, int ny);

private:
    static void triangulateQuad(
        int n0, int n1, int n2, int n3,
        int zoneId, int& elemId,
        std::vector<TriangleElement>& elements);
};

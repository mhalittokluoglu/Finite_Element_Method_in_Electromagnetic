#include "Node.h"
#include <vector>
#ifndef _NODES_
#define _NODES_
struct Nodes
{
    std::vector<std::vector<Node>> nodeList;

    Nodes(const std::vector<std::vector<double>> &NodeVoltages,
            const std::vector<std::vector<BoundaryConditionValue>> &Boundaries);
    std::vector<Node> operator[](const size_t &index) const;
    size_t size() const;
    void PrintNodes() const;
    void PrintNodeVoltages() const;
    size_t getNumberOfNodes() const;
private:
    size_t NumberOfNodes;
    
};
#endif
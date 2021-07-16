#include "Nodes.h"
#include <iostream>

Nodes::Nodes(const std::vector<std::vector<double>> &NodeVoltages,
const std::vector<std::vector<BoundaryConditionValue>> &Boundaries)
{
    std::vector<Node> tempNodes;
    int i=0, j;
    int NodeNumber = 0;
    for(auto row : NodeVoltages)
    {
        j = 0;
        for (auto nodeVoltage :row)
        {
            if(Boundaries[i][j] == NOT_ON_WORKSPACE)
            {
                tempNodes.push_back(Node(-1,nodeVoltage,
                (double)j, (double)-i,Boundaries[i][j]));
                j++;
            }
            else
            {
                tempNodes.push_back(Node(NodeNumber,nodeVoltage,
                (double)j, (double)-i,Boundaries[i][j]));
                j++;
                NodeNumber++;
            }
        }
        nodeList.push_back(tempNodes);
        tempNodes.clear();
        i++;
    }
    NumberOfNodes = NodeNumber;
}
std::vector<Node> Nodes::operator[](const size_t &index) const
{
    return nodeList[index];
}
size_t Nodes::size() const
{
    return nodeList.size();
}

void Nodes::PrintNodes() const
{
    std::cout << "\nNODE NUMBERS:\n";
    for (auto row : nodeList)
    {
        for (auto node : row)
        {
            if (node.getIsOnWorkSpace())
                std::cout << node.NodeNumber << "\t";
        }
        std::cout << "\n";
    }
}
void Nodes::PrintNodeVoltages() const
{
    for (auto row : nodeList)
    {
        for (auto node : row)
        {
            printf("%.3lf\t",node.NodeVoltage);
        }
        printf("\n");
    }
}

size_t Nodes::getNumberOfNodes() const
{
    return NumberOfNodes;
}
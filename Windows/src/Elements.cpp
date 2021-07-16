#include "Elements.h"
#include <iostream>

Elements::Elements(Nodes &nodes)
{
    std::vector<Element> tempElements;
    bool *availableNodes;
    int ElementNumber = 0;
    for (size_t i = 0; i < nodes.nodeList.size(); i++)
    {
        for (size_t j = 0; j < nodes.nodeList[i].size(); j++)
        {
            if (nodes.nodeList[i][j].getIsOnWorkSpace())
            {
                availableNodes = checkNode(nodes,i,j);
                if (availableNodes[0] && availableNodes[1])
                {
                    tempElements.push_back(Element(ElementNumber, &nodes.nodeList[i][j],&nodes.nodeList[i+1][j-1], &nodes.nodeList[i+1][j]));
                    ElementNumber++;
                }
                if (availableNodes[1] && availableNodes[2])
                {
                    tempElements.push_back(Element(ElementNumber, &nodes.nodeList[i][j],&nodes.nodeList[i+1][j], &nodes.nodeList[i][j+1]));
                    ElementNumber++;
                }
                else if(availableNodes[1] && availableNodes[3])
                {
                    tempElements.push_back(Element(ElementNumber, &nodes.nodeList[i][j],&nodes.nodeList[i+1][j], &nodes.nodeList[i+1][j+1]));
                    ElementNumber++;
                }
                delete[] availableNodes;
            }
        }
        elementList.push_back(tempElements);
        tempElements.clear();
    }
    NumberOfElements = ElementNumber;
    NumberOfNodes = nodes.getNumberOfNodes();
}

bool* Elements::checkNode(Nodes &nodes,size_t i,size_t j)
{
    bool* availableNodes = new bool[4];
    bool isDownPossible = i < nodes.nodeList.size() - 1;
    bool isLeftDownPossible = false;
    
    if(isDownPossible)
        isLeftDownPossible = (j > 0) && j < nodes.nodeList[i + 1].size();
    
    bool isRightPossible = j < nodes.nodeList[i].size() - 1;
    bool isRightDownPossible = false;

    if(isDownPossible)
        isRightDownPossible = j < nodes.nodeList[i+1].size() - 1;

    if (isLeftDownPossible)
    {
        isLeftDownPossible = nodes.nodeList[i+1][j-1].getIsOnWorkSpace();
    }
    if (isRightPossible)
    {
        isRightPossible = nodes.nodeList[i][j+1].getIsOnWorkSpace();
    }
    if (isRightDownPossible)
    {
        isRightDownPossible = nodes.nodeList[i+1][j+1].getIsOnWorkSpace();
    }
    if (isDownPossible)
    {
        isDownPossible = nodes.nodeList[i+1][j].getIsOnWorkSpace();
    }

    availableNodes[0] = isLeftDownPossible;
    availableNodes[1] = isDownPossible;
    availableNodes[2] = isRightPossible;
    availableNodes[3] = isRightDownPossible;
    return availableNodes;
}

void Elements::PrintElements() const
{
    std::cout << "\n - - - - - - - - - - - - - - - - - - - - - - \n";
    std::cout << "ELEMENTS:\n";
    std::cout << "\n - - - - - - - - - - - - - - - - - - - - - - \n";
    for (auto row : elementList)
    {
        for(const auto &e: row)
        {
            int _NodeNumbers[3];
            e.getNodeNumbers(_NodeNumbers);
            std::cout << e.getElementNumber() << ": " <<
            _NodeNumbers[0]<< ", " <<_NodeNumbers[1]
            << ", " << _NodeNumbers[2] << std::endl;
        }
        std::cout << "\n - - - - - - - - - - - - - - - - - - - - - - \n";
    }
}

size_t Elements::getNumberOfElements() const
{
    return NumberOfElements;
}

size_t Elements::getNumberOfNodes() const
{
    return NumberOfNodes;
}
#include "../include/Element.h"

Element::Element(int ElementNo, Node *node1, Node *node2, Node *node3)
{
    ElementNumber = ElementNo;
    threeNode.nodes[0] = node1;
    threeNode.nodes[1] = node2;
    threeNode.nodes[2] = node3;
}

int Element::getElementNumber() const
{
    return ElementNumber;
}

void Element::getNodeNumbers(int *NodeNumbers) const
{
    NodeNumbers[0] = threeNode.nodes[0]->NodeNumber;
    NodeNumbers[1] = threeNode.nodes[1]->NodeNumber;
    NodeNumbers[2] = threeNode.nodes[2]->NodeNumber;
}

ThreeNodesStruct Element::getNodes() const
{
    return threeNode;
}
#include "Node.h"
#ifndef _ELEMENT_
#define _ELEMENT_

struct ThreeNodesStruct
{
private:
    Node* nodes[3];

    friend class Element;
    friend class Elements;
    friend class FiniteElements;
};

struct Element
{
private:
    ThreeNodesStruct threeNode;
    int ElementNumber;
public:
    Element(int ElementNo, Node *node1, Node *node2, Node *node3);
    int getElementNumber() const;
    void getNodeNumbers(int *NodeNumbers) const;
    ThreeNodesStruct getNodes() const;

};

#endif
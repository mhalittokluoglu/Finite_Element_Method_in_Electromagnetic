#include "Nodes.h"
#include "Element.h"
#include <vector>
#ifndef _ELEMENTS_
#define _ELEMENTS
class Elements
{
private:
    std::vector<std::vector<Element>> elementList;
    size_t NumberOfElements;
    size_t NumberOfNodes;
    bool* checkNode(Nodes &nodeList,size_t i,size_t j);
public:
    Elements(Nodes &nodes);
    void PrintElements() const;
    size_t getNumberOfElements() const;
    size_t getNumberOfNodes() const;
    

    friend class FiniteElements;
};
#endif
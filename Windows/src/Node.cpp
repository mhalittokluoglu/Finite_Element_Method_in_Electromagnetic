#include "Node.h"

Node::Node(int _NodeNumber, double _NodeVoltage,double x, double y,
BoundaryConditionValue _BoundaryCondition = NOT_ON_WORKSPACE)
{
    isOnWorkSpace = true;
    NodeNumber = _NodeNumber;
    NodeVoltage = _NodeVoltage;
    this->x = x;
    this->y = y;
    BoundaryCondition  =_BoundaryCondition;
    if(BoundaryCondition == NOT_ON_WORKSPACE)
    {
        isOnWorkSpace = false;
    } 
}

void Node::setBoundary(BoundaryConditionValue condition)
{
    isOnWorkSpace = true;
    BoundaryCondition = condition;
    if(BoundaryCondition == NOT_ON_WORKSPACE)
    {
        isOnWorkSpace  = false;
    }
}

BoundaryConditionValue Node::getBoundaryCondition() const
{
    return BoundaryCondition;
}

bool Node::getIsOnWorkSpace() const
{
    return isOnWorkSpace;
}
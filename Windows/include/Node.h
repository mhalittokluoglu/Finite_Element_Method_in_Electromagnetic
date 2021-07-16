#ifndef _NODE_
#define _NODE_
enum BoundaryConditionValue
{
    BOUNDARY,
    NOT_BOUNDARY,
    NOT_ON_WORKSPACE
};
struct Node 
{
    int NodeNumber;
    double NodeVoltage;
    double x;
    double y;
    Node(int _NodeNumber, double _NodeVoltage, double x, double y, 
     BoundaryConditionValue _BoundaryCondition);
private:
    BoundaryConditionValue BoundaryCondition;
    bool isOnWorkSpace;
public:
    void setBoundary(BoundaryConditionValue condition);
    BoundaryConditionValue getBoundaryCondition() const;
    bool getIsOnWorkSpace() const;

    friend class Element;
};
#endif
#include "Elements.h"

#ifndef _FINITEELEMENTS_
#define _FINITEELEMENTS_

class FiniteElements
{
private:
    std::vector<std::vector<double>> K;
    std::vector<double> b;
    std::vector<double> X;
    bool isSolved;
    size_t NumberOfElements;
    size_t NumberOfNodes;
    Nodes *nodesPtr;
    Elements *elementsPtr;

private:
    double vNorm(const std::vector<double> &Z, const int N);
    void prod(const std::vector<std::vector<double>> &_K,
            const std::vector<double> &_X, std::vector<double> &AX,
            const int control, const int N);
    void setNodeVoltages();

public:
    FiniteElements(Elements *elements, Nodes *nodes);
    void Solve(double Tolerance);
    void PrintK_Matrix() const;
    void PrintK_Diagonal_Matrix() const;
    void PrintResults() const;
};

#endif
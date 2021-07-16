#include <iostream>
#include <vector>
#include "FiniteElements.h"
void set_initial_conditions(std::vector<std::vector<double>> &InitialNodeVoltages, std::vector<std::vector<BoundaryConditionValue>> &RegionOfWork)
{
    for (int i = 0; i < InitialNodeVoltages.size(); i++)
    {
        for (int j = 0; j < InitialNodeVoltages[0].size(); j++)
        {
            if (j == 0 && i != InitialNodeVoltages.size() - 1) // setting the "a" line 
            {
                InitialNodeVoltages[i][j] = 20.0;
                RegionOfWork[i][j] = BOUNDARY;
            }
            else if (i == j) // setting the "c" line
            {
                InitialNodeVoltages[i][j] = 10.0;
                RegionOfWork[i][j] = BOUNDARY;
            }
            else if(i == InitialNodeVoltages.size()-1)
            {
                InitialNodeVoltages[i][j] = 0.0;
                RegionOfWork[i][j] = BOUNDARY;
            }
            else if(i > j)
            {
                RegionOfWork[i][j] = NOT_BOUNDARY;
            }
            else
            {
                RegionOfWork[i][j] = NOT_ON_WORKSPACE;
            }
        }
    }
}
void PrintRegionOfWork(const std::vector<std::vector<BoundaryConditionValue>> &RegionOfWork)
{
    std::cout << " ------------------------------------------------------------------\n";
    std::cout << "\nRegion of Interest: 1-> Internal node, 0-> Boundary node, 2-> Outside\n\n";
    for (auto row : RegionOfWork)
    {
        for (auto x : row)
        {
            std::cout << (int)x << "  ";
        }
        std::cout << "\n";
    }
}

void PrintInitialNodeVoltages(const std::vector<std::vector<double>> &InitialNodeVoltages,
const std::vector<std::vector<BoundaryConditionValue>> &RegionOfWork)
{
    std::cout << " ------------------------------------------------------------------\n";
    std::cout << "\n Initial NODE VOLTAGES:\n\n";

    for (size_t i = 0; i < InitialNodeVoltages.size(); i++)
    {
        for (size_t j = 0; j < InitialNodeVoltages[i].size(); j++)
        {
            if (RegionOfWork[i][j] != NOT_ON_WORKSPACE)
            std::cout << InitialNodeVoltages[i][j] << "\t";
        }
        std::cout << "\n";
    }
}

int main()
{
    int rows = 13, columns = 13;
    std::vector<std::vector<double>> InitialNodeVoltages;
    std::vector<std::vector<BoundaryConditionValue>> RegionOfWork;


    InitialNodeVoltages.resize(rows,std::vector<double>(columns,0));
    RegionOfWork.resize(rows,std::vector<BoundaryConditionValue>(columns,NOT_ON_WORKSPACE));

    set_initial_conditions(InitialNodeVoltages,RegionOfWork);
    PrintRegionOfWork(RegionOfWork);
    PrintInitialNodeVoltages(InitialNodeVoltages, RegionOfWork);

    // Create Nodes
    Nodes nodes = Nodes(InitialNodeVoltages,RegionOfWork);
    nodes.PrintNodes();

    // Create Elements From Nodes
    Elements _elements = Elements(nodes);
    _elements.PrintElements();

    // Apply Finite Element Method
    FiniteElements finiteElement = FiniteElements(&_elements, &nodes);
    finiteElement.PrintK_Diagonal_Matrix();
    
    // Solve
    double Tolerance = 1E-6; // Max error which is allowed
    finiteElement.Solve(Tolerance);
    finiteElement.PrintResults();

}
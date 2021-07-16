#include "../include/FiniteElements.h"
#include <iostream>
#include <cmath>

FiniteElements::FiniteElements(Elements *elements, Nodes *nodes)
{
    isSolved = false;
    nodesPtr = nodes;
    elementsPtr = elements;
    NumberOfElements = elementsPtr->getNumberOfElements();
    NumberOfNodes = elementsPtr->getNumberOfNodes();
    K.resize(NumberOfNodes,std::vector<double>(NumberOfNodes,0.0));
    b.resize(NumberOfNodes,0.0);
    double be[3], ce[3], deltae, Ke[3][3];
    for (auto rows: elementsPtr->elementList)
    {
        for (const auto &e: rows)
        {
            be[0] = e.getNodes().nodes[1]->y - e.getNodes().nodes[2]->y;
            be[1] = e.getNodes().nodes[2]->y - e.getNodes().nodes[0]->y;
            be[2] = e.getNodes().nodes[0]->y - e.getNodes().nodes[1]->y;

            ce[0] = e.getNodes().nodes[2]->x - e.getNodes().nodes[1]->x;
            ce[1] = e.getNodes().nodes[0]->x - e.getNodes().nodes[2]->x;
            ce[2] = e.getNodes().nodes[1]->x - e.getNodes().nodes[0]->x;

            deltae=0.5*(be[0]*ce[1]-be[1]*ce[0]);

            for (int i=0;i<3;i++) 
            {
                for (int j=0;j<3;j++) 
                {
                    Ke[i][j] = (be[i]*be[j]+ce[i]*ce[j])/(4.0*deltae);
                }
            }
            for (int i=0;i<3;i++)
            {
                for (int j=0;j<3;j++)
                {
                    int _NodeNumbers[3];
                    e.getNodeNumbers(_NodeNumbers);
                    K[_NodeNumbers[i]][_NodeNumbers[j]] += Ke[i][j];
                }
            }
        }
    }

    for (auto rows: nodesPtr->nodeList)
    {
        for(const auto &n: rows)
        {
            if (n.getBoundaryCondition() == BOUNDARY)
            {
                b[n.NodeNumber] = n.NodeVoltage;
                K[n.NodeNumber][n.NodeNumber] = 1.0;
                for (int j = 0; j < NumberOfNodes; j++)
                {
                    if (j != n.NodeNumber)
                    {
                        b[j] -= K[j][n.NodeNumber] * n.NodeVoltage;
                        K[n.NodeNumber][j] = 0.0;
                        K[j][n.NodeNumber] = 0.0;
                    }
                }
            }
        }
    }
}

void FiniteElements::Solve(double Tolerance)
{
    std::vector<double> P;
    std::vector<double> R;
    std::vector<double> AX;
    double vnrm, alpha, beta, RSS;
    int ITMAX = 300, iter = 0;
    AX.resize(NumberOfNodes,0.0);
    P.resize(NumberOfNodes,0.0);
    R.resize(NumberOfNodes,0.0);
    vnrm = vNorm(b,NumberOfNodes);
    X.resize(NumberOfNodes,0.0);
    prod(K, X, AX, 0, NumberOfNodes);
    for (int i = 0; i < NumberOfNodes; i++)
    {
        R[i] = b[i] - AX[i];
    }
    prod(K, R, AX, 1, NumberOfNodes);
    beta = 1.0 / vNorm(AX,NumberOfNodes);
    for (int i = 0; i < NumberOfNodes; i++)
    {
        P[i] = beta*AX[i];
    }
    while (iter <= ITMAX)
    {
        prod(K,P,AX,0,NumberOfNodes);
        alpha = 1.0 / vNorm(AX, NumberOfNodes);
        for (int i=0; i < NumberOfNodes; i++)
        {
            X[i] = X[i] + alpha * P[i];
            R[i] = R[i] - alpha * AX[i];
        }
        prod(K,R,AX,1,NumberOfNodes);
        beta = 1.0 / vNorm(AX,NumberOfNodes);
        
        for (int i = 0; i < NumberOfNodes; i++)
        {
            P[i] += beta * AX[i];
        }
        iter++;
        RSS = sqrt(vNorm(R,NumberOfNodes)/vnrm);
        if(RSS <= Tolerance)
        {
            std::cout << "\nConvergence achieved, iteration: " << iter << std::endl;
            isSolved = true;
            setNodeVoltages();
            return;
        }
        else
        {
            if(iter == ITMAX)
            {
                std::cout << "\nMaximum iteration exceeded. No convergence." << std::endl;
                exit(0);
            }
        }
    }

}

double FiniteElements::vNorm(const std::vector<double> &Z, const int N)
{
    double result;
    result = 0.0;
    for (int i =0; i < N; i++)
    {
        result += pow(Z[i],2);
    }
    return result;
}
void FiniteElements::prod(const std::vector<std::vector<double>> &_K,
const std::vector<double> &_X, std::vector<double> &AX,
const int control, const int N)
{
    for(int i = 0; i < N; i++)
    {
        AX[i] = 0.0;
    }
    if (control ==1)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                AX[i] += _K[j][i]*_X[j];
            }
        }
    }
    else
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                AX[j] += _K[j][i]*_X[i];
            }
        }
    }
}
void FiniteElements::setNodeVoltages()
{
    for (int i = 0; i < nodesPtr->nodeList.size(); i++)
    {
        for (int j = 0; j < nodesPtr->nodeList[i].size(); j++)
        {
            if(nodesPtr->nodeList[i][j].NodeNumber >= 0)
                nodesPtr->nodeList[i][j].NodeVoltage = X[nodesPtr->nodeList[i][j].NodeNumber];
        }
    }
}

void FiniteElements::PrintK_Matrix() const
{
    for(int i = 0; i < NumberOfNodes; i++)
    {
        for (int j = 0; j < NumberOfNodes; j++)
        {
            std::cout << "K[" << i << "][" << j << "] = " << K[i][j] << std::endl;
        }
    }
}

void FiniteElements::PrintK_Diagonal_Matrix() const
{
    for(int i = 0; i < NumberOfNodes; i++)
    {
        for (int j = 0; j < NumberOfNodes; j++)
        {
            if(i==j)
                std::cout << "K[" << i << "][" << j << "] = " << K[i][j] << std::endl;
        }
    }
}
void FiniteElements::PrintResults() const
{
    if (isSolved)
    {
        for (int i = 0; i < nodesPtr->nodeList.size(); i++)
        {
            for (int j = 0; j < nodesPtr->nodeList[i].size(); j++)
            {
                if(nodesPtr->nodeList[i][j].getIsOnWorkSpace())
                    printf("%.3lf\t",nodesPtr->nodeList[i][j].NodeVoltage);
            }
            printf("\n");
        }
    }
}
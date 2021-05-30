#include <iostream>
#include <vector>
#include <cmath>

class Node
{
private:
    bool isNumbered,isBoundary;
    double NodeVoltage,x,y;
    int NodeNumber;
public:
    Node(double _NodeVoltage,double _x, double _y, bool _isBoundary = false, int _NodeNumber = -1, bool _isNumbered = false)
    {
        NodeVoltage = _NodeVoltage;
        x = _x;
        y = _y;
        isBoundary = _isBoundary;
        NodeNumber = _NodeNumber;
        isNumbered = _isNumbered;
        NodeNumber++;
    }
    void setNodeVoltage(double _VoltageValue)
    {
        NodeVoltage = _VoltageValue;
    }
    double getNodeVoltage() const
    {
        return NodeVoltage;
    }
    void setNodeNumber(int _NodeNumber)
    {
        NodeNumber = _NodeNumber;
        set_isNumbered(true);
    }
    int getNodeNumber() const
    {
        return NodeNumber;
    }
    void set_isNumbered(bool _isNumbered) 
    {
        isNumbered = _isNumbered;
    }
    bool get_isNumbered() const
    {
        return isNumbered;
    }
    void set_isBoundary(bool _isBoundary) 
    {
        isBoundary = _isBoundary;
    }
    bool get_isBoundary() const
    {
        return isBoundary;
    }
    double get_xPosition() const
    {
        return x;
    }
    double get_yPosition() const
    {
        return y;
    }
};


class Nodes
{
protected:
    std::vector<std::vector<Node>> nodes;
    int NumberOfNodes;
public:
    Nodes(std::vector<std::vector<double>> &ShapeValues,std::vector<std::vector<int>> &RegionOfWork)
    {
        std::vector<Node> tempNode;
        for(int i = 0; i < ShapeValues.size();i++)
        {
            for(int j = 0; j < ShapeValues[i].size();j++)
            {
                if (RegionOfWork[i][j] == 2)
                {
                    tempNode.push_back(Node(ShapeValues[i][j],(double)j,(double)-i,true));
                }
                else if(RegionOfWork[i][j] == 1)
                {
                    tempNode.push_back(Node(ShapeValues[i][j],(double)j,(double)-i,false));
                }
            }
            if (tempNode.size() > 0)
            {
                nodes.push_back(tempNode);
                tempNode.clear();
            }
        }

        int iprev = -1, jprev = -1;
        int i = 0, j = 0;
        int NodeNumber = 0;
        nodes[0][0].setNodeNumber(NodeNumber);
        NodeNumber++;
        while(i!= iprev || j != jprev)
        {
            iprev = i;
            jprev = j;
            NearBoundaryNodes(i,j,NodeNumber);
        }
        for (i = 0; i < nodes.size(); i++)
        {
            for (j = 0; j < nodes[i].size(); j++)
            {
                if (!nodes[i][j].get_isNumbered())
                {
                    nodes[i][j].setNodeNumber(NodeNumber);
                    NodeNumber++;
                }
            }
        }
        NumberOfNodes = NodeNumber;
    }
    void NearBoundaryNodes(int &i, int &j,int &NodeNumber)
    {
        bool isUpperExist = false;
        bool isLowerExist = false;
        bool isLeftExist = false;
        bool isRightExist = false;
        bool isUpperLeftExist = false;
        bool isUpperRightExist = false;
        bool isLowerRightExist = false;
        bool isLowerLeftExist = false;

        if (i - 1 >= 0)
            isUpperExist = true;
        if (i + 1 < nodes.size())
            isLowerExist = true;
        if (j - 1 >= 0)
            isLeftExist = true;
        if (j + 1 < nodes[i].size())
            isRightExist = true;
        if (i - 1 >= 0 && j - 1 >= 0)
            isUpperLeftExist = true;
        if (i - 1 >= 0 && j + 1 < nodes[i-1].size())
            isUpperRightExist = true;
        if (i + 1 < nodes.size() && j - 1 >= 0)
            isLowerLeftExist = true;
        if (i + 1 < nodes.size() && j + 1 < nodes[i+1].size())
            isLowerRightExist = true;

        if(isUpperRightExist)
        {
            if (nodes[i - 1][j + 1].get_isBoundary() && !nodes[i - 1][j + 1].get_isNumbered())
            {
                i = i - 1;
                j = j + 1;
                nodes[i][j].setNodeNumber(NodeNumber);
                NodeNumber++;
                return;
            }
        }
        if(isRightExist)
        {
            if (nodes[i][j + 1].get_isBoundary() && !nodes[i][j + 1].get_isNumbered())
            {
                j = j + 1;
                nodes[i][j].setNodeNumber(NodeNumber);
                NodeNumber++;
                return;
            }
        }
        if(isLowerRightExist)
        {
            if (nodes[i + 1][j + 1].get_isBoundary() && !nodes[i + 1][j + 1].get_isNumbered())
            {
                i = i + 1;
                j = j + 1;
                nodes[i][j].setNodeNumber(NodeNumber);
                NodeNumber++;
                return;
            }
        }
        if(isLowerExist)
        {
            if (nodes[i + 1][j].get_isBoundary() && !nodes[i + 1][j].get_isNumbered())
            {
                i = i + 1;
                nodes[i][j].setNodeNumber(NodeNumber);
                NodeNumber++;
                return;
            }
        }
        if(isLowerLeftExist)
        {
            if (nodes[i + 1][j - 1].get_isBoundary() && !nodes[i + 1][j - 1].get_isNumbered())
            {
                i = i + 1;
                j = j - 1;
                nodes[i][j].setNodeNumber(NodeNumber);
                NodeNumber++;
                return;
            }
        }
        if(isLeftExist)
        {
            if (nodes[i][j - 1].get_isBoundary() && !nodes[i][j - 1].get_isNumbered())
            {
                j = j - 1;
                nodes[i][j].setNodeNumber(NodeNumber);
                NodeNumber++;
                return;
            }
        }
        if(isUpperLeftExist)
        {
            if (nodes[i - 1][j - 1].get_isBoundary() && !nodes[i - 1][j - 1].get_isNumbered())
            {
                i = i - 1;
                j = j - 1;
                nodes[i][j].setNodeNumber(NodeNumber);
                NodeNumber++;
                return;
            }
        }
        if(isUpperExist)
        {
            if (nodes[i - 1][j].get_isBoundary() && !nodes[i - 1][j].get_isNumbered())
            {
                i = i - 1;
                nodes[i][j].setNodeNumber(NodeNumber);
                NodeNumber++;
                return;
            }
        }
        else
        {
            return;
        }
    }
    void printNodeNumbers() const
    {
        for(auto row: nodes)
        {
            for(auto &n: row)
            {
                std::cout << n.getNodeNumber() << "\t";
            }
            std::cout << std::endl;
        }
    }
    Node getNodeByNumber(int _NodeNumber)
    {
        for(int i = 0; i < nodes.size(); i++)
        {
            for(int j = 0; j < nodes[i].size(); j++)
            {
                if(nodes[i][j].getNodeNumber() == _NodeNumber )
                {
                    return nodes[i][j];
                }
            }
        }
        return Node(0.0,-1,0);
    }
};

struct Element
{
    int NodeNumbers[3];
    int ElementNumber;
    Element(int _ElementNumber,int _NodeNumber1, int _NodeNumber2, int _NodeNumber3)
    {
        NodeNumbers[0] = _NodeNumber1;
        NodeNumbers[1] = _NodeNumber2;
        NodeNumbers[2] = _NodeNumber3;
        ElementNumber = _ElementNumber;
    }
    void PrintElement() const
    {
        std::cout << "ElementNo: " << ElementNumber << ", Nodes: " << NodeNumbers[0] << ", " << NodeNumbers[1] << ", " << NodeNumbers[2] << std::endl;
    }
};



class FiniteElement : public Nodes
{
private:
    int NumberOfElements;
    bool isSolved;
    std::vector<std::vector<double>> K;
    std::vector<double> X;
    std::vector<double> b;
    std::vector<std::vector<Element>> Elements;
public:
    FiniteElement(std::vector<std::vector<double>> &ShapeValues, std::vector<std::vector<int>> &RegionOfWork) 
        : Nodes(ShapeValues,RegionOfWork)
    {
        isSolved = false;
        int ElementNumber = 0;
        std::vector<Element> tempElement;
        for (int i = 0; i < nodes.size(); i++)
        {
            bool isbehind = false;
            bool isforward = false;
            for (int j = 0; j < nodes[i].size(); j++)
            {
                if(i+1 < nodes.size())
                {
                    int k;    
                    if (j == 0)
                    {
                        k = 0;
                        while(k<nodes[i+1].size())
                        {
                            if ((int)(nodes[i][j].get_xPosition()-nodes[i+1][k].get_xPosition()) == 1)
                            {
                                isbehind = true;
                                isforward = false;
                                break;
                            }
                            else if((int)(nodes[i][j].get_xPosition()-nodes[i+1][k].get_xPosition()) == -1)
                            {
                                isbehind = false;
                                isforward = true;
                                break;
                            }
                            else if((int)(nodes[i][j].get_xPosition()-nodes[i+1][k].get_xPosition()) == 0)
                            {
                                isbehind = false;
                                isforward = false;
                                break;
                            }
                            k++;
                        }                   
                    }
                    if(isbehind)
                    {
                        if (j == 0)
                        {
                            if(j+k+1 < nodes[i+1].size())
                            {
                                tempElement.push_back(Element(ElementNumber, nodes[i][j].getNodeNumber(), nodes[i + 1][j + k].getNodeNumber(), nodes[i + 1][j + k + 1].getNodeNumber()));
                                ElementNumber++;
                            }
                            if(j+k+2 < nodes[i+1].size())
                            {
                                tempElement.push_back(Element(ElementNumber, nodes[i][j].getNodeNumber(), nodes[i + 1][j + k + 1].getNodeNumber(), nodes[i + 1][j + k + 2].getNodeNumber()));
                                ElementNumber++;
                                if (j + 1 < nodes[i].size())
                                {
                                    tempElement.push_back(Element(ElementNumber, nodes[i][j].getNodeNumber(), nodes[i + 1][j + k + 2].getNodeNumber(), nodes[i][j + 1].getNodeNumber()));
                                    ElementNumber++;
                                }
                            }
                        }
                        else if(j+k+2 < nodes[i+1].size())
                        {
                            tempElement.push_back(Element(ElementNumber, nodes[i][j].getNodeNumber(), nodes[i + 1][j + k + 1].getNodeNumber(), nodes[i + 1][j + k + 2].getNodeNumber()));
                            ElementNumber++;
                            if (j + 1 < nodes[i].size())
                                {
                                    tempElement.push_back(Element(ElementNumber, nodes[i][j].getNodeNumber(), nodes[i + 1][j + k + 2].getNodeNumber(), nodes[i][j + 1].getNodeNumber()));
                                    ElementNumber++;
                                }
                        }
                        else if(j+1 == nodes[i].size()-1)
                        {
                            tempElement.push_back(Element(ElementNumber, nodes[i][j].getNodeNumber(), nodes[i + 1][j + k + 1].getNodeNumber(), nodes[i][j + 1].getNodeNumber()));
                            ElementNumber++;
                        }
                    }
                    else if(isforward)
                    {
                        if(j == 0)
                        {
                            if (j + 1 < nodes[i].size())
                            {
                                tempElement.push_back(Element(ElementNumber, nodes[i][j].getNodeNumber(), nodes[i + 1][j + k].getNodeNumber(), nodes[i][j + 1].getNodeNumber()));
                                ElementNumber++;
                            }
                        }
                        else if(j + k < nodes[i+1].size())
                        {
                            tempElement.push_back(Element(ElementNumber, nodes[i][j].getNodeNumber(), nodes[i + 1][j + k - 1].getNodeNumber(), nodes[i + 1][j + k].getNodeNumber()));
                            ElementNumber++;
                            if (j + 1 < nodes[i].size())
                                {
                                    tempElement.push_back(Element(ElementNumber, nodes[i][j].getNodeNumber(), nodes[i + 1][j + k].getNodeNumber(), nodes[i][j + 1].getNodeNumber()));
                                    ElementNumber++;
                                }
                        }
                        else if(j+1 == nodes[i].size()-1)
                        {
                            tempElement.push_back(Element(ElementNumber, nodes[i][j].getNodeNumber(), nodes[i + 1][j + k - 1 ].getNodeNumber(), nodes[i][j + 1].getNodeNumber()));
                            ElementNumber++;
                        }
                    }
                    else
                    {
                        if (j + 1 < nodes[i + 1].size())
                        {
                            tempElement.push_back(Element(ElementNumber, nodes[i][j].getNodeNumber(), nodes[i + 1][j].getNodeNumber(), nodes[i + 1][j + 1].getNodeNumber()));
                            ElementNumber++;
                            if (j + 1 < nodes[i].size())
                            {
                                tempElement.push_back(Element(ElementNumber, nodes[i][j].getNodeNumber(), nodes[i + 1][j + 1].getNodeNumber(), nodes[i][j + 1].getNodeNumber()));
                                ElementNumber++;
                            }
                        }
                        else if(j + 1 == nodes[i].size()-1)
                        {
                            tempElement.push_back(Element(ElementNumber, nodes[i][j].getNodeNumber(), nodes[i + 1][j].getNodeNumber(), nodes[i][j + 1].getNodeNumber()));
                            ElementNumber++;
                        }
                    }
                }
            }
            Elements.push_back(tempElement);
            tempElement.clear();
        }
        NumberOfElements = ElementNumber;

        K.resize(NumberOfNodes,std::vector<double>(NumberOfNodes,0.0));
        b.resize(NumberOfNodes,0.0);
        double be[3], ce[3], deltae, Ke[3][3];
        for (auto rows: Elements)
        {
            for (const auto &e: rows)
            {
                be[0] = getNodeByNumber(e.NodeNumbers[1]).get_yPosition() - getNodeByNumber(e.NodeNumbers[2]).get_yPosition();
                be[1] = getNodeByNumber(e.NodeNumbers[2]).get_yPosition() - getNodeByNumber(e.NodeNumbers[0]).get_yPosition();
                be[2] = getNodeByNumber(e.NodeNumbers[0]).get_yPosition() - getNodeByNumber(e.NodeNumbers[1]).get_yPosition();
                ce[0] = getNodeByNumber(e.NodeNumbers[2]).get_xPosition() - getNodeByNumber(e.NodeNumbers[1]).get_xPosition();
                ce[1] = getNodeByNumber(e.NodeNumbers[0]).get_xPosition() - getNodeByNumber(e.NodeNumbers[2]).get_xPosition();
                ce[2] = getNodeByNumber(e.NodeNumbers[1]).get_xPosition() - getNodeByNumber(e.NodeNumbers[0]).get_xPosition();

                deltae=0.5*(be[0]*ce[1]-be[1]*ce[0]);
                for (int i=0;i<3;i++) {
                    for (int j=0;j<3;j++) {
                        Ke[i][j] = (be[i]*be[j]+ce[i]*ce[j])/(4.0*deltae);
                    }
                }
                for (int i=0;i<3;i++) {
                    for (int j=0;j<3;j++) {
                        K[e.NodeNumbers[i]][e.NodeNumbers[j]] += Ke[i][j];
                    }
                }
            }
        }
        for (auto rows: nodes)
        {
            for (const auto &n: rows)
            {
                if (n.get_isBoundary())
                {
                    b[n.getNodeNumber()] = n.getNodeVoltage();
                    K[n.getNodeNumber()][n.getNodeNumber()] = 1.0;
                    for (int j = 0; j < NumberOfNodes; j++)
                    {
                        if (j != n.getNodeNumber())
                        {
                            b[j] -= K[j][n.getNodeNumber()]*n.getNodeVoltage();
                            K[n.getNodeNumber()][j] = 0.0;
                            K[j][n.getNodeNumber()] = 0.0;
                        }
                    }
                }
            }
        }

    }
    void PrintElements() const
    {
        for (auto rows: Elements)
        {
            for (auto &e: rows)
            {
                e.PrintElement();
            }
        }
    }
    void PrintK_Matrix()const
    {
        for(int i = 0; i < NumberOfNodes; i++)
        {
            for (int j = 0; j < NumberOfNodes; j++)
            {
                std::cout << "K[" << i << "][" << j << "] = " << K[i][j] << std::endl;
            }
        }
    }
    void PrintK_Diagonal_Matrix()const
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

    void Solve(double Tolerance)
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
                std::cout << "Convergence achieved, iteration: " << iter << std::endl;
                isSolved = true;
                setNodeVoltages();
                return;
            }
            else
            {
                if(iter == ITMAX)
                {
                    std::cout << "ITMAX exceeded. No convergence." << std::endl;
                    exit(0);
                }
            }

        }

    }
    double vNorm(const std::vector<double> &Z, const int N)
    {
        double result;
        result = 0.0;
        for (int i =0; i < N; i++)
        {
            result += pow(Z[i],2);
        }
        return result;
    }
    void prod(const std::vector<std::vector<double>> &_K, const std::vector<double> &_X, std::vector<double> &AX, const int control, const int N)
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
    void setNodeVoltages()
    {
        for (int i = 0; i < nodes.size(); i++)
        {
            for (int j = 0; j < nodes[i].size(); j++)
            {
                nodes[i][j].setNodeVoltage(X[nodes[i][j].getNodeNumber()]);
            }
        }
    }
    void PrintResults() const
    {
        if (isSolved)
        {
            for (int i = 0; i < nodes.size(); i++)
            {
                for (int j = 0; j < nodes[i].size(); j++)
                {
                    printf("%.3lf\t",nodes[i][j].getNodeVoltage());
                }
                printf("\n");
            }
        }
        std::cout << "\nNodes Voltages According to their node numbers:\n";
        for (int i = 0; i < X.size(); i++)
        {
            std::cout << "X[" << i << "]: " << X[i] << std::endl;
        }
    }
};

void set_initial_conditions(std::vector<std::vector<double>> &grid, std::vector<std::vector<int>> &RegionOfWork)
{
    for (int i = 0; i < grid.size(); i++)
    {
        for (int j = 0; j < grid[0].size(); j++)
        {
            if (j == 0 && i != grid.size() - 1) // setting the "a" line 
            {
                grid[i][j] = 20.0;
                RegionOfWork[i][j] = 2;
            }
            else if (i == j) // setting the "c" line
            {
                grid[i][j] = 10.0;
                RegionOfWork[i][j] = 2;
            }
            else if(i == grid.size()-1)
            {
                grid[i][j] = 0.0;
                RegionOfWork[i][j] = 2;
            }
            else if(i > j)
            {
                RegionOfWork[i][j] = 1;
            }
        }
    }
}

void PrintShapeValues(const std::vector<std::vector<double>> &grid)
{
    for (int i = 0; i < grid.size(); i++)
    {
        for (int j = 0; j < grid[0].size(); j++)
        {
            if (i < j) // Print needed elements
            {
                continue;
            }
            printf("%.3lf\t", grid[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
void PrintRegionOfWork(const std::vector<std::vector<int>> &grid)
{
    for (int i = 0; i < grid.size(); i++)
    {
        for (int j = 0; j < grid[0].size(); j++)
        {
            if (i < j) // Print needed elements
            {
                continue;
            }
            printf("%d ", grid[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main()
{
    int rows = 13, columns = 13;
    std::vector<std::vector<double>> ShapeValues;
    std::vector<std::vector<int>> RegionOfWork;

    ShapeValues.resize(rows,std::vector<double>(columns,0));
    RegionOfWork.resize(rows,std::vector<int>(columns,0));

    set_initial_conditions(ShapeValues,RegionOfWork);
    // PrintShapeValues(ShapeValues);
    // PrintRegionOfWork(RegionOfWork);
    
    FiniteElement e(ShapeValues,RegionOfWork);
    e.printNodeNumbers();
    e.Solve(1E-6);
    e.PrintResults();
    // e.PrintElements();
    // e.PrintK_Matrix();
    // e.PrintK_Diagonal_Matrix();
}
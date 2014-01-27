#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

class Node {
public:
    double state;
    int count;
    
    Node(double state, int count);
    virtual ~Node();
private:
};

class Distribution {
public:
    Distribution(int maxNodes);
    virtual ~Distribution();   
    
    void addNode(double state, int count);    
    void clear();
    
    void update(double* states, int length);
    
    int getNumNodes() const;
    double getState(int nodeId) const;
    int getCount(int nodeId) const;
private:
    int numNodes;
    int maxNodes;
    Node** nodes;
};
#endif
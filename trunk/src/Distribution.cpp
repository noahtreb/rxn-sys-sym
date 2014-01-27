#include "Distribution.h"
#include <algorithm>

using namespace std;

Node::Node(double state, int count) {
    this->state = state;
    this->count = count;
}

Node::~Node() {
    
}

Distribution::Distribution(int maxNodes) {    
    this->nodes = new Node*[maxNodes];
    
    for (int i = 0; i < maxNodes; i++) {
        this->nodes[i] = new Node(-1, 0);
    }
    
    this->maxNodes = maxNodes;
    this->numNodes = 0;
}

Distribution::~Distribution() {
    for (int i = 0; i < this->maxNodes; i++) {
        delete this->nodes[i];
    }
    
    delete[] this->nodes;
}

void Distribution::addNode(double state, int count) {
    this->nodes[this->numNodes]->state = state;
    this->nodes[this->numNodes]->count = count;
    this->numNodes++;
}

void Distribution::clear() {
    this->numNodes = 0;
}

void Distribution::update(double* states, int length) {
    this->numNodes = 0;
    
    std::sort(states.begin(), states.end());
    
    int count = 0;
    double lastState = states[0];
    for (int i = 0; i < length; i++) {
        if (states[i] == lastState) {
            count++;
        } else {
            this->addNode(lastState, count);
            
            lastState = states[i];
            count = 1;
        }
    }
}

int Distribution::getNumNodes const {
    return this->numNodes;
}

double Distribution::getState(int nodeId) const {
    return this->nodes[nodeId]->state;
}

int Distribution::getCount(int nodeId) const {
    return this->nodes[nodeId]->count;
}
#include "System.h"
#include "PriorityQueue.h"
#include "Reaction.h"
#include "Species.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <random>

using namespace std;

System::System(double vol, const int numRxns, Reaction** rxns, const int numSpecies, 
        Species** species, double* initSpeciesState, bool* speciesStateChanges) : 
        numRxns(numRxns), numSpecies(numSpecies) {    
    this->rxns = rxns; 
    this->rxnPq = new PriorityQueue(numRxns, true);
    
    this->vol = vol;
    this->volRatio = 1;
    
    this->time = 0;
    
    this->species = species;
    this->speciesState = new double[numSpecies];
    this->speciesStateChanges = new bool[numSpecies];
    for (int i = 0; i < numSpecies; i++) {
        this->speciesState[i] = initSpeciesState[i];
        this->speciesStateChanges[i] = speciesStateChanges[i];
        this->species[i]->state = &speciesState[i];
    }
    
    for (int i = 0; i < numRxns; i++) {
        this->rxns[i]->updateProp(this->speciesState, this->volRatio);
        this->rxns[i]->oldProp = this->rxns[i]->prop;
    }
}

System::System(const System& other) : numRxns(other.numRxns), numSpecies(other.numSpecies) {
    // Does not copy the Mersenne twister.
    
    this->rxns = new Reaction*[this->numRxns];
    for (int i = 0; i < this->numRxns; i++) {
        this->rxns[i] = new Reaction(*(other.rxns[i]));
    }
    
    this->rxnPq = new PriorityQueue(*(other.rxnPq));
    
    this->vol = other.vol;
    this->volRatio = other.volRatio;
    
    this->time = other.time;
    
    this->species = new Species*[this->numSpecies];    
    this->speciesState = new double[this->numSpecies];
    this->speciesStateChanges = new bool[this->numSpecies];
    for (int i = 0; i < numSpecies; i++) {
        this->species[i] = new Species(*(other.species[i]));
        this->speciesState[i] = other.speciesState[i];
        this->speciesStateChanges[i] = other.speciesStateChanges[i];
    }
}

System::~System() {
    for (int i = 1; i < this->numRxns; i++) {
        delete this->rxns[i];
    }
    
    delete[] this->rxns;
    delete this->rxnPq;
    delete this->speciesState;
    delete this->rng;
}

void System::init(int seed) {
    this->rng = new std::mt19937(seed);    
    initFwd(); 
}

void System::initFwd() {
    this->rxnPq->clearQueue();
    
    for (int i = 0; i < this->numRxns; i++) {        
        if (this->rxns[i]->prop > 0) {
            this->rxnPq->addNode(this->rxns[i], - log(1.0 * (*this->rng)() / this->rng->max()) / this->rxns[i]->prop + this->time);
        } else {
            this->rxnPq->addNode(this->rxns[i], DBL_MAX);
        }
    }
}

void System::initRev() {
    this->rxnPq->clearQueue();
    
    for (int i = 0; i < this->numRxns; i++) {        
        if (this->rxns[i]->prop > 0) {
            this->rxnPq->addNode(this->rxns[i], log(1.0 * (*this->rng)() / this->rng->max()) / this->rxns[i]->prop + this->time);
        } else {
            this->rxnPq->addNode(this->rxns[i], -DBL_MAX);
        }
    }    
}

void System::setRxnTimes(Reaction* execRxn, double execRxnTime, int dir) {     
    bool updateMask[this->numRxns];    
    for (int i = 0; i < this->numRxns; i++) {
        updateMask[i] = false;
    }
    
    int rxnId;
    for (int i = 0; i < execRxn->numDeps; i++) {
        rxnId = execRxn->deps[i];
        updateMask[rxnId] = true;        
        this->rxns[rxnId]->updateProp(this->speciesState, this->volRatio);
    }
    
    execRxn->updateProp(this->speciesState, this->volRatio);
    if (execRxn->prop > 0) {
        this->rxnPq->updateNodeByNodeId(0, execRxnTime - dir * log(1.0 * (*this->rng)() / this->rng->max()) / execRxn->prop);
    }
    
    Reaction* rxn;
    int numNodes = this->rxnPq->getNumNodes();
    for (int i = 0; i < numNodes; i++) {
        rxn = this->rxnPq->getRxnByNodeId(i);
        if (updateMask[rxn->id]) {
            if (rxn->prop > 0) {
                this->rxnPq->updateNodeByNodeId(i, rxn->oldProp / rxn->prop * (this->rxnPq->getTimeByNodeId(i) - execRxnTime) + execRxnTime);                 
            } else {
                this->rxnPq->updateNodeByNodeId(i, dir * DBL_MAX);
            }
        }
    }
    
    this->rxnPq->heapSort(); 
}

void System::execRxn(int dir) {
    Reaction* rxn = this->rxnPq->getNextRxn();
    double time = this->rxnPq->getNextTime();
    
    bool negState = false;
    
    for (int i = 0; i < rxn->numStoichSpecies; i++) {
        if (this->speciesStateChanges[rxn->stoichSpeciesIds[i]]) {
            this->speciesState[rxn->stoichSpeciesIds[i]] += dir * rxn->stoichCoeffs[i];

            if (this->speciesState[rxn->stoichSpeciesIds[i]] < 0) {
                negState = true;
            }
        }
    }
    
    if (!negState) {
        this->setRxnTimes(rxn, time, dir);
        this->updateTime(time);
    } else {
        for (int i = 0; i < rxn->numStoichSpecies; i++) {
            if (this->speciesStateChanges[rxn->stoichSpeciesIds[i]]) {
                this->speciesState[rxn->stoichSpeciesIds[i]] -= dir * rxn->stoichCoeffs[i];
            }
        }
        
        this->rxnPq->updatePqByNodeId(0, dir * DBL_MAX);
    }
}

void System::updateTime(double newTime) {
    this->time = newTime;
}

/*
int System::findClosestTimePt(double time) const {
    int index = floor((time-this->startTime)/this->timeStep);
    
    if ((time-this->timePts[index])*(time-this->timePts[index]) < (time-this->timePts[index+1])*(time-this->timePts[index+1])) {
        return index;
    } else {
        return index + 1;
    }
}*/
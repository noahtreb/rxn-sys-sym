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
        Species** species, double* initSpeciesState) : 
        numRxns(numRxns), numSpecies(numSpecies) {    
    this->rxns = rxns; 
    this->rxnPq = new PriorityQueue(numRxns, true);
    
    this->vol = vol;
    this->volRatio = 1;
    
    this->time = 0;
    
    this->species = species;
    this->speciesState = new double[numSpecies];
    for (int i = 0; i < numSpecies; i++) {
        this->speciesState[i] = initSpeciesState[i];
        this->species[i]->statePtr = &speciesState[i];
    }
    
    for (int i = 0; i < numRxns; i++) {
        this->rxns[i]->updateProp(this->speciesState, this->volRatio);
        this->rxns[i]->oldProp = this->rxns[i]->prop;
    }
}

System::System(const System& other) : numRxns(other.numRxns), numSpecies(other.numSpecies) {
    // Does not copy the Mersenne twister.    
    
    this->species = new Species*[this->numSpecies];    
    this->speciesState = new double[this->numSpecies];
    for (int i = 0; i < numSpecies; i++) {
        this->species[i] = new Species(*(other.species[i]));
        this->speciesState[i] = other.speciesState[i];
        this->species[i]->statePtr = &this->speciesState[i];
    }
    
    this->rxns = new Reaction*[this->numRxns];
    for (int i = 0; i < this->numRxns; i++) {
        this->rxns[i] = new Reaction(*(other.rxns[i]));
        
        for (int j = 0; j < this->rxns[i]->numStoichSpecies; j++) {
            this->rxns[i]->stoichSpecies[j] = this->species[this->rxns[i]->stoichSpecies[j]->id];
        }
        
        for (int j = 0; j < this->rxns[i]->numRateSpecies; j++) {
            this->rxns[i]->rateSpecies[j] = this->species[this->rxns[i]->rateSpecies[j]->id];
        }
    }    
    
    this->rxnPq = new PriorityQueue(*(other.rxnPq));    
    for (int i = 0; i < this->rxnPq->getNumNodes(); i++) {
        int rxnId = this->rxnPq->getRxnIdByNodeId(i);
        this->rxnPq->updateRxnByNodeId(i, this->rxns[rxnId]);
    }
    
    this->vol = other.vol;
    this->volRatio = other.volRatio;
    
    this->time = other.time;
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

void System::setRxnTimes(Reaction* execRxn, double execRxnTime, bool fwd) {     
    int dir;
    if (fwd) {
        dir = 1;
    } else {
        dir = -1;
    }
    
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

void System::execRxn(bool fwd) {
    Reaction* rxn = this->rxnPq->getNextRxn();
    double time = this->rxnPq->getNextTime();
    
    int dir;
    if (fwd) {
        dir = 1;
    } else {
        dir = -1;
    }
    
    bool negState = false;
    bool boundBreach = false;
    
    for (int i = 0; i < rxn->numStoichSpecies; i++) {
        if (rxn->stoichSpecies[i]->stateChanges) {
            rxn->stoichSpecies[i]->statePtr[0] += dir * rxn->stoichCoeffs[i];

            if (rxn->stoichSpecies[i]->statePtr[0] < 0) {
                negState = true;
            }
            
            if ((fwd && rxn->stoichSpecies[i]->stateBoundedFwd) || (!fwd && rxn->stoichSpecies[i]->stateBoundedRev)) {
                Species* species = rxn->stoichSpecies[i];
                if (species->statePtr[0] < species->stateLowerBound || species->statePtr[0] > species->stateUpperBound) {
                    boundBreach = true;
                }
            }
        }
    }
    
    if (!negState && !boundBreach) {
        this->setRxnTimes(rxn, time, fwd);
        this->updateTime(time);
    //} else if (boundBreach) {
        //***TODO: implement behavior
    } else {
        for (int i = 0; i < rxn->numStoichSpecies; i++) {
            if (rxn->stoichSpecies[i]->stateChanges) {
                rxn->stoichSpecies[i]->statePtr[0] -= dir * rxn->stoichCoeffs[i];
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
#ifndef SYSTEM_H
#define SYSTEM_H

#include <random>

class Reaction;
class PriorityQueue;

class System{
public:
    const int numRxns;
    const int numSpecies;
    bool* speciesStateChanges;
    
    const int numTimePts;
    const double startTime;
    const double endTime;
    const double timeStep;
    double* timePts;
    
    double** avFwdSpeciesState;
    double** avRevSpeciesState;
    
    double* speciesState;
    double vol;
    double volRatio;
    
    double time;    
    
    Reaction** rxns;
    PriorityQueue* rxnPq;
    
    std::mt19937* rng;
    
    System(double vol, const int numRxns, Reaction** rxns, const int numSpecies, double* initSpeciesState, bool* speciesStateChanges, const double startTime, const double endTime, const double timeStep, const int numTimePts, double* timePts, double** avFwdSpeciesState, double** avRevSpeciesState);
    System(const System& other);
    virtual ~System();
    
    void init(int seed);
    void initFwd();
    void initRev();
    
    void setRxnTimes(Reaction* execRxn, double execRxnTime, int dir);
    void execRxn(int dir);
    void updateTime(double newTime);
private:
    int findClosestTimePt(double time) const;
};

#endif
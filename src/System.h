#ifndef SYSTEM_H
#define SYSTEM_H

#include <random>

class Species;
class Reaction;
class PriorityQueue;

class System{
public:
    const int numRxns;
    const int numSpecies;
    
    double* speciesState;
    double vol;
    double volRatio;
    
    double time;    
    
    Species** species;
    Reaction** rxns;
    PriorityQueue* rxnPq;
    
    std::mt19937* rng;
    
    System(double vol, const int numRxns, Reaction** rxns, const int numSpecies, Species** species, double* initSpeciesState);
    System(const System& other);
    virtual ~System();
    
    void init(int seed);
    void initFwd();
    void initRev();
    
    void setRxnTimes(Reaction* execRxn, double execRxnTime, bool fwd);
    void execRxn(bool fwd);
    void updateTime(double newTime);
private:
    //int findClosestTimePt(double time) const;
};

#endif
#ifndef REACTION_H
#define REACTION_H

class Reaction{
public:
    const int id;
    
    const int numStoichSpecies;
    const int* const stoichSpeciesIds;
    const int* const stoichCoeffs;
    
    const int rateLaw;
    const int numRateConsts;
    
    const int numRateSpecies;
    const int* const rateSpeciesIds;
    
    const int numDeps;
    const int* const deps;
    
    double prop;
    double oldProp;
    
    Reaction(const int id, const int numStoichSpecies, const int* const stoichSpeciesIds, const int* const stoichCoeffs, const int rateLaw, const int numRateConsts, double* rateConsts, const int numRateSpecies, const int* const rateSpeciesIds, const int numDeps, const int* const deps, double vol);
    Reaction(const Reaction& other);
    virtual ~Reaction();
    
    void updateProp(double* speciesState, double volRatio);
    
    void print() const;
private:    
    double* rateConsts;
};

#endif
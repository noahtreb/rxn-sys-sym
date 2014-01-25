#include "Reaction.h"
#include <stdio.h>

using namespace std;

Reaction::Reaction(const int id,
        const int numStoichSpecies, 
        const int* const stoichSpeciesIds, 
        const int* const stoichCoeffs, 
        const int rateLaw, 
        const int numRateConsts, 
        double* rateConsts, 
        const int numRateSpecies, 
        const int* const rateSpeciesIds, 
        const int numDeps,
        const int* const deps,
        double vol) : 
        id(id),
        numStoichSpecies(numStoichSpecies),
        stoichSpeciesIds(stoichSpeciesIds), 
        stoichCoeffs(stoichCoeffs), 
        rateLaw(rateLaw), 
        numRateConsts(numRateConsts),
        rateConsts(rateConsts),
        numRateSpecies(numRateSpecies), 
        rateSpeciesIds(rateSpeciesIds),
        numDeps(numDeps),
        deps(deps) {
    this->prop = 0;   
    this->oldProp = 0;
    
    switch(rateLaw) { // Convert to mesoscopic units.
        case 3: // Second-order (bimolecular) mass action            
            //this->rateConsts[0] = this->rateConsts[0] / vol / 6.022e23;
            break;
    }
}

Reaction::Reaction(const Reaction& other) : 
        id(other.id),
        numStoichSpecies(other.numStoichSpecies),
        stoichSpeciesIds(other.stoichSpeciesIds), 
        stoichCoeffs(other.stoichCoeffs), 
        rateLaw(other.rateLaw), 
        numRateConsts(other.numRateConsts),
        rateConsts(other.rateConsts),
        numRateSpecies(other.numRateSpecies), 
        rateSpeciesIds(other.rateSpeciesIds),
        numDeps(other.numDeps),
        deps(other.deps) {
    this->prop = other.prop;
    this->oldProp = other.oldProp;
}

Reaction::~Reaction() {
    delete[] this->stoichSpeciesIds;
    delete[] this->stoichCoeffs;
    delete[] this->rateConsts;
    delete[] this->rateSpeciesIds;
}

void Reaction::updateProp(double* speciesState, double volRatio) {
    this->oldProp = this->prop;    
    this->prop = this->rateConsts[0];
    
    for (int i = 0; i < this->numStoichSpecies; i++) {
        if (this->stoichCoeffs[i] < 0) {
            for (int j = 0; j < -this->stoichCoeffs[i]; j++) {
                this->prop *= 1.0 * (speciesState[this->stoichSpeciesIds[i]] - j);// / (j + 1);
            }
        }
    }
    
    /*switch(this->rateLaw) {
        case 2: // First-order mass action
            this->prop = this->rateConsts[0] * speciesState[this->rateSpeciesIds[0]];
            break;
        case 3: // Second-order (bimolecular) mass action
            this->prop = this->rateConsts[0] * speciesState[this->rateSpeciesIds[0]] * speciesState[this->rateSpeciesIds[1]] / volRatio;
            break;
    }*/
}

void Reaction::print() const {
    int stoichCoeff;
    int stoichSpeciesId;
    
    for (int i = 0; i < this->numStoichSpecies; i++) {
        stoichCoeff = this->stoichCoeffs[i];
        stoichSpeciesId = this->stoichSpeciesIds[i];
        
        if (i == 0 && stoichCoeff > 0) {            
            fprintf(stdout, " -> ");
        }
        
        if (stoichCoeff < 0) {
            fprintf(stdout, "%i(%i)", -stoichCoeff, stoichSpeciesId);
        } else {            
            fprintf(stdout, "%i(%i)", stoichCoeff, stoichSpeciesId);
        }
        
        if (i != this->numStoichSpecies - 1) {
            if (stoichCoeff * this->stoichCoeffs[i + 1] > 0) {
                fprintf(stdout, " + ");
            } else {
                fprintf(stdout, " -> ");
            }
        }
        
        if (i == this->numStoichSpecies - 1 && stoichCoeff < 0) {
            fprintf(stdout, " -> ");
        }
    }
    
    fprintf(stdout, "\n    ID: %i, Rate Constant: %e, Propensity: %e, Old Propensity: %e\n", this->id, this->rateConsts[0], this->prop, this->oldProp);
}
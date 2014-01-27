#include "Species.h"
#include <stdlib.h>

using namespace std;

Species::Species(const int id,
        double state, 
        const bool stateChanges,
        const double stateUpperBound, 
        const double stateLowerBound, 
        const bool stateBoundedFwd, 
        const bool stateBoundedRev) : 
        id(id),
        stateChanges(stateChanges),
        stateUpperBound(stateUpperBound),
        stateLowerBound(stateLowerBound),
        stateBoundedFwd(stateBoundedFwd),
        stateBoundedRev(stateBoundedRev) {
    this->state = state;
}

Species::Species(const Species& other) :
        id(other.id),
        stateChanges(other.stateChanges),
        stateUpperBound(other.stateUpperBound),
        stateLowerBound(other.stateLowerBound),
        stateBoundedFwd(other.stateBoundedFwd),
        stateBoundedRev(other.stateBoundedRev) {
    this->state = other.state;
}

Species::~Species() {
    
}
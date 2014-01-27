#ifndef SPECIES_H
#define SPECIES_H

class Species {
public:
    const int id;
    double* state;
    const bool stateChanges;
    const double stateUpperBound;
    const double stateLowerBound;
    const bool stateBoundedFwd;
    const bool stateBoundedRev;
    
    Species(const int id, double* state, const bool stateChanges, const double stateUpperBound, const double stateLowerBound, const bool stateBoundedFwd, const bool stateBoundedRev);
    Species(const Species& other);
    virtual ~Species();
private:
};

#endif
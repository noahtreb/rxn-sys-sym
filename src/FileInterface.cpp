#include "FileInterface.h"
#include "Species.h"
#include "System.h"
#include "Reaction.h"
#include <stdio.h>
#include <stdlib.h>
#include <netcdfcpp.h>

using namespace std;

FileInterface::FileInterface(string fileName) {
    this->fileName = fileName;
}

FileInterface::~FileInterface() {
    
}

System* FileInterface::readFileData(int& numTrials, double& startTime, double& endTime, int& numTimePts, double& timeStep, double& stoppingTol, int& numDataSavePts, int*& dataSavePts) const {
    NcFile file(this->fileName.c_str(), NcFile::ReadOnly);    
    if (!file.is_valid()) {
        fprintf(stderr, "Error: %s could not be opened.\n", this->fileName.c_str());
        abort();
    }
    
    file.get_var("stoppingTol")->get(&stoppingTol, 1);
    
    file.get_var("numDataSavePts")->get(&numDataSavePts, 1);
    dataSavePts = new int[numDataSavePts];
    file.get_var("dataSavePts")->get(dataSavePts, 1, numDataSavePts);
    
    file.get_var("numTrials")->get(&numTrials, 1);
    
    file.get_var("tStart")->get(&startTime, 1);
    file.get_var("tEnd")->get(&endTime, 1);
    file.get_var("timePts")->get(&numTimePts, 1);    
    timeStep = (endTime - startTime) / (numTimePts - 1);  
        
    int numSpecies;
    file.get_var("numSpecies")->get(&numSpecies, 1); 
    
    int numRxns;
    file.get_var("numRxns")->get(&numRxns, 1);
    
    int maxStoichSpecies = (int) file.get_dim("maxRxnSpecies")->size();    
    int maxRxnRateConsts = (int) file.get_dim("maxRxnRateConsts")->size();    
    int maxRxnRateSpecies = (int) file.get_dim("maxRxnRateSpecies")->size();    
    int maxRxnDeps = (int) file.get_dim("maxRxnDeps")->size();
    
    double vol;
    file.get_var("initVol")->get(&vol, 1);
    
    double* speciesInitState = new double[numSpecies];
    file.get_var("speciesInitState")->get(speciesInitState, 1, numSpecies);  
        
    int* speciesStateChangesTemp = new int[numSpecies];
    file.get_var("speciesStateChanges")->get(speciesStateChangesTemp, 1, numSpecies);
    
    bool* speciesStateChanges = new bool[numSpecies];
    for (int i = 0; i < numSpecies; i++) {
        if (speciesStateChangesTemp[i] == 1) {
            speciesStateChanges[i] = true;
        } else {
            speciesStateChanges[i] = false;
        }
    }
    
    int* speciesStateBoundedFwdTemp = new int[numSpecies];
    file.get_var("speciesStateBoundedFwd")->get(speciesStateBoundedFwdTemp, 1, numSpecies);
    
    bool* speciesStateBoundedFwd = new bool[numSpecies];
    for (int i = 0; i < numSpecies; i++) {
        if (speciesStateBoundedFwdTemp[i] == 1) {
            speciesStateBoundedFwd[i] = true;
        } else {
            speciesStateBoundedFwd[i] = false;
        }
    }
    
    int* speciesStateBoundedRevTemp = new int[numSpecies];
    file.get_var("speciesStateBoundedRev")->get(speciesStateBoundedRevTemp, 1, numSpecies);
    
    bool* speciesStateBoundedRev = new bool[numSpecies];
    for (int i = 0; i < numSpecies; i++) {
        if (speciesStateBoundedRevTemp[i] == 1) {
            speciesStateBoundedRev[i] = true;
        } else {
            speciesStateBoundedRev[i] = false;
        }
    }
        
    double* speciesStateLowerBounds = new double[numSpecies];
    file.get_var("speciesStateLowerBounds")->get(speciesStateLowerBounds, 1, numSpecies);
    
    double* speciesStateUpperBounds = new double[numSpecies];
    file.get_var("speciesStateUpperBounds")->get(speciesStateUpperBounds, 1, numSpecies);
    
    int* numStoichSpecies = new int[numRxns];
    file.get_var("numRxnSpecies")->get(numStoichSpecies, 1, numRxns);
    
    int** rxnStoichSpeciesIds = new int*[numRxns];
    for (int i = 0; i < numRxns; i++) {
        rxnStoichSpeciesIds[i] = new int[maxStoichSpecies];
        file.get_var("rxnSpecies")->set_cur(0, i, 0);
        file.get_var("rxnSpecies")->get(rxnStoichSpeciesIds[i], 1, 1, maxStoichSpecies);
        
        for (int j = 0; j < maxStoichSpecies; j++) {
            rxnStoichSpeciesIds[i][j]--;
        }
    }
            
    int** rxnStoichCoeffs = new int*[numRxns];
    for (int i = 0; i < numRxns; i++) {
        rxnStoichCoeffs[i] = new int[maxStoichSpecies];
        file.get_var("rxnSpeciesCoeffs")->set_cur(0, i, 0);
        file.get_var("rxnSpeciesCoeffs")->get(rxnStoichCoeffs[i], 1, 1, maxStoichSpecies);
    }
    
    int* rxnRateLaws = new int[numRxns];
    file.get_var("rxnRateLaws")->get(rxnRateLaws, 1, numRxns);
    
    int* numRxnRateConsts = new int[numRxns];
    file.get_var("numRxnRateConsts")->get(numRxnRateConsts, 1, numRxns);
    
    double** rxnRateConsts = new double*[numRxns];
    for (int i = 0; i < numRxns; i++) {
        rxnRateConsts[i] = new double[maxRxnRateConsts];        
        file.get_var("rxnRateConsts")->set_cur(0, i, 0);
        file.get_var("rxnRateConsts")->get(rxnRateConsts[i], 1, 1, maxRxnRateConsts);
    }
        
    int* numRxnRateSpecies = new int[maxRxnRateSpecies];
    file.get_var("numRxnRateSpecies")->get(numRxnRateSpecies, 1, numRxns);
    
    int** rxnRateSpeciesIds = new int*[numRxns];
    for (int i = 0; i < numRxns; i++) {
        rxnRateSpeciesIds[i] = new int[maxRxnRateSpecies];
        file.get_var("rxnRateSpecies")->set_cur(0, i, 0);
        file.get_var("rxnRateSpecies")->get(rxnRateSpeciesIds[i], 1, 1, maxRxnRateSpecies);
        
        for (int j = 0; j < maxRxnRateSpecies; j++) {
            rxnRateSpeciesIds[i][j]--;
        }
    }
    
    int* numRxnDeps = new int[numRxns];
    file.get_var("numRxnDeps")->get(numRxnDeps, 1, numRxns);
    
    int** rxnDeps = new int*[numRxns];
    for (int i = 0; i < numRxns; i++) {
        rxnDeps[i] = new int[maxRxnDeps];
        file.get_var("rxnDeps")->set_cur(0, i, 0);
        file.get_var("rxnDeps")->get(rxnDeps[i], 1, 1, maxRxnDeps);
        
        for (int j = 0; j < maxRxnDeps; j++) {
            rxnDeps[i][j]--;
        }
    }
    
    Species** species = new Species*[numSpecies];
    for (int i = 0; i < numSpecies; i++) {
        species[i] = new Species(i, NULL, speciesStateChanges[i], speciesStateUpperBounds[i], speciesStateLowerBounds[i], speciesStateBoundedFwd[i], speciesStateBoundedRev[i]);
    }
    
    Reaction** rxns = new Reaction*[numRxns];    
    for (int i = 0; i < numRxns; i++) {        
        rxns[i] = new Reaction(i, species, numStoichSpecies[i], rxnStoichSpeciesIds[i], rxnStoichCoeffs[i], rxnRateLaws[i], numRxnRateConsts[i], rxnRateConsts[i], numRxnRateSpecies[i], rxnRateSpeciesIds[i], numRxnDeps[i], rxnDeps[i], vol);
    }
    
    System* sys = new System(vol, numRxns, rxns, numSpecies, species, speciesInitState);
    
    delete[] speciesInitState;
    delete[] speciesStateChangesTemp;
    delete[] speciesStateChanges;
    delete[] numStoichSpecies;    
    delete[] rxnStoichSpeciesIds;    
    delete[] rxnStoichCoeffs;    
    delete[] rxnRateLaws;    
    delete[] numRxnRateConsts;    
    delete[] rxnRateConsts;        
    delete[] numRxnRateSpecies;    
    delete[] rxnRateSpeciesIds;
    delete[] numRxnDeps;
    delete[] rxnDeps;
    
    return sys;
}

double** FileInterface::readLastDataPt(std::string varName, int numTrials, int numSpecies, int numTimePts) const {
    NcFile file(this->fileName.c_str(), NcFile::ReadOnly);
    if (!file.is_valid()) {
        fprintf(stderr, "Error: %s could not be opened.\n", this->fileName.c_str());
        abort();
    }     
    
    double** lastDataPt = new double*[numTrials];
    for (int i = 0; i < numTrials; i++) {
        lastDataPt[i] = new double[numSpecies];
    }
    
    NcVar* var = file.get_var(varName.c_str());    
    for (int i = 0; i < numTrials; i++) {
        var->set_cur(0, i, numTimePts-1, 0);
        var->get(lastDataPt[i], 1, 1, 1, numSpecies);
    }
    
    return lastDataPt;
}

double*** FileInterface::readStateData(std::string varName, int numTrials, int numSpecies, int numTimePts) const {
    NcFile file(this->fileName.c_str(), NcFile::ReadOnly);
    if (!file.is_valid()) {
        fprintf(stderr, "Error: %s could not be opened.\n", this->fileName.c_str());
        abort();
    } 
    
    double*** stateData = new double**[numTrials];
    for (int i = 0; i < numTrials; i++) {
        stateData[i] = new double*[numTimePts];
        for (int j = 0; j < numTimePts; j++) {
            stateData[i][j] = new double[numSpecies];
        }
    }
    
    NcVar* var = file.get_var(varName.c_str());    
    for (int i = 0; i < numTrials; i++) {
        for (int j = 0; j < numTimePts; j++) {
            var->set_cur(0, i, j, 0);
            var->get(stateData[i][j], 1, 1, 1, numSpecies);
        }
    }
    
    return stateData;
}

void FileInterface::overwriteStoppingTol(double stoppingTol) const {
    NcFile file(this->fileName.c_str(), NcFile::Write);
    if (!file.is_valid()) {
        fprintf(stderr, "Error: %s could not be opened.\n", this->fileName.c_str());
        abort();
    }
    
    file.get_var("stoppingTol")->put(&stoppingTol, 1);
}

void FileInterface::writeTimeData(double* time, int numTimePts) const {
    NcFile file(this->fileName.c_str(), NcFile::Write);
    if (!file.is_valid()) {
        fprintf(stderr, "Error: %s could not be opened.\n", this->fileName.c_str());
        abort();
    }
    
    file.get_var("time")->put(time, 1, numTimePts);
}

void FileInterface::writeStateData(string varName, int trial, double** data, int numSpecies, int numTimePts) const {
    NcFile file(this->fileName.c_str(), NcFile::Write);
    if (!file.is_valid()) {
        fprintf(stderr, "Error: %s could not be opened.\n", this->fileName.c_str());
        abort();
    }
    
    long dataCounts[] = {1, 1, numTimePts, 1};    
    for (int i = 0; i < numSpecies; i++) {
        long currPos[] = {0, trial, 0, i};
        file.get_var(varName.c_str())->set_cur(currPos);
        file.get_var(varName.c_str())->put(data[i], dataCounts);
    }
}

void FileInterface::writeStateData(string varName, int dataSavePtId, int trial, double** data, int numSpecies, int numTimePts) const {
    NcFile file(this->fileName.c_str(), NcFile::Write);
    if (!file.is_valid()) {
        fprintf(stderr, "Error: %s could not be opened.\n", this->fileName.c_str());
        abort();
    }
    
    long dataCounts[] = {1, 1, 1, numTimePts, 1};    
    for (int i = 0; i < numSpecies; i++) {
        long currPos[] = {0, dataSavePtId, trial, 0, i};
        file.get_var(varName.c_str())->set_cur(currPos);
        file.get_var(varName.c_str())->put(data[i], dataCounts);
    }
}

void FileInterface::writeAbsorbingCurrentData(int dataSavePtId, int absorbingCurrent) const {
    NcFile file(this->fileName.c_str(), NcFile::Write);
    if (!file.is_valid()) {
        fprintf(stderr, "Error: %s could not be opened.\n", this->fileName.c_str());
        abort();
    }
    
    file.get_var("absorbingCurrent")->set_cur(dataSavePtId);
    file.get_var("absorbingCurrent")->put(&absorbingCurrent, 1);
}
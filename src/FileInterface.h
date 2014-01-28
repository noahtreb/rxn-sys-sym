#ifndef FILEINTERFACE_H
#define FILEINTERFACE_H

#include <string>

class System;
class Reaction;

class FileInterface{
public:
    std::string fileName;
    
    FileInterface(std::string fileName);
    virtual ~FileInterface();
    
    System* readFileData(int& numTrials, double& startTime, double& endTime, int& numTimePts, double& timeStep, double& stoppingTol, int& numDataSavePts, int*& dataSavePts, int& numBoundedSpeciesStates) const;
    double** readInitDataPt(std::string varName, int numTrials, int numSpecies, int timePtId, int numTimePts) const;
    double** readDataPt(std::string varName, int dataSavePtId, int numTrials, int numSpecies, int timePtId, int numTimePts) const;
    double*** readStateData(std::string varName, int numTrials, int numSpecies, int numTimePts) const;
        
    void overwriteStoppingTol(double stoppingTol) const;
    void writeTimeData(double* time, int numTimePts) const;
    void writeInitStateData(std::string varName, int trial, double** data, int numSpecies, int numTimePts) const;    
    void writeStateData(std::string varName, int dataSavePtId, int trial, double** data, int numSpecies, int numTimePts) const;
    void writeAbsorbingCurrentData(int dataSavePtId, int absorbingCurrent) const;
private:
};

#endif
#ifndef MAIN_H
#define MAIN_H

#include "FileInterface.h"
#include "System.h"
#include <string>

void simFwd(System* sys, int numTimePts, double* time, double** state);
void simRev(System* sys, int numTimePts, double* time, double** state);
void writeStateData(System* sys, FileInterface* fi, std::string varName, int trial, double** state, int numTimePts, double* writeStart, double* writeEnd, omp_lock_t& lock);       
void averageStateData(int numTrials, int numTimePts, int numSpecies, double*** stateData, double** avStateData, double** lastAvStateData);
double calcDist(double** avStateData1, double** avStateData2);

#endif
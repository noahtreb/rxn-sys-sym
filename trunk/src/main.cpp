#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <float.h>
#include <random>
#include <string>
#include "main.h"
#include "FileInterface.h"
#include "PriorityQueue.h"
#include "System.h"
#include "Reaction.h"

using namespace std;

int main(int argc, const char* argv[]) {
    bool skipInitFwd = false;
    bool skipInitRev = false;
    bool skipRefineFwd = false;
    bool skipRefineRev = false;
    int seed;
    double stoppingTol = -1;
    string fileName;
    
    if (argc >= 4) {
        for (int i = 0; i < argc; i++) {
            int temp = 0;
            if (argv[i][0] == '-') {
                switch(argv[i][1]) {
                    case 'r':
                        seed = atoi(argv[i+1]);
                        i++;
                        break;
                    case 's':
                        temp = atoi(argv[i+1]);
                        if (temp == 1) {
                            skipInitFwd = true;
                        }
                        i++;
                        break;
                    case 't':
                        temp = atoi(argv[i+1]);
                        if (temp == 1) {
                            skipInitRev = true;
                        }
                        i++;
                        break;
                    case 'u':
                        temp = atoi(argv[i+1]);
                        if (temp == 1) {
                            skipRefineFwd = true;
                        }
                        i++;
                        break;
                    case 'v':
                        temp = atoi(argv[i+1]);
                        if (temp == 1) {
                            skipRefineRev = true;
                        }
                        i++;
                        break;
                    case 'e':
                        stoppingTol = atof(argv[i+1]);
                        i++;
                        break;
                }
            } else {
                fileName = argv[i];
            }
        }
    } else {
        fprintf(stderr, "Usage: ./RxnSysSim <fileName> -r <seed> [-s <skipInitFwd>] [-t <skipInitRev>] [-u <skipRefineFwd>]\n");
        fprintf(stderr, "                   [-v <skipRefineRev>] [-e <stoppingTol>] \n");
        fprintf(stderr, "            <fileName> Name of the netCDF file to read from and write to.\n");
        fprintf(stderr, "    -r          <seed> Integer used to seed the random number generator.\n");
        fprintf(stderr, "Optional arguments:\n");
        fprintf(stderr, "    -s   <skipInitFwd> Use this flag to indicate that preexisting data from the 'initFwdData' variable\n");
        fprintf(stderr, "                       should be used in the forward refinement algorithm.\n");
        fprintf(stderr, "    -t   <skipInitRev> Use this flag to indicate that preexisting data from the 'initRevData' variable\n");
        fprintf(stderr, "                       should be used in the reverse refinement algorithm.\n");
        fprintf(stderr, "    -u <skipRefineFwd> Use this flag to bypass the forward refinement algorithm.\n");
        fprintf(stderr, "    -v <skipRefineRev> Use this flag to bypass the reverse refinement algorithm.\n");
        fprintf(stderr, "    -e   <stoppingTol> Maximum allowable L2 distance between successive probability distributions.\n\n");
        abort();
    }
    
    FileInterface* fi = new FileInterface(fileName);
    std::mt19937 rng(seed);   
    
    int numTrials;
    double startTime;
    double endTime;
    int numTimePts;
    double timeStep;
    int numDataSavePts;
    int* dataSavePts;
    
    if (stoppingTol > 0) {
        fi->overwriteStoppingTol(stoppingTol);
    }
    
    System* masterSys = fi->readFileData(numTrials, startTime, endTime, numTimePts, timeStep, stoppingTol, numDataSavePts, dataSavePts);
            
    double* time = new double[numTimePts];  
    double** state = new double*[masterSys->numSpecies];
    for (int i = 0; i < masterSys->numSpecies; i++) {
        state[i] = new double[numTimePts];
    }
    
    double** avFwdSpeciesState = new double*[numTimePts];    
    double** avRevSpeciesState = new double*[numTimePts];
    double** lastAvFwdSpeciesState = new double*[numTimePts];
    double** lastAvRevSpeciesState = new double*[numTimePts];
    
    for (int i = 0; i < numTimePts; i++) {
        avFwdSpeciesState[i] = new double[masterSys->numSpecies];
        avRevSpeciesState[i] = new double[masterSys->numSpecies];
        lastAvFwdSpeciesState[i] = new double[masterSys->numSpecies];
        lastAvRevSpeciesState[i] = new double[masterSys->numSpecies];
                
        for (int j = 0; j < masterSys->numSpecies; j++) {
            avFwdSpeciesState[i][j] = 0;
            avRevSpeciesState[i][j] = 0;
            lastAvFwdSpeciesState[i][j] = 0;
            lastAvRevSpeciesState[i][j] = 0;
        }
    }       
    
    for (int i = 0; i < numTimePts; i++) {
        time[i] = startTime + timeStep * i;
    }    
    fi->writeTimeData(time, numTimePts);
    
    masterSys->timePts = time;
    masterSys->avFwdSpeciesState = avFwdSpeciesState;
    masterSys->avRevSpeciesState = avRevSpeciesState;
            
    string varName;
    double** lastFwdStatePt;
    double** lastRevStatePt;
    double*** initFwdData;
    double*** initRevData;
        
    omp_lock_t lock;
    omp_init_lock(&lock);
        
    System* sys;
    double initStart[numTrials];
    double initEnd[numTrials];
    double writeStart[numTrials];
    double writeEnd[numTrials];
    double revWriteStart[numTrials];
    double revWriteEnd[numTrials];
    double fwdWriteStart[numTrials];
    double fwdWriteEnd[numTrials];
    
    double progStart;
       
    if (!skipInitFwd) {                
        progStart = omp_get_wtime();
        
        #pragma omp parallel for default(shared) firstprivate(state) private(sys)
        for (int i = 0; i < numTrials; i++) {        
            initStart[i] = omp_get_wtime();        
            sys = new System(*masterSys);
            sys->init(rng());        
            initEnd[i] = omp_get_wtime();
            
            varName = "initFwdData";
            simFwd(sys, numTimePts, time, state);
            writeStateData(sys, fi, varName, i, state, numTimePts, revWriteStart, revWriteEnd, lock);
        }
    }
    
    if (!skipInitRev) {
        masterSys->rxnPq->minHeap = false;
        masterSys->time = time[numTimePts - 1];
        lastFwdStatePt = fi->readLastDataPt("initFwdData", numTrials, masterSys->numSpecies, numTimePts);     
        
        #pragma omp parallel for default(shared) firstprivate(state) private(sys)
        for (int i = 0; i < numTrials; i++) {
            sys = new System(*masterSys);
            sys->init(rng());
            
            for (int j = 0; j < sys->numSpecies; j++) {
                sys->speciesState[j] = lastFwdStatePt[i][j];
            }
            delete[] lastFwdStatePt[i];

            for (int j = 0; j < sys->numRxns; j++) {
                sys->rxns[j]->updateProp(sys->speciesState, sys->volRatio);
            }
            sys->initRev();

            varName = "initRevData";
            simRev(sys, numTimePts, time, state);
            writeStateData(sys, fi, varName, i, state, numTimePts, revWriteStart, revWriteEnd, lock);
        }
    }
    
    abort();
    
    if (false){//!skipRefineFwd) {
        initFwdData = fi->readStateData("initFwdState", numTrials, masterSys->numSpecies, numTimePts);         
        averageStateData(numTrials, numTimePts, masterSys->numSpecies, initFwdData, avFwdSpeciesState, lastAvFwdSpeciesState);
        
        int i = 0;
        while (true) {
            #pragma omp parallel for default(shared) firstprivate(state) private(sys)
            for (int j = 0; j < numTrials; j++) {        
                initStart[j] = omp_get_wtime();        
                sys = new System(*masterSys);
                sys->init(rng());        
                initEnd[j] = omp_get_wtime();

                varName = "fwdData";
                simFwd(sys, numTimePts, time, state);
                writeStateData(sys, fi, varName, j, state, numTimePts, revWriteStart, revWriteEnd, lock);
            }
        }
    }
    
    if (!skipRefineRev) {
        if (skipInitRev) {
            lastRevStatePt = fi->readLastDataPt("initRevState", numTrials, masterSys->numSpecies, numTimePts);            
            initRevData = fi->readStateData("initRevState", numTrials, masterSys->numSpecies, numTimePts);
        }
        
        // ***Simulate backward in time with current.
    }    
    
    abort();
    
    #pragma omp parallel for default(shared) firstprivate(state) private(sys)
    for (int i = 0; i < numTrials; i++) {        
        initStart[i] = omp_get_wtime();        
        sys = new System(*masterSys);
        sys->init(rng());        
        initEnd[i] = omp_get_wtime();
        
        if (skipInitRev) {
            for (int j = 0; j < sys->numSpecies; j++) {
                sys->speciesState[j] = lastRevStatePt[i][j];
            }
            delete[] lastRevStatePt[i];
            
            for (int j = 0; j < sys->numRxns; j++) {
                sys->rxns[j]->updateProp(sys->speciesState, sys->volRatio);
            }            
            sys->initFwd();
            
            revWriteStart[i] = writeEnd[i];
            revWriteEnd[i] = writeEnd[i];
        } else {
            if (skipInitFwd) {
                sys->time = time[numTimePts - 1];

                for (int j = 0; j < sys->numSpecies; j++) {
                    sys->speciesState[j] = lastFwdStatePt[i][j];
                }            
                delete[] lastFwdStatePt[i];

                for (int j = 0; j < sys->numRxns; j++) {
                    sys->rxns[j]->updateProp(sys->speciesState, sys->volRatio);
                }
                sys->initRev();

                writeStart[i] = initEnd[i];
                writeEnd[i] = initEnd[i];
            } else {
                simFwd(sys, numTimePts, time, state);
                writeStateData(sys, fi, varName, i, state, numTimePts, writeStart, writeEnd, lock);
            }
            
            varName = "revState";
            simRev(sys, numTimePts, time, state);
            writeStateData(sys, fi, varName, i, state, numTimePts, revWriteStart, revWriteEnd, lock);
        }
        
        if (skipRefineFwd) {
            fwdWriteStart[i] = revWriteEnd[i];
            fwdWriteEnd[i] = revWriteEnd[i];
        } else {            
            varName = "fwdState";
            simFwd(sys, numTimePts, time, state);
            writeStateData(sys, fi, varName, i, state, numTimePts, fwdWriteStart, fwdWriteEnd, lock);
        }
    }
    
    double progEnd = omp_get_wtime();
    
    double avInitTime = 0;
    double avExecTime = 0;
    double avWriteTime = 0;
    double avEwTime = 0;
    double avRevExecTime = 0;
    double avRevWriteTime = 0;
    double avRevEwTime = 0;
    double avFwdExecTime = 0;
    double avFwdWriteTime = 0;
    double avFwdEwTime = 0;
    double avTrialTime = 0;
    
    for (int i = 0; i < numTrials; i++) {
        avInitTime += initEnd[i] - initStart[i];
        
        avExecTime += writeStart[i] - initEnd[i];
        avWriteTime += writeEnd[i] - writeStart[i];
        avEwTime += writeEnd[i] - initEnd[i];

        avRevExecTime += revWriteStart[i] - writeEnd[i];
        avRevWriteTime += revWriteEnd[i] - revWriteStart[i];
        avRevEwTime += revWriteEnd[i] - writeEnd[i];
        
        avFwdExecTime += fwdWriteStart[i] - revWriteEnd[i];
        avFwdWriteTime += fwdWriteEnd[i] - fwdWriteStart[i];
        avFwdEwTime += fwdWriteEnd[i] - revWriteEnd[i];
        
        avTrialTime += fwdWriteEnd[i] - initStart[i];
    }
    avInitTime /= numTrials;
    avExecTime /= numTrials;
    avWriteTime /= numTrials;
    avEwTime /= numTrials;
    avRevExecTime /= numTrials;
    avRevWriteTime /= numTrials;
    avRevEwTime /= numTrials;
    avFwdExecTime /= numTrials;
    avFwdWriteTime /= numTrials;
    avFwdEwTime /= numTrials;
    avTrialTime /= numTrials;
    
    fprintf(stdout, "\nAverage system initialization time:    %e s\n\n", avInitTime);
    
    fprintf(stdout, "First forward simulation:\n");
    fprintf(stdout, "    Average simulation execution time: %e s\n", avExecTime);
    fprintf(stdout, "    Average file write execution time: %e s\n", avWriteTime);
    fprintf(stdout, "    Average execution time:            %e s\n\n", avEwTime);
    
    fprintf(stdout, "First reverse simulation:\n");
    fprintf(stdout, "    Average simulation execution time: %e s\n", avRevExecTime);
    fprintf(stdout, "    Average file write execution time: %e s\n", avRevWriteTime);
    fprintf(stdout, "    Average execution time:            %e s\n\n", avRevEwTime);
    
    fprintf(stdout, "Second forward simulation:\n");
    fprintf(stdout, "    Average simulation execution time: %e s\n", avFwdExecTime);
    fprintf(stdout, "    Average file write execution time: %e s\n", avFwdWriteTime);
    fprintf(stdout, "    Average execution time:            %e s\n\n", avFwdEwTime);
    
    fprintf(stdout, "Total program execution time:          %e s\n\n", progEnd - progStart);
        
    if (skipInitRev) {
        delete[] lastRevStatePt;
    } else if (skipInitFwd) {
        delete[] lastFwdStatePt;
    }
    
    delete masterSys;    
    delete fi;  
    
    return 0;
}

void simFwd(System* sys, int numTimePts, double* time, double** state) {
    for (int i = 0; i < numTimePts; i++) {
        while(sys->rxnPq->getNextTime() <= time[i]) {
            sys->execRxn(1);
        }
        
        sys->updateTime(time[i]);

        for (int j = 0; j < sys->numSpecies; j++) {
            state[j][i] = sys->speciesState[j];
        }
    }
}

void simRev(System* sys, int numTimePts, double* time, double** state) {
    for (int i = numTimePts - 1; i >= 0; i--) {
        double nextTime = sys->rxnPq->getNextTime();
        
        while(nextTime >= time[i] && nextTime != -DBL_MAX) {
            sys->execRxn(-1);
            nextTime = sys->rxnPq->getNextTime();
        }
        
        sys->updateTime(time[i]);
        
        for (int j = 0; j < sys->numSpecies; j++) {
            state[j][i] = sys->speciesState[j];
        }
    }
}

void writeStateData(System* sys, FileInterface* fi, string varName, int trial, double** state, int numTimePts, double* writeStart, double* writeEnd, omp_lock_t& lock) {        
    writeStart[trial] = omp_get_wtime();
    omp_set_lock(&lock);
    fi->writeStateData(varName, trial, state, sys->numSpecies, numTimePts);
    omp_unset_lock(&lock);        
    writeEnd[trial] = omp_get_wtime();
}

void averageStateData(int numTrials, int numTimePts, int numSpecies, double*** stateData, double** avStateData, double** lastAvStateData) {
    for (int i = 0; i < numTimePts; i++) {
        for (int j = 0; j < numSpecies; j++) {
            lastAvStateData[i][j] = avStateData[i][j];
            avStateData[i][j] = 0;
        }
    }
    
    for (int i = 0; i < numTrials; i++) {
        for (int j = 0; j < numTimePts; j++) {
            for (int k = 0; k < numSpecies; k++) {
                avStateData[j][k] += stateData[i][j][k];
            }
        }
    }
    
    for (int i = 0; i < numTimePts; i++) {
        for (int j = 0; j < numSpecies; j) {
            avStateData[i][j] /= numTrials;
        }
    }
}

double calcDist(double** avStateData1, double** avStateData2) {
    
}
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <float.h>
#include <random>
#include <string>
#include "main.h"
#include "Distribution.h"
#include "FileInterface.h"
#include "PriorityQueue.h"
#include "System.h"
#include "Reaction.h"
#include "Species.h"

using namespace std;

int main(int argc, const char* argv[]) {
    bool skipInitFwd = false;
    bool skipInitRev = false;
    bool skipRefine = false;
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
                            skipRefine = true;
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
        fprintf(stderr, "Usage: ./RxnSysSim <fileName> -r <seed> [-s <skipInitFwd>] [-t <skipInitRev>] [-u <skipRefine>]\n");
        fprintf(stderr, "                   [-e <stoppingTol>] \n");
        fprintf(stderr, "            <fileName> Name of the netCDF file to read from and write to.\n");
        fprintf(stderr, "    -r          <seed> Integer used to seed the random number generator.\n");
        fprintf(stderr, "Optional arguments:\n");
        fprintf(stderr, "    -s   <skipInitFwd> Use this flag to indicate that preexisting data from the 'initFwdData' variable\n");
        fprintf(stderr, "                       should be used in the forward refinement algorithm.\n");
        fprintf(stderr, "    -t   <skipInitRev> Use this flag to indicate that preexisting data from the 'initRevData' variable\n");
        fprintf(stderr, "                       should be used in the reverse refinement algorithm.\n");
        fprintf(stderr, "    -u    <skipRefine> Use this flag to bypass the forward and reverse refinement algorithms.\n");
        fprintf(stderr, "    -e   <stoppingTol> Maximum allowable L2 distance between successive probability distributions.\n\n");
        abort();
    }
    
    FileInterface* fi = new FileInterface(fileName);
    std::mt19937 rng(seed);   
        
    if (stoppingTol > 0) {
        fi->overwriteStoppingTol(stoppingTol);
    }
    
    int numTrials;
    double startTime;
    double endTime;
    int numTimePts;
    double timeStep;
    int numDataSavePts;
    int* dataSavePts;
    int numBoundedSpeciesStates;
    
    System* masterSys = fi->readFileData(numTrials, startTime, endTime, numTimePts, timeStep, stoppingTol, numDataSavePts, dataSavePts, numBoundedSpeciesStates);
            
    double* time = new double[numTimePts];  
    
    int numThreads = omp_get_max_threads();
    
    double*** state = new double**[numThreads];
    for (int i = 0; i < numThreads; i++) {
        state[i] = new double*[masterSys->numSpecies];
        for (int j = 0; j < masterSys->numSpecies; j++) {
            state[i][j] = new double[numTimePts];
        }
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
    
    int j = 0;
    int* distSpeciesKey = new int[numBoundedSpeciesStates];
    int* speciesDistKey = new int[masterSys->numSpecies];
    for (int i = 0; i < masterSys->numSpecies; i++) {
        if (masterSys->species[i]->stateBounded) {
            distSpeciesKey[j] = i;
            speciesDistKey[i] = j;
            j++;
        } else {
            speciesDistKey[i] = -1;
        }
    }
    
    Distribution** revDistsPrev = new Distribution*[numBoundedSpeciesStates];
    for (int i = 0; i < numBoundedSpeciesStates; i++) {
        revDistsPrev[i] = new Distribution(masterSys->species[distSpeciesKey[i]], numTrials, rng());
    }
    
    Distribution** revDists = new Distribution*[numBoundedSpeciesStates];
    for (int i = 0; i < numBoundedSpeciesStates; i++) {
        revDists[i] = new Distribution(masterSys->species[distSpeciesKey[i]], numTrials, rng());
    }
    
    Distribution** fwdDistsPrev = new Distribution*[numBoundedSpeciesStates];
    for (int i = 0; i < numBoundedSpeciesStates; i++) {
        fwdDistsPrev[i] = new Distribution(masterSys->species[distSpeciesKey[i]], numTrials, rng());
    }
    
    Distribution** fwdDists = new Distribution*[numBoundedSpeciesStates];
    for (int i = 0; i < numBoundedSpeciesStates; i++) {
        fwdDists[i] = new Distribution(masterSys->species[distSpeciesKey[i]], numTrials, rng());
    }
    
    if (!skipInitFwd) {                
        progStart = omp_get_wtime();
        
        for (int i = 0; i < numBoundedSpeciesStates; i++) {
            revDists[i]->addNode(revDists[i]->species->state, numTrials);
        }
        
        #pragma omp parallel for default(shared) private(sys)
        for (int i = 0; i < numTrials; i++) {        
            int threadId = omp_get_thread_num();            
            
            initStart[i] = omp_get_wtime();        
            sys = new System(*masterSys);
            sys->seed(rng());
            sys->initFwd();
            initEnd[i] = omp_get_wtime();
            
            varName = "initFwdData";
            simFwd(sys, numTimePts, time, state[threadId], revDists, speciesDistKey);
            writeStateData(-1, sys, fi, varName, i, state[threadId], numTimePts, lock);
        }
    }
    
    if (!skipInitRev) {
        masterSys->rxnPq->minHeap = false;
        masterSys->time = time[numTimePts - 1];
        lastFwdStatePt = fi->readInitDataPt("initFwdData", numTrials, masterSys->numSpecies, numTimePts - 1, numTimePts);
        
        for (int i = 0; i < numBoundedSpeciesStates; i++) {             
            fwdDists[i]->update(lastFwdStatePt, numTrials);
        }
        
        #pragma omp parallel for default(shared) private(sys)
        for (int i = 0; i < numTrials; i++) {    
            int threadId = omp_get_thread_num(); 
            
            sys = new System(*masterSys);
            sys->seed(rng());            
            
            for (int j = 0; j < sys->numSpecies; j++) {
                sys->species[j]->state = lastFwdStatePt[i][j];
            }
            delete[] lastFwdStatePt[i];

            for (int j = 0; j < sys->numRxns; j++) {
                sys->rxns[j]->updateProp(sys->volRatio);
            }
            sys->initRev();

            varName = "initRevData";
            simRev(sys, numTimePts, time, state[threadId], fwdDists, speciesDistKey);
            writeStateData(-1, sys, fi, varName, i, state[threadId], numTimePts, lock);
        }
    }  
    
    if (!skipRefine) {        
        for (int i = 0; i < numDataSavePts; i++) {
            masterSys->rxnPq->minHeap = true;
            masterSys->time = time[0];
            if (i == 0) {
                lastRevStatePt = fi->readInitDataPt("initRevData", numTrials, masterSys->numSpecies, 0, numTimePts);
            } else {
                lastRevStatePt = fi->readDataPt("revData", i - 1, numTrials, masterSys->numSpecies, 0, numTimePts);
            }

            for (int i = 0; i < numBoundedSpeciesStates; i++) {             
                revDists[i]->update(lastRevStatePt, numTrials);
            }
            
            #pragma omp parallel for default(shared) private(sys)
            for (int j = 0; j < numTrials; j++) {                
                int threadId = omp_get_thread_num(); 

                sys = new System(*masterSys);
                sys->seed(rng());            

                for (int k = 0; k < sys->numSpecies; k++) {
                    sys->species[k]->state = lastRevStatePt[j][k];
                }
                delete[] lastRevStatePt[j];

                for (int k = 0; k < sys->numRxns; k++) {
                    sys->rxns[k]->updateProp(sys->volRatio);
                }
                sys->initFwd();

                varName = "fwdData";
                simFwd(sys, numTimePts, time, state[threadId], revDists, speciesDistKey);
                writeStateData(i, sys, fi, varName, j, state[threadId], numTimePts, lock);
            }
            
            masterSys->rxnPq->minHeap = false;
            masterSys->time = time[numTimePts - 1];
            lastFwdStatePt = fi->readDataPt("fwdData", i, numTrials, masterSys->numSpecies, numTimePts - 1, numTimePts);

            for (int i = 0; i < numBoundedSpeciesStates; i++) {             
                fwdDists[i]->update(lastFwdStatePt, numTrials);
            }
            
            #pragma omp parallel for default(shared) private(sys)
            for (int j = 0; j < numTrials; j++) {                
                int threadId = omp_get_thread_num(); 

                sys = new System(*masterSys);
                sys->seed(rng());            

                for (int k = 0; k < sys->numSpecies; k++) {
                    sys->species[k]->state = lastFwdStatePt[j][k];
                }
                delete[] lastFwdStatePt[j];

                for (int k = 0; k < sys->numRxns; k++) {
                    sys->rxns[k]->updateProp(sys->volRatio);
                }
                sys->initRev();

                varName = "revData";
                simRev(sys, numTimePts, time, state[threadId], revDists, speciesDistKey);
                writeStateData(i, sys, fi, varName, j, state[threadId], numTimePts, lock);
            }
        }
    }
    
    abort();
    
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

void simFwd(System* sys, int numTimePts, double* time, double** state, Distribution** dists, int* speciesDistKey) {
    int boundBreachSpeciesId;
    
    for (int i = 0; i < numTimePts; i++) {
        while(sys->rxnPq->getNextTime() <= time[i]) {
            boundBreachSpeciesId = sys->execRxn(true);
            
            if (boundBreachSpeciesId >= 0) {
                i = 0;
                sys->updateTime(time[i]);
                sys->species[boundBreachSpeciesId]->state = dists[speciesDistKey[boundBreachSpeciesId]]->sample();
                sys->initFwd();
            }
        }
        
        sys->updateTime(time[i]);

        for (int j = 0; j < sys->numSpecies; j++) {
            state[j][i] = sys->species[j]->state;
        }
    }
}

void simRev(System* sys, int numTimePts, double* time, double** state, Distribution** dists, int* speciesDistKey) {
    int boundBreachSpeciesId;
    
    for (int i = numTimePts - 1; i >= 0; i--) {
        double nextTime = sys->rxnPq->getNextTime();
        
        while(nextTime >= time[i] && nextTime != -DBL_MAX) {
            boundBreachSpeciesId = sys->execRxn(false);
            nextTime = sys->rxnPq->getNextTime();
            
            if (boundBreachSpeciesId >= 0) {
                i = numTimePts - 1;
                sys->updateTime(time[i]);
                sys->species[boundBreachSpeciesId]->state = dists[speciesDistKey[boundBreachSpeciesId]]->sample();
                sys->initRev();
            }
        }
        
        sys->updateTime(time[i]);
        
        for (int j = 0; j < sys->numSpecies; j++) {
            state[j][i] = sys->species[j]->state;
        }
    }
}

void writeStateData(int dataSavePtId, System* sys, FileInterface* fi, string varName, int trial, double** state, int numTimePts, omp_lock_t& lock) {
    omp_set_lock(&lock);
    if (dataSavePtId < 0) { 
        fi->writeInitStateData(varName, trial, state, sys->numSpecies, numTimePts);
    } else {        
        fi->writeStateData(varName, dataSavePtId, trial, state, sys->numSpecies, numTimePts);
    }
    omp_unset_lock(&lock);
}

/*
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
}*/
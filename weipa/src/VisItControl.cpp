
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <weipa/VisItControl.h>
#include <weipa/EscriptDataset.h>

#ifdef USE_VISIT
#include <weipa/VisItData.h>

#include <VisItControlInterface_V2.h>

#include <cstring>
#include <iostream>

#define VISIT_COMMAND_PROCESS 0
#define VISIT_COMMAND_SUCCESS 1
#define VISIT_COMMAND_FAILURE 2
#endif

namespace weipa {

namespace VisItControl {
    
bool initialized = false;

#ifdef USE_VISIT

weipa::VisItData_ptr visitData(new VisItData());
int mpiRank = 0;
int mpiSize = 1;
bool runFlag = true;
bool connected = false;

// Helper function for processVisItCommand()
static void broadcastSlaveCommand(int* command)
{
#if WEIPA_HAVE_MPI
    MPI_Bcast(command, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
}

// Processes commands from the VisIt viewer on all processors
static int processVisItCommand()
{
    if (mpiRank == 0) {
        int success = VisItProcessEngineCommand();
        int command;

        if (success) {
            command = VISIT_COMMAND_SUCCESS;
            broadcastSlaveCommand(&command);
            return 1;
        } else {
            command = VISIT_COMMAND_FAILURE;
            broadcastSlaveCommand(&command);
            return 0;
        }
    } else {
        // Note: Only through the SlaveProcessCallback() function
        // can rank 0 send VISIT_COMMAND_PROCESS to non-zero ranks
        while (1) {
            int command=VISIT_COMMAND_SUCCESS;
            broadcastSlaveCommand(&command);
            switch (command) {
                case VISIT_COMMAND_PROCESS:
                    VisItProcessEngineCommand();
                    break;
                case VISIT_COMMAND_SUCCESS:
                    return 1;
                case VISIT_COMMAND_FAILURE:
                    return 0;
            }
        }
    }
}

// Callback involved in command communication
static void slaveProcessCallback()
{
    int command = VISIT_COMMAND_PROCESS;
    broadcastSlaveCommand(&command);
}

#if WEIPA_HAVE_MPI
static int broadcastIntCallback(int* value, int sender)
{
    return MPI_Bcast(value, 1, MPI_INT, sender, MPI_COMM_WORLD);
}

static int broadcastStringCallback(char* str, int len, int sender)
{
    return MPI_Bcast(str, len, MPI_CHAR, sender, MPI_COMM_WORLD);
}
#endif

// Processes a control command
static void controlCommandCallback(const char* cmd, const char* sdata,
                                   void* cb_data)
{
#ifdef _DEBUG
    std::cout << "Control Command: " << cmd << std::endl;
#endif

    if (strstr(cmd, "pause")) {
        runFlag = !runFlag;
    } else if (strstr(cmd, "update_plots")) {
        VisItUpdatePlots();
    }
}

// Callback that returns metadata for this dataset
visit_handle VisItGetMetaData(void* cbdata)
{
#ifdef _DEBUG
    std::cout << "VisItGetMetaData()" << std::endl;
#endif
    return visitData->getSimMetaData();
}

// Callback that returns the domain list for this dataset
visit_handle VisItGetDomainList(const char* name, void* cbdata)
{
#ifdef _DEBUG
    std::cout << "VisItGetDomainList(" << name << ")" << std::endl;
#endif
    return visitData->getDomainList();
}

// Callback that returns mesh data
visit_handle VisItGetMesh(int domain, const char* name, void* cbdata)
{
#ifdef _DEBUG
    std::cout << "VisItGetMesh(" << domain << ", '" << name << "')" << std::endl;
#endif
    if (mpiRank != domain) {
#ifdef _DEBUG
        std::cout << "I don't have data for domain " << domain << std::endl;
#endif
        return VISIT_INVALID_HANDLE;
    }
    return visitData->getMesh(name);
}

// Callback that returns variable data
visit_handle VisItGetVariable(int domain, const char* name, void* cbdata)
{
#ifdef _DEBUG
    std::cout << "VisItGetVariable(" << domain << ", '" << name << "')" << std::endl;
#endif
    if (mpiRank != domain) {
#ifdef _DEBUG
        std::cout << "I don't have data for domain " << domain << std::endl;
#endif
        return VISIT_INVALID_HANDLE;
    }
    return visitData->getVariable(name);
}

#endif // USE_VISIT

/**************************************************************************************************/

// Initializes libsim and the communication to VisIt
bool initialize(const std::string& simFile, const std::string& comment)
{
#ifdef USE_VISIT
    if (connected) {
        VisItDisconnect();
        connected = false;
        std::cout << "VisIt client disconnected." << std::endl;
    }

    if (!initialized) {
        //VisItOpenTraceFile("/tmp/simV2trace.txt");
        if (!VisItSetupEnvironment()) {
            std::cerr << "Error setting up VisIt environment" << std::endl;
            return false;
        }

#if WEIPA_HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

        // Install callback functions for global communication
        VisItSetBroadcastIntFunction(broadcastIntCallback);
        VisItSetBroadcastStringFunction(broadcastStringCallback);

        // Tell libsim whether the simulation is parallel
        VisItSetParallel(mpiSize > 1);
        VisItSetParallelRank(mpiRank);
#endif
        std::vector<std::string> cmdNames;
        cmdNames.push_back("pause");
        cmdNames.push_back("update_plots");
        visitData->setCommandNames(cmdNames);
        initialized = true;
    }

    if (mpiRank == 0) {
        std::string filename(simFile);
        if (filename.rfind(".sim2") != filename.length()-5) {
            filename += ".sim2";
        }
        VisItInitializeSocketAndDumpSimFile("escript", comment.c_str(),
                NULL, NULL, NULL, filename.c_str());
    }
#endif // USE_VISIT
    return initialized;
}

// Main entry point that checks for client input and publishes new data
bool publishData(EscriptDataset_ptr dataset)
{
    if (!initialized || dataset->getConvertedDomain().size()==0) {
        return false;
    }

#ifdef USE_VISIT
    int visitState = 0, err = 0;

    if (connected) {
        visitData->publishData(dataset);
        visitData->setSimulationStatus(runFlag);
        VisItTimeStepChanged();
    }

    do {
        // Get input from VisIt or timeout so the simulation can continue
        if (mpiRank == 0) {
            int blocking = (connected && !runFlag) ? 1 : 0;
            visitState = VisItDetectInput(blocking, -1);
        }

#if WEIPA_HAVE_MPI
        MPI_Bcast(&visitState, 1, MPI_INT, 0, dataset->getMPIComm());
#endif

        // visitState values:
        //   0: there was no input from VisIt so nothing to do
        //   1: VisIt is trying to connect
        //   2: VisIt sent a command which should be processed
        // < 0: an unrecoverable error occurred
        if (visitState < 0) {
            std::cerr << "VisItControl: Can't recover from error! State="
                << visitState << std::endl; fflush(stderr);
            err = 1;
        } else if (visitState == 1) {
            if (VisItAttemptToCompleteConnection() == VISIT_OKAY) {
                std::cout << "VisIt client connected." << std::endl;
                // publish latest data
                visitData->publishData(dataset);
                visitData->setSimulationStatus(runFlag);
                void* cbdata = NULL;
                VisItSetCommandCallback(controlCommandCallback, cbdata);
                VisItSetSlaveProcessCallback(slaveProcessCallback);

                VisItSetGetMetaData(VisItGetMetaData, cbdata);
                VisItSetGetMesh(VisItGetMesh, cbdata);
                VisItSetGetVariable(VisItGetVariable, cbdata);
                VisItSetGetDomainList(VisItGetDomainList, cbdata);
                VisItTimeStepChanged();
                connected = true;
            } else {
                char *errorString = VisItGetLastError();
                if (mpiRank == 0) {
                    if (strlen(errorString) > 0) {
                        std::cerr << errorString << std::endl;
                    } else {
                        std::cerr << "VisIt failed to connect successfully."
                            << std::endl;
                    }
                }
                free(errorString);
            }
        } else if (visitState == 2) {
            if (!processVisItCommand()) {
                VisItDisconnect();
                connected = false;
                std::cout << "VisIt client disconnected." << std::endl;
            }
        }
    } while (visitState != 0);

    return err==0;
#else // USE_VISIT

    return false;
#endif
}

} // namespace VisItControl

} // namespace weipa


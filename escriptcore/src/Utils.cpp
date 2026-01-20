
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "Data.h"
#include "DataVector.h"
#include "EsysMPI.h"
#include "FileWriter.h"
#include "Utils.h"

#include <cstring>
#include <fstream>
#ifndef _WIN32
# include <unistd.h>
#endif // _WIN32
#include <string>
#include <cstdlib>

#include <boost/python.hpp>
#include <boost/scoped_array.hpp>

#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#include "DataTypes.h"
#endif

#ifdef ESYS_HAVE_NETCDF4
#include <ncFile.h>
#endif

// #define RELEASE "5.7"

namespace bp = boost::python;

#ifndef OVERLORDPATH
#define OVERLORDPATH ""
#endif

namespace escript {

int getSvnVersion()
{
#ifdef SVN_VERSION
    return SVN_VERSION;
#else
    return 0;
#endif
}

std::string getGitVersion()
{
#ifdef RELEASE
    return RELEASE;
#else
    char chardate[100]={-1};
    FILE *p = popen("git show -s --format=%ct","r");
    if(p != NULL) 
    {
        while(fgets(chardate, sizeof(chardate), p) != NULL)
        {}
    }
    pclose(p);
    std::string answer = "";
    int counter=0;
    while(1)
    {
        if(chardate[counter] == '\000')
            break;
        answer+=chardate[counter];
        counter++;
    }
    return answer;
#endif
}


// This is probably not very robust, but it works on Savanna today and is
// useful for performance analysis
int get_core_id()
{
    int processor_num=-1;
#ifdef CORE_ID1
    FILE *fp;
    int i, count_spaces=0;
    char fname[100];
    char buf[1000];

    sprintf(fname, "/proc/%d/stat", getpid());
    fp = fopen(fname, "r");
    if (fp == NULL) return(-1);
    fgets(buf, 1000, fp);
    fclose(fp);

    for (i=strlen(buf)-1; i>=0; i--) {
        if (buf[i] == ' ') count_spaces++;
        if (count_spaces == 4) break;
    }
    processor_num = atoi(&buf[i+1]);
#endif
    return processor_num;
}

void setNumberOfThreads(int num_threads)
{
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#endif
}

int getNumberOfThreads()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

int getMPISizeWorld()
{
    int mpi_num = 1;
#ifdef ESYS_MPI
    ensureMPIInitialized();
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_num);
#endif
    return mpi_num;
}

int getMPIRankWorld()
{
    int mpi_iam = 0;
#ifdef ESYS_MPI
    ensureMPIInitialized();
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_iam);
#endif
    return mpi_iam;
}

int getMPIWorldMax(int val)
{
#ifdef ESYS_MPI
    ensureMPIInitialized();
    int val2 = val;
    int out = val;
    MPI_Allreduce(&val2, &out, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
    int out = val;
#endif
    return out;
}

int getMPIWorldSum(int val)
{
#ifdef ESYS_MPI
    ensureMPIInitialized();
    int val2 = val;
    int out = 0;
    MPI_Allreduce(&val2, &out, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
    int out = val;
#endif
    return out;
}

void printParallelThreadCnt()
{
    char hname[64];

#ifdef HAVE_GETHOSTNAME
    gethostname(hname, 64);
    hname[63] = '\0';
#else
    strcpy(hname, "unknown host");
#endif

    const int mpi_num = getMPISizeWorld();
    const int mpi_iam = getMPIRankWorld();

#pragma omp parallel
    {
        int omp_iam=0, omp_num=1;
#ifdef _OPENMP
        omp_iam = omp_get_thread_num();
        omp_num = omp_get_num_threads();
#endif
#pragma omp critical (printthrdcount)
        printf("printParallelThreadCounts: MPI=%03d/%03d OpenMP=%03d/%03d "
               "running on %s core %d\n", mpi_iam, mpi_num, omp_iam, omp_num,
               hname, get_core_id());
    }
}

#define CHILD_FAIL 2
#define CHILD_COMPLETE 4

#ifndef _WIN32
#ifdef ESYS_MPI
#include <sys/socket.h>
#include <sys/select.h>
#include <errno.h>
#include <arpa/inet.h>
int prepareSocket(unsigned short *port, int *key)
{
    if (getMPIRankWorld() != 0)
        return 0;
    int sfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sfd < 0) {
        perror("socket creation failure");
        return -1;
    }
    int opt = 1;
    if (setsockopt(sfd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(int)) < 0) {
        perror("socket option setting failure");
        close(sfd);
        return -1;
    }

    struct sockaddr_in addr;
    addr.sin_family = AF_INET;
    addr.sin_port = htons(0);
    addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    memset(&addr.sin_zero, 0, sizeof(addr.sin_zero));

    if (bind(sfd, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
        perror("bind failure");
        close(sfd);
        return -1;
    }

    if (listen(sfd, SOMAXCONN) < 0) {
        perror("listen failure");
        close(sfd);
        return -1;
    }

    struct sockaddr actual;
    unsigned int size = sizeof(actual);
    if (getsockname(sfd, &actual, &size) < 0) {
        perror("failed when determining bound port number");
        close(sfd);
        return -1;
    }

    //if size > sizeof(actual), some info was truncated, but certainly not port

    *port = ntohs(((struct sockaddr_in *) &actual)->sin_port);

    unsigned int seed = time(NULL) % UINT_MAX;
    *key = rand_r(&seed);
    return sfd;
}

int check_data(unsigned int max, fd_set *all, fd_set *valid, int key, int sfd)
{
    int ret = 0;
    for (int i = 0; i <= max; i++) {
        if (i == sfd)
            continue;
        if (FD_ISSET(i, all)) {
            int provided = 0;
            if (recv(i, &provided, sizeof(int), MSG_WAITALL) == sizeof(int)
                    && provided == key) {
                char deadspace[1024];
                while ((provided = recv(i, deadspace, 1024, 0))) {
                    if (provided == -1) {
                        if (errno == EAGAIN || errno == EWOULDBLOCK) {
                            continue;
                        } else {
                            perror("connection failure");
                            return CHILD_FAIL;
                        }
                    }
                }
                return CHILD_COMPLETE;
            } else {
                FD_CLR(i, all);
                close(i);
            }
        }
    }
    return ret;
}

void close_all(unsigned int maxfd, fd_set *all)
{
    for (int i = 0; i <= maxfd; i++) {
        if (FD_ISSET(i, all))
            close(i);
    }
}

int waitForCompletion(int sfd, int key)
{
    if (getMPIRankWorld() != 0)
        return 0;
    int timeout = 10; //max of 10 seconds with no communication

    fd_set all, valid;
    FD_ZERO(&all);
    FD_ZERO(&valid);
    FD_SET(sfd, &all);
    time_t last_good_time = time(NULL);
    unsigned int maxfd = sfd;
    while (time(NULL) - last_good_time < timeout) {
        struct timeval timer = {1,0}; //1 sec, 0 usec
        int count = select(maxfd + 1, &all, NULL, NULL, &timer);
        if (count == -1) { //error
            if (errno == EINTR) {
                continue; //just a signal, continue as we were
            } else {
                perror("socket operation error");
                close(sfd);
                return -1;
            }
        } else if (FD_ISSET(sfd, &all)) { //new connection
            int connection = accept(sfd, NULL, NULL);
            if (connection > maxfd)
                maxfd = connection;
            FD_SET(connection, &all);
            FD_CLR(connection, &valid);
            time(&last_good_time);
            count--;
        }
        if (count > 0) { //something to read, either connection key or state
            int res = check_data(maxfd, &all, &valid, key, sfd);
            if (res == CHILD_FAIL) {
                return -1;
            } else if (res == CHILD_COMPLETE) {
                close_all(maxfd, &all);
                return 0;
            }
        }
    }
    close_all(maxfd, &all);
    fprintf(stderr, "Connection to child process timed out\n");
    return -1;
}
#endif //ESYS_MPI
#endif //not _WIN32

int runMPIProgram(bp::list args)
{
#ifdef ESYS_MPI
    ensureMPIInitialized();
    unsigned short port = 0;
    int key = 0;
    int sock = prepareSocket(&port, &key);
    if (getMPIWorldSum(sock) < 0)
        return -1;
    char portstr[20] = {'\0'}, keystr[20] = {'\0'};
    sprintf(portstr, "%d", port);
    sprintf(keystr, "%d", key);

    int nargs = bp::extract<int>(args.attr("__len__")());
    // make room for port, key and terminator
    char** c_args = new char*[nargs + 3];
    c_args[0] = portstr;
    c_args[1] = keystr;
    std::vector<std::string> cpp_args(nargs);
    // skip command name in argument list
    for (int i = 0; i < nargs; i++) {
        cpp_args[i] = bp::extract<std::string>(args[i]);
        c_args[i+2] = const_cast<char*>(cpp_args[i].c_str());
    }
    c_args[nargs+2] = NULL;

    MPI_Info info;
    MPI_Info_create(&info);
    // force the gmsh process to run on this host as well for network comm
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int temp = MPI_MAX_PROCESSOR_NAME;
    MPI_Get_processor_name(hostname, &temp);
    char hoststr[] = "host"; //for warnings
    MPI_Info_set(info, hoststr, hostname);
    MPI_Comm intercomm;
    int errors;
    char c_cmd[] = OVERLORDPATH"escript-overlord";
    MPI_Comm_spawn(c_cmd, c_args, 1, info, 0, MPI_COMM_WORLD, &intercomm, &errors);
    MPI_Info_free(&info);
    delete[] c_args;
    if (errors != MPI_SUCCESS)
        return errors;
    return getMPIWorldMax(waitForCompletion(sock, key));
#else //#ifdef ESYS_MPI
    std::string cmd;
    int nargs = bp::extract<int>(args.attr("__len__")());
    for (int i = 0; i < nargs; i++) {
        cmd += bp::extract<std::string>(args[i]);
        cmd += " ";
    }
    return system(cmd.c_str());
#endif //#ifdef ESYS_MPI/else
}
#undef CHILD_COMPLETE
#undef CHILD_FAIL

double getMachinePrecision()
{
    return std::numeric_limits<double>::epsilon();
}

double getMaxFloat()
{
    return std::numeric_limits<double>::max();
}

void MPIBarrierWorld()
{
#ifdef ESYS_MPI
    ensureMPIInitialized();
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void saveDataCSV(const std::string& filename, bp::dict arg,
                 const std::string& sep, const std::string& csep, bool refid, bool append)
{
    bp::list keys = arg.keys();
    int numdata = bp::extract<int>(arg.attr("__len__")());

    bool hasmask = arg.has_key("mask");
    Data mask;
    if (hasmask) {
        mask = bp::extract<escript::Data>(arg["mask"]);
        if (mask.getDataPointRank() != 0) {
            throw DataException("saveDataCSV: mask must be scalar.");
        }
        keys.remove("mask");
        numdata--;
    }
    if (numdata < 1) {
        throw DataException("saveDataCSV: no data to save specified.");
    }

    std::vector<int> step(numdata);
    std::vector<std::string> names(numdata);
    std::vector<Data> data(numdata);
    std::vector<int> fstypes(numdata+int(hasmask)); // FunctionSpace types for each data for interpolation

    keys.sort(); // to get some predictable order to things

    // We need to interpret the samples correctly even if they are different
    // types. For this reason, we should iterate over samples...
    for (int i=0; i<numdata; ++i) {
        names[i] = bp::extract<std::string>(keys[i]);
        data[i] = bp::extract<escript::Data>(arg[keys[i]]);
        fstypes[i] = data[i].getFunctionSpace().getTypeCode();
        step[i] = (data[i].actsExpanded() ? DataTypes::noValues(data[i].getDataPointShape()) : 0);

        if (i > 0) {
            if (data[i].getDomain()!=data[i-1].getDomain()) {
                throw DataException("saveDataCSV: all data must be on the same domain.");
            }
        }
        if (data[i].isComplex()) {
            throw DataException("saveDataCSV: complex values must be separated into components before calling this.");
        }
    }

    if (hasmask) {
        if (mask.getDomain() != data[0].getDomain())
            throw DataException("saveDataCSV: mask domain must be the same as the data domain.");
        fstypes[numdata] = mask.getFunctionSpace().getTypeCode();
    }

    int bestfnspace = 0;
    if (!data[0].getDomain()->commonFunctionSpace(fstypes, bestfnspace)) {
        throw DataException("saveDataCSV: FunctionSpaces of data are incompatible");
    }


    // now we interpolate all data to the same type
    FunctionSpace best(data[0].getDomain(), bestfnspace);
    for (int i=0; i<numdata; ++i) {
        data[i] = data[i].interpolate(best);
    }
    if (hasmask)
        mask = mask.interpolate(best);

    // these must be the same for all data
    int numsamples = data[0].getNumSamples();
    int dpps = data[0].getNumDataPointsPerSample();
    std::ostringstream os;
    bool first = true;

    if (data[0].getDomain()->getMPIRank() == 0) {
        for (int i=0; i<numdata; ++i) {
            const DataTypes::ShapeType& s=data[i].getDataPointShape();
            switch (data[i].getDataPointRank()) {
                case 0:
                    if (!first) {
                        os << sep;
                    } else {
                        first=false;
                        if(refid){
                            os << "Ref_ID" << sep;
                        }
                    }
                    os << names[i];
                    break;

                case 1:
                    for (int j=0; j<s[0]; ++j) {
                        if (!first) {
                            os << sep;
                        } else {
                            first=false;
                            if(refid){
                                os << "Ref_ID" << sep;
                            }
                        }
                        os << names[i] << csep << j;
                    }
                    break;

                case 2:
                    for (int j=0; j<s[0]; ++j) {
                        for (int k=0; k<s[1]; ++k) {
                            if (!first) {
                                os << sep;
                            } else {
                                first=false;
                                if(refid){
                                    os << "Ref_ID" << sep;
                                }
                            }
                            os << names[i] << csep << k << csep << j;
                        }
                    }
                    break;

                case 3:
                    for (int j=0; j<s[0]; ++j) {
                        for (int k=0; k<s[1]; ++k) {
                            for (int l=0; l<s[2]; ++l) {
                                if (!first) {
                                    os << sep;
                                } else {
                                    first=false;
                                    if(refid){
                                        os << "Ref_ID" << sep;
                                    }
                                }
                                os << names[i] << csep << k << csep << j
                                   << csep << l;
                            }
                        }
                    }
                    break;

                case 4:
                    for (int j=0; j<s[0]; ++j) {
                        for (int k=0; k<s[1]; ++k) {
                            for (int l=0; l<s[2]; ++l) {
                                for (int m=0; m<s[3]; ++m) {
                                    if (!first) {
                                        os << sep;
                                    } else {
                                        first=false;
                                        if(refid){
                                            os << "Ref_ID" << sep;
                                        }
                                    }
                                    os << names[i] << csep << k << csep << j
                                       << csep << l << csep << m;
                                }
                            }
                        }
                    }
                    break;

                default:
                    throw DataException("saveDataCSV: Illegal rank");
            }
        }
        os << std::endl;
    }

    const double* masksample = NULL;

    // does the mask act expanded?
    // Are there mask value for each point in the sample?
    bool expandedmask = false;
    // do we output this row?
    bool wantrow = true;
    if (hasmask && mask.actsExpanded()) {
        expandedmask=true;
    }
    os.setf(std::ios_base::scientific, std::ios_base::floatfield);
    os.precision(15);

    // errors prior to this point will occur on all processes anyway
    // so there is no need to explicitly notify other ranks
    int error = 0;
    std::string localmsg;
    try {
        std::vector<int> offset(numdata);
        std::vector<const DataTypes::real_t*> samples(numdata);
#ifdef ESYS_MPI
        // Loop over the MPI ranks
        for(int token = 0; token < data[0].getDomain()->getMPISize(); token++)
        {
            data[0].getDomain()->MPIBarrier();
            if(data[0].getDomain()->getMPIRank() != token)
                continue;
#endif

	    	const DataTypes::real_t onlyreal=0;
            for (int i=0; i<numsamples; ++i) {
                if (!best.ownSample(i)) {
                    continue;
                }
                wantrow = true;
                for (int d=0; d<numdata; ++d) {
                    samples[d] = data[d].getSampleDataRO(i, onlyreal);
                }
                if (hasmask) {
                    masksample = mask.getSampleDataRO(i, onlyreal);
                    if (!expandedmask) {
                        // mask controls whole sample
                        if (masksample[0] <= 0) {
                            wantrow = false;
                        }
                    }
                }
                for (int j=0; j<dpps; ++j) {
                    // now we need to check if this point is masked off
                    if (expandedmask) {
                        // masks are scalar so the relevant value is at [j]
                        wantrow = (masksample[j]>0);
                    }
                    if (wantrow) {
                        bool needsep = false;
                        // If necessary, write the element id to the output stream
                        if(refid){
                            os << best.getReferenceIDOfSample(i) << sep;
                        }
                        for (int d = 0; d < numdata; ++d) {
                            DataTypes::pointToStream(os, samples[d],
                                data[d].getDataPointShape(), offset[d], needsep, sep);
                            needsep = true;
                            offset[d] += step[d];
                        }
                        os << std::endl;
                    }
                }
                for (int d = 0; d<numdata; d++) {
                    offset[d]=0;
                }
            }
#ifdef ESYS_MPI
        }
#endif
    } catch (EsysException* e) {
        error=1;
        if (data[0].getDomain()->getMPISize()==1) {
            throw;
        } else {
            localmsg=e->what();
        }
    } catch (...) {
        error=1;
        if (data[0].getDomain()->getMPISize()==1) {
            throw;
        }
    }

#ifdef ESYS_MPI
    MPI_Comm com = data[0].get_MPIComm();
    int rerror = 0;
    MPI_Allreduce(&error, &rerror, 1, MPI_INT, MPI_MAX, com);
    error = rerror;
#else
    MPI_Comm com = MPI_COMM_NULL;
#endif

    if (error) {
        if (localmsg.empty()) {
            throw DataException("saveDataCSV: error building output");
        }
        throw DataException(std::string("saveDataCSV:")+localmsg);
    }

    // at this point os will contain the text to be written
    FileWriter fw(com);
    if (!fw.openFile(filename, 0, false, append)) {
        throw DataException("saveDataCSV: unable to open file for writing");
    }
    error = !fw.writeOrdered(os);
    fw.close();

    data[0].getDomain()->MPIBarrier();

    if (error)
        throw DataException("saveDataCSV: Error writing to file");
}

#ifdef ESYS_HAVE_BOOST_NUMPY
boost::python::list getNumpy(boost::python::dict arg)
{
    // Initialise boost python numpy
    bp::numpy::initialize();

    // Extract the key information and the number of data ppoints
    bp::list keys = arg.keys();
    int numdata = bp::extract<int>(arg.attr("__len__")());

    // Process the mask information
    bool hasmask = arg.has_key("mask");
    Data mask;
    if (hasmask) {
        mask = bp::extract<escript::Data>(arg["mask"]);

        // Possible error: Mask is not scalar
        if (mask.getDataPointRank() != 0) {
            throw DataException("getNumpy: mask must be scalar.");
        }

        keys.remove("mask");
        numdata--;
    }

    // Possible error: User forgot to pass some data
    if (numdata < 1) {
        throw DataException("getNumpy: no data to save specified.");
    }

    // Initialise vectors
    std::vector<int> step(numdata);
    std::vector<std::string> names(numdata);
    std::vector<Data> data(numdata);
    std::vector<int> fstypes(numdata+int(hasmask)); // FunctionSpace types for each data for interpolation

    // We need to interpret the samples correctly even if they are different
    // types. For this reason, we should iterate over samples...
    for (int i = 0; i < numdata; ++i) {
        names[i] = bp::extract<std::string>(keys[i]);
        data[i] = bp::extract<escript::Data>(arg[keys[i]]);
        fstypes[i] = data[i].getFunctionSpace().getTypeCode();
        step[i] = (data[i].actsExpanded() ? DataTypes::noValues(data[i].getDataPointShape()) : 0);

        // Possible error: The data are on different domains
        if (i > 0) {
            if (data[i].getDomain()!=data[i-1].getDomain()) {
                throw DataException("getNumpy: all data must be on the same domain.");
            }
        }
    }

    // Check for complex data
    bool have_complex = false;
    for (int i = 0; i < numdata; ++i){
        if(data[i].isComplex()){
            have_complex = true;
            break;
        }
    }

    if (hasmask) {
        //Possible error: Mask domain is different to the data domain
        if (mask.getDomain() != data[0].getDomain())
            throw DataException("getNumpy: mask domain must be the same as the data domain.");

        // Get the mask fs code
        fstypes[numdata] = mask.getFunctionSpace().getTypeCode();
    }

    // Possible error: Functionspaces are incompatible so interpolation isn't possible
    int bestfnspace = 0;
    if (!data[0].getDomain()->commonFunctionSpace(fstypes, bestfnspace)) {
        throw DataException("getNumpy: FunctionSpaces of data are incompatible");
    }

    // Interpolate the data onto the same function space
    FunctionSpace best(data[0].getDomain(), bestfnspace);
    for (int i = 0; i < numdata; ++i) {
        data[i] = data[i].interpolate(best);
    }
    if (hasmask)
        mask = mask.interpolate(best);

    // This are needed below
    const DataTypes::real_t onlyreal = 0;
    const DataTypes::cplx_t gotcomplex = 0.0;

    // Work out how big the ndarrays have to be
    int arraylength = 0;
    int numsamples = data[0].getNumSamples();
    int dpps = data[0].getNumDataPointsPerSample();
    if(hasmask){
        for(int i = 0; i < numsamples * dpps; i++){
            arraylength += (bool) *mask.getSampleDataRO(i, onlyreal);
        }
    } else {
        arraylength = dpps * numsamples;
    }

    //Work out how many rows each array should have
    std::vector<int> spaces(numdata);
    signed int total = 0;
    spaces[0] = 0;
    for(int i = 0; i < numdata; i++){
        total += data[i].getShapeProduct();
        spaces[i+1] = total;
    }

    // Initialise the numpy ndarray
    bp::tuple arrayshape = bp::make_tuple(total, arraylength);
    bp::numpy::dtype datatype = bp::numpy::dtype::get_builtin<double>();
    if(have_complex){
        datatype = bp::numpy::dtype::get_builtin<std::complex<double>>();
    }
    bp::numpy::ndarray dataArray = bp::numpy::zeros(arrayshape, datatype);

    // Initialise variables
    const double* masksample = NULL;
    bool expandedmask = false;              // does the mask act expanded?
    bool wantrow = true;                    // do we output this row?
    if (hasmask && mask.actsExpanded()) {
        expandedmask=true;
    }

    // Make everything more complicated.
    if(have_complex){
        for(int i = 0; i < numdata; i++){
            if(!data[i].isComplex()){
                data[i].complicate();
            }
        }
    }

    int error = 0;
    // int maskcounter[numdata];
    std::vector<int> maskcounter(numdata);
    for(int i = 0; i < numdata; i++) maskcounter[i]=0;
    std::string localmsg;
    try {
        std::vector<int> offset(numdata);
        std::vector<const DataTypes::real_t*> samplesR(numdata);
        std::vector<const DataTypes::cplx_t*> samplesC(numdata);

        for (int i = 0; i < numsamples; ++i) {

            // If this MPI process does not own the sample then continue
            if (!best.ownSample(i)) {
                continue;
            }

            // Get the sample data
            for (int d = 0; d < numdata; ++d) {
                if(data[d].isComplex()){
                    samplesC[d] = data[d].getSampleDataRO(i, gotcomplex);
                } else {
                    samplesR[d] = data[d].getSampleDataRO(i, onlyreal);
                }
            }

            // Work out if we want to get this row
            wantrow = true;
            if (hasmask) {
                masksample = mask.getSampleDataRO(i, onlyreal);
                if (!expandedmask) {
                    // mask controls whole sample
                    if (masksample[0] <= 0) {
                        wantrow = false;
                    }
                }
            }

            // Loop over data points in the sample
            for (int j = 0; j < dpps; ++j) {

                // now we need to check if this point is masked off
                if (expandedmask) {
                    // masks are scalar so the relevant value is at [j]
                    wantrow = (masksample[j]>0);
                }

                // If we want the row then add it to the array
                if (wantrow) {
                    for (int d = 0; d < numdata; ++d) {
                        if(have_complex){
                            DataTypes::pointToNumpyArray(dataArray, samplesC[d],
                                data[d].getDataPointShape(), offset[d], spaces[d], hasmask ? maskcounter[d] : i+j*numsamples);
                        } else {

                            DataTypes::pointToNumpyArray(dataArray, samplesR[d],
                                data[d].getDataPointShape(), offset[d], spaces[d], hasmask ? maskcounter[d] : i+j*numsamples);
                        }

                        offset[d] += step[d];
                        maskcounter[d]++;
                    }
                }
            }

            // Debug output
            // std::cout << "Current array:\n" << bp::extract<char const *>(bp::str(dataArray)) << std::endl;

            // Reset
            for (int d = 0; d < numdata; d++) {
                offset[d] = 0;
            }

        }
    } catch (EsysException* e) {
        error=1;
        if (data[0].getDomain()->getMPISize()==1) {
            throw;
        } else {
            localmsg=e->what();
        }
    } catch (...) {
        error=1;
        if (data[0].getDomain()->getMPISize()==1) {
            throw;
        }
    }

#ifdef ESYS_MPI
    MPI_Comm com = data[0].get_MPIComm();
    int rerror = 0;
    MPI_Allreduce(&error, &rerror, 1, MPI_INT, MPI_MAX, com);
    error = rerror;
#endif

    if (error) {
        if (localmsg.empty()) {
            throw DataException("getNumpy: error building output");
        }
        throw DataException(std::string("getNumpy:")+localmsg);
    }

    // MPI Barrier
    data[0].getDomain()->MPIBarrier();

    // Other, unknown errors
    if (error)
        throw DataException("getNumpy: Unknown error.");

    // Put everything into a list to return to python
    bp::list answer;
    for(int i = 0; i < numdata; i++){
        bp::numpy::ndarray temp = bp::numpy::zeros(bp::make_tuple(data[i].getShapeProduct(), arraylength), datatype);
        for(int j = 0; j < data[i].getShapeProduct(); j++){
            temp[j] = dataArray[spaces[i]+j];
        }
        // std::cout << "Temp array:\n" << bp::extract<char const *>(bp::str(temp)) << std::endl;
        answer.extend(names[i]);
        answer.extend(temp);
    }

    // Print out the ndarray to the console - used during debugging
    // std::cout << "Finished array:\n" << bp::extract<char const *>(bp::str(dataArray)) << std::endl;

    return answer;
}
#else
void getNumpy(bp::dict arg){
    throw DataException("getNumpy: Error - Please recompile escripts with the boost numpy library");
}
#endif

#ifdef ESYS_HAVE_BOOST_NUMPY
boost::python::numpy::ndarray convertToNumpy(escript::Data data)
{
    // Initialise boost numpy
    // Py_Initialize();
    boost::python::numpy::initialize();

    // Check to see if we have complex data
    bool have_complex = data.isComplex();

    // Work out how many data points there are
    int numsamples = data.getNumSamples();
    int dpps = data.getNumDataPointsPerSample();
    // int numdata = data.getNumDataPoints();

    // Work out the data point shape
    std::vector<int> shape = data.getDataPointShape();
    if(shape.size() == 0){ // If we have scalar data, the shape will be ()
        shape.push_back(1);
    }

    //Work out how many rows each array should have
    // int spaces = data.getShapeProduct();
        
    // We need to interpret the samples correctly even if they are different
    // types. For this reason, we should iterate over samples...
    int step = (data.actsExpanded() ? DataTypes::noValues(data.getDataPointShape()) : 0);

    // Initialise the numpy ndarray
    int arraylength = dpps * numsamples;
    bp::tuple arrayshape = bp::make_tuple(shape[0], arraylength);
    bp::numpy::dtype datatype = bp::numpy::dtype::get_builtin<double>();
    if(have_complex){
        datatype = bp::numpy::dtype::get_builtin<std::complex<double>>();
    }
    bp::numpy::ndarray dataArray = bp::numpy::zeros(arrayshape, datatype);

    // Initialise variables
    const DataTypes::cplx_t* samplesC=nullptr;
    const DataTypes::real_t* samplesR=nullptr;

    // This is needed below in getSampleDataRO
    const DataTypes::cplx_t onlycomplex=(0);
    const DataTypes::real_t onlyreal=0;
    
    int offset=0;
    for (int i=0; i<numsamples; ++i) {
        if(have_complex)
            samplesC = data.getSampleDataRO(i, onlycomplex);
        else
            samplesR = data.getSampleDataRO(i, onlyreal);

        for (int j=0; j<dpps; ++j) {
            if(have_complex){
                DataTypes::pointToNumpyArray(dataArray, samplesC,
                    data.getDataPointShape(), offset, 0, i+j*numsamples);
            } else {
                DataTypes::pointToNumpyArray(dataArray, samplesR,
                    data.getDataPointShape(), offset, 0, i+j*numsamples);
            }
            offset += step;
        }
        offset=0;
    }

    // Print out the ndarray to the console - used during debugging
    // std::cout << "Finished array:\n" << bp::extract<char const *>(bp::str(dataArray)) << std::endl;

    return dataArray;
}
#else
void convertToNumpy(escript::Data data){
    throw DataException("getNumpy: Error - Please recompile escripts with the boost numpy library");
}
#endif

void resolveGroup(bp::object obj)
{
    int len=0;
    try {
        len=bp::extract<int>(obj.attr("__len__")());
    }
    catch(...)
    {
        // tell python the error isn't there anymore
        PyErr_Clear();
        throw DataException("resolveGroup: sequence object expected.");
    }
    std::vector<DataLazy*> dats;
    std::vector<Data*> dp;
    for (int i=0; i<len; ++i) {
        Data* p=0;
        try {
            p = bp::extract<Data*>(obj[i]);
        } catch(...) {
            PyErr_Clear();
            throw DataException("resolveGroup: only accepts Data objects.");
        }
        if (p->isLazy()) {
            dats.push_back(dynamic_cast<DataLazy*>(p->borrowData()));
            dp.push_back(p);
        }
    }
    if (!dats.empty()) {
        dats[0]->resolveGroupWorker(dats);
    }

    // all the data will be identities now but still lazy
    // convert it to ready
    for (int i=dp.size()-1; i>=0; --i)
        dp[i]->resolve();
}

#ifdef ESYS_HAVE_NETCDF4
bool openNcFile(netCDF::NcFile& file, const std::string& filename)
{
    try {
        file.open(filename, netCDF::NcFile::read);
        return !file.isNull();
    } catch (...) {
        return false;
    }
}
#endif

} // end of namespace

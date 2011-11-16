
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <ripley/RipleyDomain.h>
#include <ripley/IndexList.h>
#include <ripley/Util.h>
extern "C" {
#include "esysUtils/blocktimer.h"
}
#include <escript/Data.h>
#include <escript/DataFactory.h>

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

#ifdef USE_PARMETIS
#include <parmetis.h>
#endif

#include <boost/python/import.hpp>
#include <boost/python/tuple.hpp>
#include <iomanip>

using namespace std;

namespace ripley {

//
// define the static constants
RipleyDomain::FunctionSpaceNamesMapType RipleyDomain::m_functionSpaceTypeNames;

RipleyDomain::RipleyDomain(const string& name, dim_t numDim, Esys_MPIInfo* mpiInfo) :
    m_name(name)
{
    setFunctionSpaceTypeNames();
    m_mpiInfo = Esys_MPIInfo_getReference(mpiInfo);
    m_nodes.reset(new NodeFile(numDim, m_mpiInfo));
    if (numDim==3) {
        m_elements.reset(new ElementFile(Hex8, m_mpiInfo));
        m_faceElements.reset(new ElementFile(Rec4, m_mpiInfo));
    } else if (numDim==2) {
        m_elements.reset(new ElementFile(Rec4, m_mpiInfo));
        m_faceElements.reset(new ElementFile(Line2, m_mpiInfo));
    } else
        throw RipleyException("RipleyDomain: Invalid number of dimensions");
    m_points.reset(new ElementFile(Point1, m_mpiInfo));
    m_fullFullPattern = NULL;
    m_fullReducedPattern = NULL;
    m_reducedFullPattern = NULL;
    m_reducedReducedPattern = NULL;
}

RipleyDomain::~RipleyDomain()
{
    Esys_MPIInfo_free(m_mpiInfo);
}

int RipleyDomain::getMPISize() const
{
    return m_mpiInfo->size;
}

int RipleyDomain::getMPIRank() const
{
    return m_mpiInfo->rank;
}

void RipleyDomain::MPIBarrier() const
{
#ifdef ESYS_MPI
    MPI_Barrier(m_mpiInfo->comm);
#endif
}

bool RipleyDomain::onMasterProcessor() const
{
    return m_mpiInfo->rank == 0;
}

#ifdef ESYS_MPI
MPI_Comm
#else
unsigned int
#endif
RipleyDomain::getMPIComm() const
{
#ifdef ESYS_MPI
    return m_mpiInfo->comm;
#else
    return 0;
#endif
}

void RipleyDomain::write(const string& filename) const
{
    if (m_mpiInfo->size > 1)
        throw RipleyException("write: only single processor runs are supported.");

    // open file
    ofstream f(filename.c_str());
    if (!f.is_open()) {
        stringstream msg;
        msg << "write: Opening file " << filename << " for writing failed.";
        throw RipleyException(msg.str());
    }

    // write header
    f << m_name << endl;

    f.precision(15);
    f.setf(ios::scientific, ios::floatfield);

    // write nodes
    dim_t numDim = getDim();
    f << numDim << "D-Nodes " << m_nodes->getNumNodes() << endl;
    for (dim_t i = 0; i < m_nodes->getNumNodes(); i++) {
        f << m_nodes->getIdVector()[i] << " "
            << m_nodes->getGlobalDegreesOfFreedom()[i] << " "
            << m_nodes->getTagVector()[i];
        for (dim_t j = 0; j < numDim; j++)
            f << " " << setw(20) << m_nodes->getCoordinates()[INDEX2(j, i, numDim)];
        f << endl;
    }

    // write all element types
    const ElementFile_ptr eFiles[3] = { m_elements, m_faceElements, m_points };

    for (size_t k=0; k<3; k++) {
        f << eFiles[k]->getName() << " " << eFiles[k]->getNumElements() << endl;
        dim_t NN = eFiles[k]->getNumNodes();
        for (dim_t i = 0; i < eFiles[k]->getNumElements(); i++) {
            f << eFiles[k]->getIdVector()[i] << " "
                << eFiles[k]->getTagVector()[i];
            for (dim_t j = 0; j < NN; j++)
                f << " " << m_nodes->getIdVector()[eFiles[k]->getNodes()[INDEX2(j, i, NN)]];
            f << endl;
        }
    }

    // write tags
    if (m_tagMap.size() > 0) {
        f << "Tags" << endl;
        TagMap::const_iterator it;
        for (it=m_tagMap.begin(); it!=m_tagMap.end(); it++) {
            f << it->first << " " << it->second << endl;
        }
    }
    f.close();
#ifdef Ripley_TRACE
    cout << "mesh " << m_name << " has been written to file " << filename << endl;
#endif
}

void RipleyDomain::Print_Mesh_Info(const bool full) const
{
    cout << "Ripley_PrintMesh_Info running on CPU " <<
        m_mpiInfo->rank << " of " << m_mpiInfo->size << endl;
    cout << "\tMesh name '" << m_name << "'" << endl;

    // write nodes
    int numDim = getDim();
    cout << "\tNodes: " << numDim << "D-Nodes " << m_nodes->getNumNodes() << endl;
    if (full) {
        cout << "\t     Id   Tag  gDOF   gNI grDfI  grNI:  Coordinates" << endl;
        cout.precision(15);
        cout.setf(ios::scientific, ios::floatfield);
        for (int i = 0; i < m_nodes->getNumNodes(); i++) {
            cout << "\t  " << setw(5) << m_nodes->getIdVector()[i] << " "
               << setw(5) << m_nodes->getTagVector()[i] << " "
               << setw(5) << m_nodes->getGlobalDegreesOfFreedom()[i] << " "
               << setw(5) << m_nodes->getGlobalNodesIndex()[i] << " "
               << setw(5) << m_nodes->getGlobalReducedDOFIndex()[i] << " "
               << setw(5) << m_nodes->getGlobalReducedNodesIndex()[i] << ": ";
            for (int j = 0; j < numDim; j++)
                cout << " " << m_nodes->getCoordinates()[INDEX2(j, i, numDim)];
            cout << endl;
        }
    }

    // write all element types
    const char *eNames[3] = { "Elements", "Face elements", "Points" };
    const ElementFile_ptr eFiles[3] = { m_elements, m_faceElements, m_points };

    for (size_t k=0; k<3; k++) {
        int mine = 0, overlap = 0;
        for (dim_t i = 0; i < eFiles[k]->getNumElements(); i++) {
            if (eFiles[k]->getOwnerVector()[i] == m_mpiInfo->rank)
                mine++;
            else
                overlap++;
        }
        cout << "\t" << eNames[k] << ": " << eFiles[k]->getName()
            << " " << eFiles[k]->getNumElements()
            << " (TypeId=" << eFiles[k]->getTypeId() << ") owner="
            << mine << " overlap=" << overlap << endl;
        int NN = eFiles[k]->getNumNodes();
        if (full) {
            cout << "\t     Id   Tag Owner Color:  Nodes" << endl;
            for (int i = 0; i < eFiles[k]->getNumElements(); i++) {
                cout << "\t  " << setw(5) << eFiles[k]->getIdVector()[i] << " "
                   << setw(5) << eFiles[k]->getTagVector()[i] << " "
                   << setw(5) << eFiles[k]->getOwnerVector()[i] << " "
                   << setw(5) << eFiles[k]->getColorVector()[i] << ": ";
                for (int j = 0; j < NN; j++)
                    cout << " " << setw(5) << m_nodes->getIdVector()[eFiles[k]->getNodes()[INDEX2(j, i, NN)]];
                cout << endl;
            }
        }
    }


    // write tags
    if (m_tagMap.size() > 0) {
        cout << "\tTags:" << endl;
        TagMap::const_iterator it;
        for (it=m_tagMap.begin(); it!=m_tagMap.end(); it++) {
            cout << "\t  " << setw(5) << it->second << " "
                << it->first << endl;
        }
    }
}

void RipleyDomain::dump(const string& fileName) const
{
#ifdef USE_NETCDF
    NcDim* ncdim = NULL;
    int mpi_size = m_mpiInfo->size;
    int mpi_rank = m_mpiInfo->rank;
    int numDim = getDim();

/* Incoming token indicates it's my turn to write */
#ifdef ESYS_MPI
    MPI_Status status;
    if (mpi_rank>0) {
        int dummy;
        MPI_Recv(&dummy, 0, MPI_INT, mpi_rank-1, 81800, m_mpiInfo->comm, &status);
    }
#endif

    char *newFileName = Esys_MPI_appendRankToFileName(fileName.c_str(),
                                                      mpi_size, mpi_rank);

    // NetCDF error handler
    NcError err(NcError::verbose_nonfatal);
    // Create the file.
    NcFile dataFile(newFileName, NcFile::Replace);
    const string msgPrefix("dump: netCDF operation failed - ");

    // check if writing was successful
    if (!dataFile.is_valid())
        throw RipleyException(msgPrefix+"Open file for output");

    const size_t numTags = m_tagMap.size();

    if (numTags > 0)
        if (! (ncdim = dataFile.add_dim("dim_Tags", numTags)) )
            throw RipleyException(msgPrefix+"add_dim(dim_Tags)");

    // Attributes: MPI size, MPI rank, Name, order, reduced_order
    if (!dataFile.add_att("mpi_size", mpi_size))
        throw RipleyException(msgPrefix+"add_att(mpi_size)");
    if (!dataFile.add_att("mpi_rank", mpi_rank))
        throw RipleyException(msgPrefix+"add_att(mpi_rank)");
    if (!dataFile.add_att("Name", m_name.c_str()))
        throw RipleyException(msgPrefix+"add_att(Name)");
    if (!dataFile.add_att("numDim", numDim))
        throw RipleyException(msgPrefix+"add_att(numDim)");
    if (!dataFile.add_att("order", 1))
        throw RipleyException(msgPrefix+"add_att(order)");
    if (!dataFile.add_att("reduced_order", 1))
        throw RipleyException(msgPrefix+"add_att(reduced_order)");
    if (!dataFile.add_att("num_Tags", static_cast<int>(numTags)))
        throw RipleyException(msgPrefix+"add_att(num_Tags)");

    m_nodes->dumpToNetCDF(dataFile);
    m_elements->dumpToNetCDF(dataFile, "Elements");
    m_faceElements->dumpToNetCDF(dataFile, "FaceElements");
    m_points->dumpToNetCDF(dataFile, "Points");

    // // // // // TagMap // // // // //
    if (numTags > 0) {

        // Copy tag keys into temp array
        vector<int> Tags_keys;
        TagMap::const_iterator it;
        for (it=m_tagMap.begin(); it!=m_tagMap.end(); it++) {
            Tags_keys.push_back(it->second);
        }

        // Tags_keys
        NcVar *ncVar;
        if (! (ncVar = dataFile.add_var("Tags_keys", ncInt, ncdim)) )
            throw RipleyException(msgPrefix+"add_var(Tags_keys)");
        if (!ncVar->put(&Tags_keys[0], numTags))
            throw RipleyException(msgPrefix+"put(Tags_keys)");

        // Tags_names_*
        // This is an array of strings, it should be stored as an array but
        // instead I have hacked in one attribute per string because the NetCDF
        // manual doesn't tell how to do an array of strings
        int i = 0;
        for (it=m_tagMap.begin(); it!=m_tagMap.end(); it++) {
            stringstream name;
            name << "Tags_name_" << i;
            if (!dataFile.add_att(name.str().c_str(), it->first.c_str()) )
                throw RipleyException(msgPrefix+"add_att(Tags_names_XX)");
            i++;
        }
    }

    // Send token to next MPI process so he can take his turn
#ifdef ESYS_MPI
    if (mpi_rank<mpi_size-1) {
        int dummy = 0;
        MPI_Send(&dummy, 0, MPI_INT, mpi_rank+1, 81800, m_mpiInfo->comm);
    }
#endif

    // netCDF file is closed by destructor of NcFile object
#else
    throw RipleyException("dump: not configured with netCDF. Please contact your installation manager.");
#endif  /* USE_NETCDF */
}

string RipleyDomain::getDescription() const
{
    return "RipleyMesh";
}

string RipleyDomain::functionSpaceTypeAsString(int functionSpaceType) const
{
    FunctionSpaceNamesMapType::iterator loc;
    loc=m_functionSpaceTypeNames.find(functionSpaceType);
    if (loc!=m_functionSpaceTypeNames.end())
        return loc->second;

    return "Invalid function space type code";
}

bool RipleyDomain::isValidFunctionSpaceType(int functionSpaceType) const
{
    FunctionSpaceNamesMapType::iterator loc;
    loc=m_functionSpaceTypeNames.find(functionSpaceType);
    return (loc!=m_functionSpaceTypeNames.end());
}

void RipleyDomain::setFunctionSpaceTypeNames()
{
    if (m_functionSpaceTypeNames.size() > 0)
        return;

    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(DegreesOfFreedom,"Ripley_DegreesOfFreedom"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(ReducedDegreesOfFreedom,"Ripley_ReducedDegreesOfFreedom"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(Nodes,"Ripley_Nodes"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(ReducedNodes,"Ripley_Reduced_Nodes"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(Elements,"Ripley_Elements"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(ReducedElements,"Ripley_Reduced_Elements"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(FaceElements,"Ripley_Face_Elements"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(ReducedFaceElements,"Ripley_Reduced_Face_Elements"));
    m_functionSpaceTypeNames.insert(FunctionSpaceNamesMapType::value_type(Points,"Ripley_Points"));
}

int RipleyDomain::getContinuousFunctionCode() const
{
    return Nodes;
}

int RipleyDomain::getReducedContinuousFunctionCode() const
{
    return ReducedNodes;
}

int RipleyDomain::getFunctionCode() const
{
    return Elements;
}

int RipleyDomain::getReducedFunctionCode() const
{
    return ReducedElements;
}

int RipleyDomain::getFunctionOnBoundaryCode() const
{
    return FaceElements;
}

int RipleyDomain::getReducedFunctionOnBoundaryCode() const
{
    return ReducedFaceElements;
}

int RipleyDomain::getFunctionOnContactZeroCode() const
{
    throw RipleyException("Ripley does not support contact elements");
}

int RipleyDomain::getReducedFunctionOnContactZeroCode() const
{
    throw RipleyException("Ripley does not support contact elements");
}

int RipleyDomain::getFunctionOnContactOneCode() const
{
    throw RipleyException("Ripley does not support contact elements");
}

int RipleyDomain::getReducedFunctionOnContactOneCode() const
{
    throw RipleyException("Ripley does not support contact elements");
}

int RipleyDomain::getSolutionCode() const
{
    return DegreesOfFreedom;
}

int RipleyDomain::getReducedSolutionCode() const
{
    return ReducedDegreesOfFreedom;
}

int RipleyDomain::getDiracDeltaFunctionsCode() const
{
    return Points;
}

//
// returns the spatial dimension of the mesh
//
int RipleyDomain::getDim() const
{
    return m_nodes->getNumDim();
}

//
// returns the number of data points summed across all MPI processes
//
int RipleyDomain::getNumDataPointsGlobal() const
{
    return m_nodes->getGlobalNumNodes();
}

//
// returns the number of data points per sample and the number of samples
// needed to represent data on a parts of the mesh.
//
pair<int,int> RipleyDomain::getDataShape(int functionSpaceCode) const
{
    int numDataPointsPerSample=0;
    int numSamples=0;
    switch (functionSpaceCode) {
        case Nodes:
            numDataPointsPerSample=1;
            numSamples=m_nodes->getNumNodes();
            break;
        case ReducedNodes:
            numDataPointsPerSample=1;
            numSamples=m_nodes->getNumReducedNodes();
            break;
        case Elements:
            numSamples=m_elements->getNumElements();
            numDataPointsPerSample=m_elements->getNumLocalDim()+1;
            break;
        case ReducedElements:
            numSamples=m_elements->getNumElements();
            numDataPointsPerSample=(m_elements->getNumLocalDim()==0)?0:1;
            break;
        case FaceElements:
            numDataPointsPerSample=m_faceElements->getNumLocalDim()+1;
            numSamples=m_faceElements->getNumElements();
            break;
        case ReducedFaceElements:
            numDataPointsPerSample=(m_faceElements->getNumLocalDim()==0)?0:1;
            numSamples=m_faceElements->getNumElements();
            break;
        case Points:
            numDataPointsPerSample=1;
            numSamples=m_points->getNumElements();
            break;
        case DegreesOfFreedom:
            numDataPointsPerSample=1;
            numSamples=m_nodes->getNumDegreesOfFreedom();
            break;
        case ReducedDegreesOfFreedom:
            numDataPointsPerSample=1;
            numSamples=m_nodes->getNumReducedDegreesOfFreedom();
            break;
        default:
            stringstream temp;
            temp << "Invalid function space type: " << functionSpaceCode << " for domain: " << getDescription();
            throw RipleyException(temp.str());
    }
    return pair<int,int>(numDataPointsPerSample,numSamples);
}

//
// adds linear PDE of second order into a given stiffness matrix and right hand side
//
void RipleyDomain::addPDEToSystem(
        escript::AbstractSystemMatrix& mat, escript::Data& rhs,
        const escript::Data& A, const escript::Data& B, const escript::Data& C,
        const escript::Data& D, const escript::Data& X, const escript::Data& Y,
        const escript::Data& d, const escript::Data& y,
        const escript::Data& d_contact, const escript::Data& y_contact,
        const escript::Data& d_dirac,const escript::Data& y_dirac) const
{
    if (!d_contact.isEmpty() || !y_contact.isEmpty())
        throw RipleyException("Ripley does not support contact elements");

    SystemMatrixAdapter* smat=dynamic_cast<SystemMatrixAdapter*>(&mat);
    if (smat==NULL)
        throw RipleyException("Ripley only accepts its own system matrices");
/*
    escriptDataC _rhs=rhs.getDataC();
    escriptDataC _A =A.getDataC();
    escriptDataC _B=B.getDataC();
    escriptDataC _C=C.getDataC();
    escriptDataC _D=D.getDataC();
    escriptDataC _X=X.getDataC();
    escriptDataC _Y=Y.getDataC();
    escriptDataC _d=d.getDataC();
    escriptDataC _y=y.getDataC();
    escriptDataC _d_dirac=d_dirac.getDataC();
    escriptDataC _y_dirac=y_dirac.getDataC();

    Ripley_Mesh* mesh=m_ripleyMesh.get();
    Ripley_Assemble_PDE(m_nodes,mesh->Elements,smat->getPaso_SystemMatrix(), &_rhs, &_A, &_B, &_C, &_D, &_X, &_Y);
    checkPasoError();

    Ripley_Assemble_PDE(m_nodes,mesh->FaceElements, smat->getPaso_SystemMatrix(), &_rhs, 0, 0, 0, &_d, 0, &_y);
    checkPasoError();

    Ripley_Assemble_PDE(m_nodes,mesh->Points, smat->getPaso_SystemMatrix(), &_rhs, 0, 0, 0, &_d_dirac, 0, &_y_dirac);
    checkPasoError();
*/
}

void RipleyDomain::addPDEToLumpedSystem(escript::Data& mat,
        const escript::Data& D, const escript::Data& d,
        const escript::Data& d_dirac, const bool useHRZ) const
{
/*
    escriptDataC _mat=mat.getDataC();
    escriptDataC _D=D.getDataC();
    escriptDataC _d=d.getDataC();
    escriptDataC _d_dirac=d_dirac.getDataC();

    Ripley_Mesh* mesh=m_ripleyMesh.get();
    Ripley_Assemble_LumpedSystem(m_nodes, mesh->Elements, &_mat, &_D, useHRZ);
    checkPasoError();

    Ripley_Assemble_LumpedSystem(m_nodes, mesh->FaceElements, &_mat, &_d, useHRZ);
    checkPasoError();

    Ripley_Assemble_LumpedSystem(m_nodes, mesh->FaceElements, &_mat, &_d_dirac, useHRZ);
    checkPasoError();
*/
}


//
// adds linear PDE of second order into the right hand side only
//
void RipleyDomain::addPDEToRHS(escript::Data& rhs, const escript::Data& X,
        const escript::Data& Y, const escript::Data& y,
        const escript::Data& y_contact, const escript::Data& y_dirac) const
{
    if (!y_contact.isEmpty())
        throw RipleyException("Ripley does not support contact elements");
/*
    Ripley_Mesh* mesh=m_ripleyMesh.get();

    escriptDataC _rhs=rhs.getDataC();
    escriptDataC _X=X.getDataC();
    escriptDataC _Y=Y.getDataC();
    escriptDataC _y=y.getDataC();
    escriptDataC _y_dirac=y_dirac.getDataC();
    Ripley_Assemble_PDE(m_nodes, mesh->Elements, 0, &_rhs, 0, 0, 0, 0, &_X, &_Y );
    checkPasoError();

    Ripley_Assemble_PDE(m_nodes, mesh->FaceElements, 0, &_rhs, 0, 0, 0, 0, 0, &_y );
    checkPasoError();

    Ripley_Assemble_PDE(m_nodes, mesh->Points, 0, &_rhs, 0, 0, 0, 0, 0, &_y_dirac );
    checkPasoError();
*/
}

//
// adds PDE of second order into a transport problem
//
void RipleyDomain::addPDEToTransportProblem(
        escript::AbstractTransportProblem& tp,
        escript::Data& source, const escript::Data& M,
        const escript::Data& A, const escript::Data& B, const escript::Data& C,
        const escript::Data& D, const escript::Data& X, const escript::Data& Y,
        const escript::Data& d, const escript::Data& y,
        const escript::Data& d_contact, const escript::Data& y_contact,
        const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
    if (!d_contact.isEmpty() || !y_contact.isEmpty())
        throw RipleyException("Ripley does not support contact elements");
    TransportProblemAdapter* tpa=dynamic_cast<TransportProblemAdapter*>(&tp);
    if (tpa==NULL)
        throw RipleyException("Ripley only accepts its own transport problems");
/*
    DataTypes::ShapeType shape;
    source.expand();
    escriptDataC _source=source.getDataC();
    escriptDataC _M=M.getDataC();
    escriptDataC _A=A.getDataC();
    escriptDataC _B=B.getDataC();
    escriptDataC _C=C.getDataC();
    escriptDataC _D=D.getDataC();
    escriptDataC _X=X.getDataC();
    escriptDataC _Y=Y.getDataC();
    escriptDataC _d=d.getDataC();
    escriptDataC _y=y.getDataC();
    escriptDataC _d_dirac=d_dirac.getDataC();
    escriptDataC _y_dirac=y_dirac.getDataC();

    Ripley_Mesh* mesh=m_ripleyMesh.get();
    Paso_TransportProblem* _tp = tpa->getPaso_TransportProblem();
    Ripley_Assemble_PDE(m_nodes, mesh->Elements,_tp->mass_matrix, &_source, 0, 0, 0, &_M, 0, 0);
    checkPasoError();

    Ripley_Assemble_PDE(m_nodes, mesh->Elements,_tp->transport_matrix, &_source, &_A, &_B, &_C, &_D, &_X, &_Y);
    checkPasoError();

    Ripley_Assemble_PDE(m_nodes, mesh->FaceElements, _tp->transport_matrix, &_source, 0, 0, 0, &_d, 0, &_y);
    checkPasoError();

    Ripley_Assemble_PDE(m_nodes, mesh->Points, _tp->transport_matrix, &_source, 0, 0, 0, &_d_dirac, 0, &_y_dirac);
    checkPasoError();
*/
}

//
// interpolates data between different function spaces
//
void RipleyDomain::interpolateOnDomain(escript::Data& target,
                                      const escript::Data& in) const
{
    const RipleyDomain& inDomain=dynamic_cast<const RipleyDomain&>(*(in.getFunctionSpace().getDomain()));
    const RipleyDomain& targetDomain=dynamic_cast<const RipleyDomain&>(*(target.getFunctionSpace().getDomain()));
    if (inDomain != *this)
        throw RipleyException("Illegal domain of interpolant");
    if (targetDomain != *this)
        throw RipleyException("Illegal domain of interpolation target");
/*
    Ripley_Mesh* mesh=m_ripleyMesh.get();
    escriptDataC _target=target.getDataC();
    escriptDataC _in=in.getDataC();

    switch (in.getFunctionSpace().getTypeCode()) {
        case Nodes:
            switch (target.getFunctionSpace().getTypeCode()) {
                case Nodes:
                case ReducedNodes:
                case DegreesOfFreedom:
                case ReducedDegreesOfFreedom:
                    Ripley_Assemble_CopyNodalData(m_nodes, &_target, &_in);
                    break;
                case Elements:
                case ReducedElements:
                    Ripley_Assemble_interpolate(m_nodes, mesh->Elements, &_in, &_target);
                    break;
                case FaceElements:
                case ReducedFaceElements:
                    Ripley_Assemble_interpolate(m_nodes, mesh->FaceElements, &_in, &_target);
                    break;
                case Points:
                    Ripley_Assemble_interpolate(m_nodes, mesh->Points, &_in, &_target);
                    break;
                default:
                    stringstream temp;
                    temp << "Interpolation on Domain: Ripley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
                    throw RipleyException(temp.str());
            }
            break;
        case ReducedNodes:
            switch (target.getFunctionSpace().getTypeCode()) {
                case Nodes:
                case ReducedNodes:
                case DegreesOfFreedom:
                case ReducedDegreesOfFreedom:
                    Ripley_Assemble_CopyNodalData(m_nodes,&_target,&_in);
                    break;
                case Elements:
                case ReducedElements:
                    Ripley_Assemble_interpolate(m_nodes,mesh->Elements,&_in,&_target);
                    break;
                case FaceElements:
                case ReducedFaceElements:
                    Ripley_Assemble_interpolate(m_nodes,mesh->FaceElements,&_in,&_target);
                    break;
                case Points:
                    Ripley_Assemble_interpolate(m_nodes,mesh->Points,&_in,&_target);
                    break;
                default:
                    stringstream temp;
                    temp << "Interpolation on Domain: Ripley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
                    throw RipleyException(temp.str());
            }
            break;
        case Elements:
            if (target.getFunctionSpace().getTypeCode()==Elements) {
                Ripley_Assemble_CopyElementData(mesh->Elements,&_target,&_in);
            } else if (target.getFunctionSpace().getTypeCode()==ReducedElements) {
                Ripley_Assemble_AverageElementData(mesh->Elements,&_target,&_in);
            } else {
                throw RipleyException("No interpolation with data on elements possible.");
            }
            break;
        case ReducedElements:
            if (target.getFunctionSpace().getTypeCode()==ReducedElements) {
                Ripley_Assemble_CopyElementData(mesh->Elements,&_target,&_in);
            } else {
                throw RipleyException("No interpolation with data on elements with reduced integration order possible.");
            }
            break;
        case FaceElements:
            if (target.getFunctionSpace().getTypeCode()==FaceElements) {
                Ripley_Assemble_CopyElementData(mesh->FaceElements,&_target,&_in);
            } else if (target.getFunctionSpace().getTypeCode()==ReducedFaceElements) {
                Ripley_Assemble_AverageElementData(mesh->FaceElements,&_target,&_in);
            } else {
                throw RipleyException("No interpolation with data on face elements possible.");
            }
            break;
        case ReducedFaceElements:
            if (target.getFunctionSpace().getTypeCode()==ReducedFaceElements) {
                Ripley_Assemble_CopyElementData(mesh->FaceElements,&_target,&_in);
            } else {
                throw RipleyException("No interpolation with data on face elements with reduced integration order possible.");
            }
            break;
        case Points:
            if (target.getFunctionSpace().getTypeCode()==Points) {
                Ripley_Assemble_CopyElementData(mesh->Points,&_target,&_in);
            } else {
                throw RipleyException("No interpolation with data on points possible.");
            }
            break;
        case DegreesOfFreedom:
            switch (target.getFunctionSpace().getTypeCode()) {
                case ReducedDegreesOfFreedom:
                case DegreesOfFreedom:
                    Ripley_Assemble_CopyNodalData(m_nodes,&_target,&_in);
                    break;

                case Nodes:
                case ReducedNodes:
                    if (getMPISize()>1) {
                        escript::Data temp=escript::Data(in);
                        temp.expand();
                        escriptDataC _in2 = temp.getDataC();
                        Ripley_Assemble_CopyNodalData(m_nodes,&_target,&_in2);
                    } else {
                        Ripley_Assemble_CopyNodalData(m_nodes,&_target,&_in);
                    }
                    break;
                case Elements:
                case ReducedElements:
                    if (getMPISize()>1) {
                        escript::Data temp=escript::Data(in, continuousFunction(*this) );
                        escriptDataC _in2 = temp.getDataC();
                        Ripley_Assemble_interpolate(m_nodes, mesh->Elements,&_in2,&_target);
                    } else {
                        Ripley_Assemble_interpolate(m_nodes, mesh->Elements,&_in,&_target);
                    }
                    break;
                case FaceElements:
                case ReducedFaceElements:
                    if (getMPISize()>1) {
                        escript::Data temp=escript::Data(in, continuousFunction(*this) );
                        escriptDataC _in2 = temp.getDataC();
                        Ripley_Assemble_interpolate(m_nodes, mesh->FaceElements,&_in2,&_target);
                    } else {
                        Ripley_Assemble_interpolate(m_nodes, mesh->FaceElements,&_in,&_target);
                    }
                    break;
                case Points:
                    if (getMPISize()>1) {
                        //escript::Data temp=escript::Data(in, continuousFunction(*this) );
                        //escriptDataC _in2 = temp.getDataC();
                    } else {
                        Ripley_Assemble_interpolate(m_nodes,mesh->Points,&_in,&_target);
                    }
                    break;
                default:
                    stringstream temp;
                    temp << "Interpolation On Domain: Ripley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
                    throw RipleyException(temp.str());
            }
            break;
        case ReducedDegreesOfFreedom:
            switch (target.getFunctionSpace().getTypeCode()) {
                case Nodes:
                    throw RipleyException("Ripley does not support interpolation from reduced degrees of freedom to mesh nodes.");
                    break;
                case ReducedNodes:
                    if (getMPISize()>1) {
                        escript::Data temp=escript::Data(in);
                        temp.expand();
                        escriptDataC _in2 = temp.getDataC();
                        Ripley_Assemble_CopyNodalData(m_nodes,&_target,&_in2);
                    } else {
                        Ripley_Assemble_CopyNodalData(m_nodes,&_target,&_in);
                    }
                    break;
                case DegreesOfFreedom:
                    throw RipleyException("Ripley does not support interpolation from reduced degrees of freedom to degrees of freedom");
                case ReducedDegreesOfFreedom:
                    Ripley_Assemble_CopyNodalData(m_nodes,&_target,&_in);
                    break;
                case Elements:
                case ReducedElements:
                    if (getMPISize()>1) {
                        escript::Data temp=escript::Data(in, reducedContinuousFunction(*this) );
                        escriptDataC _in2 = temp.getDataC();
                        Ripley_Assemble_interpolate(m_nodes, mesh->Elements,&_in2,&_target);
                    } else {
                        Ripley_Assemble_interpolate(m_nodes, mesh->Elements,&_in,&_target);
                    }
                    break;
                case FaceElements:
                case ReducedFaceElements:
                    if (getMPISize()>1) {
                        escript::Data temp=escript::Data(in, reducedContinuousFunction(*this) );
                        escriptDataC _in2 = temp.getDataC();
                        Ripley_Assemble_interpolate(m_nodes,mesh->FaceElements,&_in2,&_target);
                    } else {
                        Ripley_Assemble_interpolate(m_nodes,mesh->FaceElements,&_in,&_target);
                    }
                    break;
                case Points:
                    if (getMPISize()>1) {
                        escript::Data temp=escript::Data(in, reducedContinuousFunction(*this) );
                        escriptDataC _in2 = temp.getDataC();
                        Ripley_Assemble_interpolate(m_nodes, mesh->Points,&_in2,&_target);
                    } else {
                        Ripley_Assemble_interpolate(m_nodes, mesh->Points,&_in,&_target);
                    }
                    break;
                default:
                    stringstream temp;
                    temp << "Interpolation On Domain: Ripley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
                    throw RipleyException(temp.str());
            }
            break;
        default:
            stringstream temp;
            temp << "Interpolation On Domain: Ripley does not know anything about function space type %d" << in.getFunctionSpace().getTypeCode();
            throw RipleyException(temp.str());
    }
*/
    checkPasoError();
}

//
// copies the locations of sample points into arg:
//
void RipleyDomain::setToX(escript::Data& arg) const
{
    const RipleyDomain& argDomain=dynamic_cast<const RipleyDomain&>(*(arg.getFunctionSpace().getDomain()));
    if (argDomain != *this)
        throw RipleyException("setToX: Illegal domain of data point locations");

    // do we need to interpolate?
    if (arg.getFunctionSpace().getTypeCode()==Nodes) {
        m_nodes->assembleCoordinates(arg);
    } else {
        escript::Data tmp_data=Vector(0.0, continuousFunction(*this), true);
        m_nodes->assembleCoordinates(tmp_data);
        // this is then interpolated onto arg:
        interpolateOnDomain(arg, tmp_data);
    }
}

//
// returns the normal vectors at the location of data points as a Data object
//
void RipleyDomain::setToNormal(escript::Data& normal) const
{
/*
    const RipleyDomain& normalDomain=dynamic_cast<const RipleyDomain&>(*(normal.getFunctionSpace().getDomain()));
    if (normalDomain!=*this)
        throw RipleyException("Illegal domain of normal locations");
    Ripley_Mesh* mesh=m_ripleyMesh.get();
    escriptDataC _normal=normal.getDataC();
    switch(normal.getFunctionSpace().getTypeCode()) {
        case Nodes:
            throw RipleyException("Ripley does not support surface normal vectors for nodes");
        case ReducedNodes:
            throw RipleyException("Ripley does not support surface normal vectors for reduced nodes");
        case Elements:
            throw RipleyException("Ripley does not support surface normal vectors for elements");
        case ReducedElements:
            throw RipleyException("Ripley does not support surface normal vectors for elements with reduced integration order");
        case FaceElements:
            Ripley_Assemble_setNormal(m_nodes,mesh->FaceElements,&_normal);
            break;
        case ReducedFaceElements:
            Ripley_Assemble_setNormal(m_nodes,mesh->FaceElements,&_normal);
            break;
        case Points:
            throw RipleyException("Ripley does not support surface normal vectors for point elements");
        case DegreesOfFreedom:
            throw RipleyException("Ripley does not support surface normal vectors for degrees of freedom.");
        case ReducedDegreesOfFreedom:
            throw RipleyException("Ripley does not support surface normal vectors for reduced degrees of freedom.");
        default:
            stringstream temp;
            temp << "Normal Vectors: Ripley does not know anything about function space type " << normal.getFunctionSpace().getTypeCode();
            throw RipleyException(temp.str());
    }
*/
}

//
// interpolates data to other domain
//
void RipleyDomain::interpolateACross(escript::Data& target, const escript::Data& source) const
{
    escript::const_Domain_ptr targetDomain_p=target.getFunctionSpace().getDomain();
    const RipleyDomain* targetDomain=dynamic_cast<const RipleyDomain*>(targetDomain_p.get());
    if (targetDomain!=this)
        throw RipleyException("Illegal domain of interpolation target");

    throw RipleyException("Ripley does not allow interpolation across domains yet");
}

//
// calculates the integral of a function arg
//
void RipleyDomain::setToIntegrals(vector<double>& integrals, const escript::Data& arg) const
{
/*
    const RipleyDomain& argDomain=dynamic_cast<const RipleyDomain&>(*(arg.getFunctionSpace().getDomain()));
    if (argDomain!=*this)
        throw RipleyException("Illegal domain of integration kernel");

    double blocktimer_start = blocktimer_time();
    Ripley_Mesh* mesh=m_ripleyMesh.get();
    escriptDataC _temp;
    escript::Data temp;
    escriptDataC _arg=arg.getDataC();
    switch(arg.getFunctionSpace().getTypeCode()) {
        case Nodes:
        case ReducedNodes:
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
            temp=escript::Data( arg, escript::function(*this) );
            _temp=temp.getDataC();
            Ripley_Assemble_integrate(m_nodes,mesh->Elements,&_temp,&integrals[0]);
            break;
        case Elements:
        case ReducedElements:
            Ripley_Assemble_integrate(m_nodes,mesh->Elements,&_arg,&integrals[0]);
            break;
        case FaceElements:
        case ReducedFaceElements:
            Ripley_Assemble_integrate(m_nodes,mesh->FaceElements,&_arg,&integrals[0]);
            break;
        case Points:
            throw RipleyException("Integral of data on points is not supported.");
        default:
            stringstream temp;
            temp << "Integrals: Ripley does not know anything about function space type " << arg.getFunctionSpace().getTypeCode();
            throw RipleyException(temp.str());
    }
    blocktimer_increment("integrate()", blocktimer_start);
*/
}

//
// calculates the gradient of arg
//
void RipleyDomain::setToGradient(escript::Data& grad, const escript::Data& arg) const
{
/*
    const RipleyDomain& argDomain=dynamic_cast<const RipleyDomain&>(*(arg.getFunctionSpace().getDomain()));
    if (argDomain!=*this)
        throw RipleyException("Illegal domain of gradient argument");
    const RipleyDomain& gradDomain=dynamic_cast<const RipleyDomain&>(*(grad.getFunctionSpace().getDomain()));
    if (gradDomain!=*this)
        throw RipleyException("Illegal domain of gradient");

    Ripley_Mesh* mesh=m_ripleyMesh.get();
    escriptDataC _grad=grad.getDataC();
    escriptDataC nodeDataC;
    escript::Data temp;
    if (getMPISize()>1) {
        if( arg.getFunctionSpace().getTypeCode() == DegreesOfFreedom ) {
            temp=escript::Data(arg, continuousFunction(*this) );
            nodeDataC = temp.getDataC();
        } else if( arg.getFunctionSpace().getTypeCode() == ReducedDegreesOfFreedom ) {
            temp=escript::Data(arg, reducedContinuousFunction(*this) );
            nodeDataC = temp.getDataC();
        } else {
            nodeDataC = arg.getDataC();
        }
    } else {
        nodeDataC = arg.getDataC();
    }
    switch(grad.getFunctionSpace().getTypeCode()) {
        case Nodes:
            throw RipleyException("Gradient at nodes is not supported.");
        case ReducedNodes:
            throw RipleyException("Gradient at reduced nodes is not supported.");
        case Elements:
            Ripley_Assemble_gradient(m_nodes,mesh->Elements,&_grad,&nodeDataC);
            break;
        case ReducedElements:
            Ripley_Assemble_gradient(m_nodes,mesh->Elements,&_grad,&nodeDataC);
            break;
        case FaceElements:
            Ripley_Assemble_gradient(m_nodes,mesh->FaceElements,&_grad,&nodeDataC);
            break;
        case ReducedFaceElements:
            Ripley_Assemble_gradient(m_nodes,mesh->FaceElements,&_grad,&nodeDataC);
            break;
        case Points:
            throw RipleyException("Gradient at points is not supported");
        case DegreesOfFreedom:
            throw RipleyException("Gradient at degrees of freedom is not supported");
        case ReducedDegreesOfFreedom:
            throw RipleyException("Gradient at reduced degrees of freedom is not supported");
        default:
            stringstream temp;
            temp << "Gradient: Ripley does not know anything about function space type " << arg.getFunctionSpace().getTypeCode();
            throw RipleyException(temp.str());
    }
*/
}

//
// returns the size of elements
//
void RipleyDomain::setToSize(escript::Data& size) const
{
/*
    Ripley_Mesh* mesh=m_ripleyMesh.get();
    escriptDataC tmp=size.getDataC();
    switch(size.getFunctionSpace().getTypeCode()) {
        case Nodes:
            throw RipleyException("Size of nodes is not supported");
        case ReducedNodes:
            throw RipleyException("Size of reduced nodes is not supported");
        case Elements:
            Ripley_Assemble_getSize(m_nodes, mesh->Elements, &tmp);
            break;
        case ReducedElements:
            Ripley_Assemble_getSize(m_nodes, mesh->Elements, &tmp);
            break;
        case FaceElements:
            Ripley_Assemble_getSize(m_nodes, mesh->FaceElements, &tmp);
            break;
        case ReducedFaceElements:
            Ripley_Assemble_getSize(m_nodes, mesh->FaceElements, &tmp);
            break;
        case Points:
            throw RipleyException("Size of point elements is not supported");
        case DegreesOfFreedom:
            throw RipleyException("Size of degrees of freedom is not supported");
        case ReducedDegreesOfFreedom:
            throw RipleyException("Size of reduced degrees of freedom is not supported");
        default:
            stringstream temp;
            temp << "Element size: Ripley does not know anything about function space type " << size.getFunctionSpace().getTypeCode();
            throw RipleyException(temp.str());
    }
*/
}

//
// sets the location of nodes
//
void RipleyDomain::setNewX(const escript::Data& new_x)
{
    const RipleyDomain& newDomain=dynamic_cast<const RipleyDomain&>(*(new_x.getFunctionSpace().getDomain()));
    if (newDomain!=*this)
        throw RipleyException("Illegal domain of new node locations");

    if (new_x.getFunctionSpace() == continuousFunction(*this)) {
        m_nodes->setCoordinates(new_x);
    } else {
        escript::Data new_x_inter=escript::Data(new_x, continuousFunction(*this) );
        m_nodes->setCoordinates(new_x_inter);
    }
}

bool RipleyDomain::ownSample(int fs_code, index_t id) const
{
#ifdef ESYS_MPI
    index_t myFirstNode, myLastNode, k;
    if (fs_code == ReducedNodes) {
        myFirstNode = m_nodes->getFirstReducedNode();
        myLastNode = m_nodes->getLastReducedNode();
        k = m_nodes->getIndexForGlobalReducedNode(id);
    } else {
        myFirstNode = m_nodes->getFirstNode();
        myLastNode = m_nodes->getLastNode();
        k = m_nodes->getIndexForGlobalNode(id);
    }
    return static_cast<bool>((myFirstNode <= k) && (k < myLastNode));
#endif
    return true;
}


//
// creates a SystemMatrixAdapter stiffness matrix and initializes it with zeros
//
escript::ASM_ptr RipleyDomain::newSystemMatrix(const int row_blocksize,
        const escript::FunctionSpace& row_functionspace,
        const int column_blocksize,
        const escript::FunctionSpace& column_functionspace,
        const int type) const
{
    bool reduceRowOrder=false;
    bool reduceColOrder=false;
    // is the domain right?
    const RipleyDomain& row_domain=dynamic_cast<const RipleyDomain&>(*(row_functionspace.getDomain()));
    if (row_domain!=*this)
        throw RipleyException("Domain of row function space does not match the domain of matrix generator");
    const RipleyDomain& col_domain=dynamic_cast<const RipleyDomain&>(*(column_functionspace.getDomain()));
    if (col_domain!=*this)
        throw RipleyException("Domain of columnn function space does not match the domain of matrix generator");
    // is the function space type right
    if (row_functionspace.getTypeCode()==ReducedDegreesOfFreedom) {
        reduceRowOrder=true;
    } else if (row_functionspace.getTypeCode()!=DegreesOfFreedom)
        throw RipleyException("Illegal function space type for system matrix rows");
    if (column_functionspace.getTypeCode()==ReducedDegreesOfFreedom) {
        reduceColOrder=true;
    } else if (column_functionspace.getTypeCode()!=DegreesOfFreedom)
        throw RipleyException("Illegal function space type for system matrix columns");

    // generate matrix:
    Paso_SystemMatrixPattern* fsystemMatrixPattern=getPattern(reduceRowOrder, reduceColOrder);
    checkPasoError();
    Paso_SystemMatrix* fsystemMatrix = Paso_SystemMatrix_alloc(type,
            fsystemMatrixPattern, row_blocksize, column_blocksize, FALSE);
    checkPasoError();
    Paso_SystemMatrixPattern_free(fsystemMatrixPattern);
    SystemMatrixAdapter* sma = new SystemMatrixAdapter(fsystemMatrix,
            row_blocksize, row_functionspace, column_blocksize, column_functionspace);
    return escript::ASM_ptr(sma);
}

//
// creates a TransportProblemAdapter
//
escript::ATP_ptr RipleyDomain::newTransportProblem(const bool useBackwardEuler,
        const int blocksize, const escript::FunctionSpace& functionspace,
        const int type) const
{
    int reduceOrder=0;
    // is the domain right?
    const RipleyDomain& domain=dynamic_cast<const RipleyDomain&>(*(functionspace.getDomain()));
    if (domain!=*this)
        throw RipleyException("Domain of function space does not match the domain of transport problem generator");
    // is the function space type right
    if (functionspace.getTypeCode()==ReducedDegreesOfFreedom) {
        reduceOrder=1;
    } else if (functionspace.getTypeCode()!=DegreesOfFreedom)
        throw RipleyException("Illegal function space type for system matrix rows");

    // generate matrix
    Paso_SystemMatrixPattern* fsystemMatrixPattern = getPattern(reduceOrder, reduceOrder);
    Paso_TransportProblem* transportProblem = Paso_TransportProblem_alloc(
            useBackwardEuler, fsystemMatrixPattern, blocksize);
    checkPasoError();
    Paso_SystemMatrixPattern_free(fsystemMatrixPattern);
    escript::AbstractTransportProblem* atp=new TransportProblemAdapter(
            transportProblem, useBackwardEuler, blocksize, functionspace);
    return escript::ATP_ptr(atp);
}

// returns true if data on the function space is considered as being cell
// centered
bool RipleyDomain::isCellOriented(int functionSpaceCode) const
{
    switch(functionSpaceCode) {
        case Nodes:
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
            return false;
        case Elements:
        case FaceElements:
        case Points:
        case ReducedElements:
        case ReducedFaceElements:
            return true;
        default:
            stringstream temp;
            temp << "Ripley does not know anything about function space type " << functionSpaceCode;
            throw RipleyException(temp.str());
    }
    return false;
}

bool RipleyDomain::commonFunctionSpace(const vector<int>& fs, int& resultcode) const
{
   /* The idea is to use equivalence classes. [Types which can be interpolated back and forth]
    class 1: DOF <-> Nodes
    class 2: ReducedDOF <-> ReducedNodes
    class 3: Points
    class 4: Elements
    class 5: ReducedElements
    class 6: FaceElements
    class 7: ReducedFaceElements
    class 8: ContactElementZero <-> ContactElementOne
    class 9: ReducedContactElementZero <-> ReducedContactElementOne

   There is also a set of lines. Interpolation is possible down a line but not between lines.
   class 1 and 2 belong to all lines so aren't considered.
    line 0: class 3
    line 1: class 4,5
    line 2: class 6,7
    line 3: class 8,9

    For classes with multiple members (eg class 2) we have vars to record if
    there is at least one instance. e.g. hasnodes is true if we have at least
    one instance of Nodes.
    */
    if (fs.empty())
        return false;
    vector<int> hasclass(10);
    vector<int> hasline(4);
    bool hasnodes=false;
    bool hasrednodes=false;
    for (int i=0;i<fs.size();++i) {
        switch (fs[i]) {
            case Nodes: hasnodes=true; // no break is deliberate
            case DegreesOfFreedom:
                hasclass[1]=1;
                break;
            case ReducedNodes: hasrednodes=true; // no break is deliberate
            case ReducedDegreesOfFreedom:
                hasclass[2]=1;
                break;
            case Points:
                hasline[0]=1;
                hasclass[3]=1;
                break;
            case Elements:
                hasclass[4]=1;
                hasline[1]=1;
                break;
            case ReducedElements:
                hasclass[5]=1;
                hasline[1]=1;
                break;
            case FaceElements:
                hasclass[6]=1;
                hasline[2]=1;
                break;
            case ReducedFaceElements:
                hasclass[7]=1;
                hasline[2]=1;
                break;
            default:
                return false;
        }
    }
    int totlines=hasline[0]+hasline[1]+hasline[2]+hasline[3];

    // fail if we have more than one leaf group
    // = there are at least two branches we can't interpolate between
    if (totlines>1)
        return false;
    else if (totlines==1) {
        // we have points
        if (hasline[0]==1) {
            resultcode=Points;
        } else if (hasline[1]==1) {
            if (hasclass[5]==1) {
                resultcode=ReducedElements;
            } else {
                resultcode=Elements;
            }
        } else if (hasline[2]==1) {
            if (hasclass[7]==1) {
                resultcode=ReducedFaceElements;
            } else {
                resultcode=FaceElements;
            }
        } else { // so we must be in line3
            throw RipleyException("Programmer Error - choosing between contact elements - we should never get here");
        }
    } else { // totlines==0
        if (hasclass[2]==1) {
            // something from class 2
            resultcode=(hasrednodes ? ReducedNodes : ReducedDegreesOfFreedom);
        } else { // something from class 1
            resultcode=(hasnodes ? Nodes : DegreesOfFreedom);
        }
    }
    return true;
}

bool RipleyDomain::probeInterpolationOnDomain(int functionSpaceType_source,
                                             int functionSpaceType_target) const
{
    if (functionSpaceType_target != Nodes &&
            functionSpaceType_target != ReducedNodes &&
            functionSpaceType_target != ReducedDegreesOfFreedom &&
            functionSpaceType_target != DegreesOfFreedom &&
            functionSpaceType_target != Elements &&
            functionSpaceType_target != ReducedElements &&
            functionSpaceType_target != FaceElements &&
            functionSpaceType_target != ReducedFaceElements &&
            functionSpaceType_target != Points) {
        stringstream temp;
        temp << "Interpolation On Domain: Ripley does not know anything about function space type " << functionSpaceType_target;
        throw RipleyException(temp.str());
    }

    switch (functionSpaceType_source) {
        case Nodes:
        case DegreesOfFreedom:
            return true;
        case ReducedNodes:
        case ReducedDegreesOfFreedom:
            return (functionSpaceType_target != Nodes &&
                    functionSpaceType_target != DegreesOfFreedom);
        case Elements:
            return (functionSpaceType_target==Elements ||
                    functionSpaceType_target==ReducedElements);
        case ReducedElements:
            return (functionSpaceType_target==ReducedElements);
        case FaceElements:
            return (functionSpaceType_target==FaceElements ||
                    functionSpaceType_target==ReducedFaceElements);
        case ReducedFaceElements:
            return (functionSpaceType_target==ReducedFaceElements);
        case Points:
            return (functionSpaceType_target==Points);

        default:
            stringstream temp;
            temp << "Interpolation On Domain: Ripley does not know anything about function space type " << functionSpaceType_source;
            throw RipleyException(temp.str());
    }
    return false;
}

bool RipleyDomain::probeInterpolationACross(int functionSpaceType_source,
        const AbstractDomain& targetDomain, int functionSpaceType_target) const
{
    return false;
}

bool RipleyDomain::operator==(const AbstractDomain& other) const
{
    const RipleyDomain* temp=dynamic_cast<const RipleyDomain*>(&other);
    if (temp != NULL) {
        return (m_name==temp->m_name && m_nodes==temp->m_nodes &&
                m_elements==temp->m_elements &&
                m_faceElements==temp->m_faceElements &&
                m_points==temp->m_points);
    }
    return false;
}

bool RipleyDomain::operator!=(const AbstractDomain& other) const
{
    return !(operator==(other));
}

int RipleyDomain::getSystemMatrixTypeId(const int solver,
        const int preconditioner, const int package, const bool symmetry) const
{
    int out=Paso_SystemMatrix_getSystemMatrixTypeId(
            SystemMatrixAdapter::mapOptionToPaso(solver),
            SystemMatrixAdapter::mapOptionToPaso(preconditioner),
            SystemMatrixAdapter::mapOptionToPaso(package),
            symmetry?1:0, m_mpiInfo);
    checkPasoError();
    return out;
}

int RipleyDomain::getTransportTypeId(const int solver, const int preconditioner,
        const int package, const bool symmetry) const
{
    int out=Paso_TransportProblem_getTypeId(
            SystemMatrixAdapter::mapOptionToPaso(solver),
            SystemMatrixAdapter::mapOptionToPaso(preconditioner),
            SystemMatrixAdapter::mapOptionToPaso(package),
            symmetry?1:0, m_mpiInfo);
    checkPasoError();
    return out;
}

escript::Data RipleyDomain::getX() const
{
    return continuousFunction(*this).getX();
}

escript::Data RipleyDomain::getNormal() const
{
    return functionOnBoundary(*this).getNormal();
}

escript::Data RipleyDomain::getSize() const
{
    return escript::function(*this).getSize();
}

const int* RipleyDomain::borrowSampleReferenceIDs(int functionSpaceType) const
{
    switch (functionSpaceType) {
        case Nodes:
            return &m_nodes->getIdVector()[0];
        case ReducedNodes:
            return &m_nodes->getReducedNodesId()[0];
        case Elements:
        case ReducedElements:
            return &m_elements->getIdVector()[0];
        case FaceElements:
        case ReducedFaceElements:
            return &m_faceElements->getIdVector()[0];
        case Points:
            return &m_points->getIdVector()[0];
        case DegreesOfFreedom:
            return &m_nodes->getDegreesOfFreedomId()[0];
        case ReducedDegreesOfFreedom:
            return &m_nodes->getReducedDegreesOfFreedomId()[0];
        default:
            stringstream msg;
            msg << "borrowSampleReferenceIDs: Invalid function space type " << functionSpaceType << " for domain: " << getDescription();
            throw RipleyException(msg.str());
    }
}

int RipleyDomain::getTagFromSampleNo(int functionSpaceType, int sampleNo) const
{
    throw RipleyException("not implemented");
}


void RipleyDomain::setTags(const int functionSpaceType, const int newTag, const escript::Data& mask) const
{
    switch (functionSpaceType) {
        case Nodes:
            m_nodes->setTags(newTag, mask);
            break;
        case ReducedNodes:
            throw RipleyException("ReducedNodes does not support tags");
        case DegreesOfFreedom:
            throw RipleyException("DegreesOfFreedom does not support tags");
        case ReducedDegreesOfFreedom:
            throw RipleyException("ReducedDegreesOfFreedom does not support tags");
        case Elements:
            m_elements->setTags(newTag, mask);
            break;
        case ReducedElements:
            m_elements->setTags(newTag, mask);
            break;
        case FaceElements:
            m_faceElements->setTags(newTag, mask);
            break;
        case ReducedFaceElements:
            m_faceElements->setTags(newTag, mask);
            break;
        case Points:
            m_points->setTags(newTag, mask);
            break;
        default:
            stringstream msg;
            msg << "Ripley does not know anything about function space type " << functionSpaceType;
            throw RipleyException(msg.str());
    }
}

void RipleyDomain::setTagMap(const string& name, int tag)
{
    m_tagMap[name] = tag;
}

int RipleyDomain::getTag(const string& name) const
{
    return m_tagMap.find(name)->second;
}

bool RipleyDomain::isValidTagName(const string& name) const
{
    return (m_tagMap.find(name)!=m_tagMap.end());
}

string RipleyDomain::showTagNames() const
{
    stringstream ret;
    TagMap::const_iterator it;
    for (it=m_tagMap.begin(); it!=m_tagMap.end(); it++) {
        if (it!=m_tagMap.begin()) ret << ", ";
        ret << it->first;
    }
    return ret.str();
}

int RipleyDomain::getNumberOfTagsInUse(int functionSpaceCode) const
{
    dim_t numTags=0;
    switch (functionSpaceCode) {
        case Nodes:
            numTags=m_nodes->getNumberOfTagsInUse();
            break;
        case ReducedNodes:
            throw RipleyException("ReducedNodes does not support tags");
        case DegreesOfFreedom:
            throw RipleyException("DegreesOfFreedom does not support tags");
        case ReducedDegreesOfFreedom:
            throw RipleyException("ReducedDegreesOfFreedom does not support tags");
        case Elements:
        case ReducedElements:
            numTags=m_elements->getNumberOfTagsInUse();
            break;
        case FaceElements:
        case ReducedFaceElements:
            numTags=m_faceElements->getNumberOfTagsInUse();
            break;
        case Points:
            numTags=m_points->getNumberOfTagsInUse();
            break;
        default:
            stringstream msg;
            msg << "Ripley does not know anything about function space type " << functionSpaceCode;
            throw RipleyException(msg.str());
    }
    return numTags;
}

const int* RipleyDomain::borrowListOfTagsInUse(int functionSpaceCode) const
{
    switch (functionSpaceCode) {
        case Nodes:
            return &m_nodes->getTagsInUse()[0];
        case ReducedNodes:
            throw RipleyException("ReducedNodes does not support tags");
        case DegreesOfFreedom:
            throw RipleyException("DegreesOfFreedom does not support tags");
        case ReducedDegreesOfFreedom:
            throw RipleyException("ReducedDegreesOfFreedom does not support tags");
        case Elements:
        case ReducedElements:
            return &m_elements->getTagsInUse()[0];
        case FaceElements:
        case ReducedFaceElements:
            return &m_faceElements->getTagsInUse()[0];
        case Points:
            return &m_points->getTagsInUse()[0];
        default:
            stringstream msg;
            msg << "Ripley does not know anything about function space type " << functionSpaceCode;
            throw RipleyException(msg.str());
    }
}

bool RipleyDomain::canTag(int functionSpaceCode) const
{
    switch(functionSpaceCode) {
        case Nodes:
        case Elements:
        case ReducedElements:
        case FaceElements:
        case ReducedFaceElements:
        case Points:
            return true;
    }
    return false;
}

escript::AbstractDomain::StatusType RipleyDomain::getStatus() const
{
    return m_nodes->getStatus();
}

void RipleyDomain::prepare(bool optimize)
{
    resolveNodeIds();

    // first step is to distribute the elements according to a global
    // distribution of DOF
    // - create dense labeling for the DOFs
    dim_t newGlobalNumDOFs=m_nodes->createDenseDOFLabeling();

    // - create a distribution of the global DOFs and determine
    //   the MPI_rank controlling the DOFs on this processor
    IndexVector distribution(m_mpiInfo->size+1);
    Esys_MPIInfo_setDistribution(m_mpiInfo, 0, newGlobalNumDOFs-1, &distribution[0]);
    checkPasoError();

    // Now the mesh is redistributed according to the mpiRankOfDOF vector.
    // This will redistribute the Nodes and Elements including overlap and will
    // create an element coloring but will not create any mappings (see later
    // in this function).
    distributeByRankOfDOF(distribution);

    // At this stage we are able to start an optimization of the DOF
    // distribution using ParMetis. On return distribution is altered and new
    // DOF IDs have been assigned.
    if (optimize) {
        if (m_mpiInfo->size > 1) {
            optimizeDOFDistribution(distribution);
            distributeByRankOfDOF(distribution);
        }

        // optimize the local labeling of the degrees of freedom
        optimizeDOFLabeling(distribution);
    }

    // Rearrange elements in order to try and bring them closer to memory
    // locations of the nodes (distributed shared memory).
    optimizeElementOrdering();

/* useful DEBUG:
    {
        printf("prepare: global DOF : %d\n",newGlobalNumDOFs);
        IndexPair range = m_nodes->getGlobalIdRange();
        printf("prepare: global node id range = %d :%d\n",range.first,range.second);
        range = m_nodes->getIdRange();
        printf("prepare: local node id range = %d :%d\n",range.first,range.second);
    }
*/

    // create the global indices
    IndexVector maskReducedNodes(m_nodes->getNumNodes(), -1);
    markNodes(maskReducedNodes, 0);
    IndexVector indexReducedNodes = packMask(maskReducedNodes);

    IndexVector nodeDistribution(m_mpiInfo->size+1);
    m_nodes->createDenseNodeLabeling(nodeDistribution, distribution);
    m_nodes->createDenseReducedDOFLabeling(maskReducedNodes);
    m_nodes->createDenseReducedNodeLabeling(maskReducedNodes);

    // create the missing mappings
    m_nodes->createNodeFileMappings(indexReducedNodes, distribution,
                                    nodeDistribution);

    updateTagsInUse();
}

void RipleyDomain::optimizeElementOrdering()
{
    m_elements->optimizeOrdering();
    m_faceElements->optimizeOrdering();
    m_points->optimizeOrdering();
}

void RipleyDomain::updateTagsInUse()
{
    m_nodes->updateTagsInUse();
    m_elements->updateTagsInUse();
    m_faceElements->updateTagsInUse();
    m_points->updateTagsInUse();
}

void RipleyDomain::createColoring(const IndexVector &node_localDOF_map)
{
    m_elements->createColoring(node_localDOF_map);
    m_faceElements->createColoring(node_localDOF_map);
    m_points->createColoring(node_localDOF_map);
}

void RipleyDomain::distributeByRankOfDOF(const IndexVector &dofDistribution)
{
    RankVector mpiRankOfDOF = m_nodes->getOwnerOfDOFs(dofDistribution);

    // First the elements are redistributed according to mpiRankOfDOF.
    // At the input the Node tables refer to the local labeling of
    // the nodes while at the output they refer to the global labeling
    // which is rectified in the next step
    m_elements->distributeByRankOfDOF(mpiRankOfDOF, m_nodes->getIdVector());
    m_faceElements->distributeByRankOfDOF(mpiRankOfDOF, m_nodes->getIdVector());
    m_points->distributeByRankOfDOF(mpiRankOfDOF, m_nodes->getIdVector());

    // resolve the node IDs
    resolveNodeIds();

    // create a local labeling of the DOFs
    IndexPair dofRange = m_nodes->getDOFRange();
    dim_t len = dofRange.second - dofRange.first + 1;

    // local mask for used nodes
    IndexVector tmp_node_localDOF_mask(len, -1);
    IndexVector tmp_node_localDOF_map(m_nodes->getNumNodes(), -1);

#pragma omp parallel for schedule(static)
    for (dim_t n = 0; n < m_nodes->getNumNodes(); n++) {
#ifdef BOUNDS_CHECK
        if ((m_nodes->getGlobalDegreesOfFreedom()[n] - dofRange.first) >= len
                || (m_nodes->getGlobalDegreesOfFreedom()[n] - dofRange.first) < 0) {
            printf("BOUNDS_CHECK %s %d\n", __FILE__, __LINE__);
            exit(1);
        }
#endif
        tmp_node_localDOF_mask[m_nodes->getGlobalDegreesOfFreedom()[n] - dofRange.first] = n;
    }

    dim_t numDOFs = 0;
    for (dim_t n = 0; n < len; n++) {
        if (tmp_node_localDOF_mask[n] >= 0) {
            tmp_node_localDOF_mask[n] = numDOFs;
            numDOFs++;
        }
    }
#pragma omp parallel for
    for (dim_t n = 0; n < m_nodes->getNumNodes(); n++) {
        tmp_node_localDOF_map[n] = tmp_node_localDOF_mask[m_nodes->getGlobalDegreesOfFreedom()[n] - dofRange.first];
    }
    // create element coloring
    createColoring(tmp_node_localDOF_map);
}

/**************************************************************
   Check whether there is any node which has no vertex. In case 
   such a node exists, we don't use ParMetis since it requires
   that every node has at least 1 vertex (at line 129 of file
   "xyzpart.c" in parmetis 3.1.1, variable "nvtxs" would be 0 if 
   any node has no vertex).
 **************************************************************/
#ifdef USE_PARMETIS
static int Check_Inputs_For_Parmetis(dim_t mpiSize, dim_t rank,
                        const vector<dim_t> &distribution, MPI_Comm *comm)
{
    dim_t i, len;
    int ret_val = 1;

    if (rank == 0) {
        for (i = 0; i < mpiSize; i++) {
            len = distribution[i + 1] - distribution[i];
            if (len == 0) {
                ret_val = 0;
                break;
            }
        }
    }
    MPI_Bcast(&ret_val, 1, MPI_INTEGER, 0, *comm);
    if (ret_val == 0)
        printf("INFO: Not using ParMetis since some nodes have no vertex!\n");
    return ret_val;
}
#endif

void RipleyDomain::optimizeDOFDistribution(RankVector &distribution)
{
    const Esys_MPI_rank myRank = m_mpiInfo->rank;
    dim_t mpiSize = m_mpiInfo->size;
    dim_t dim = getDim();

    // first step is to distribute the elements according to a global X of DOF
    const index_t myFirstVertex = distribution[myRank];
    const index_t myLastVertex = distribution[myRank + 1];
    const dim_t myNumVertices = myLastVertex - myFirstVertex;
    const dim_t globalNumVertices = distribution[mpiSize];
    dim_t len = 0;
    for (dim_t p = 0; p < mpiSize; ++p)
        len = MAX(len, distribution[p+1] - distribution[p]);
    index_t *partition = TMPMEMALLOC(len, index_t);
    float *xyz = TMPMEMALLOC(myNumVertices * dim, float);
    if (!(Esys_checkPtr(partition) || Esys_checkPtr(xyz))) {
        // set the coordinates
        // it is assumed that at least one node on this processor provides a
        // coordinate
#pragma omp parallel for
        for (dim_t i = 0; i < m_nodes->getNumNodes(); ++i) {
            dim_t k = m_nodes->getGlobalDegreesOfFreedom()[i] - myFirstVertex;
            if ((k >= 0) && (k < myNumVertices)) {
                for (dim_t j = 0; j < dim; ++j)
                    xyz[k * dim + j] = (float)(m_nodes->getCoordinates()[INDEX2(j, i, dim)]);
            }
        }

        IndexMatrix indexList(myNumVertices);
        // ksteube CSR of DOF IDs
        // create the adjacency structure xadj and adjncy
        // insert contributions from element matrices into columns of indexList
        m_elements->insertIntoIndexMatrixNoMainDiagonal(indexList,
                m_nodes->getGlobalDegreesOfFreedom(),
                m_nodes->getGlobalDegreesOfFreedom(),
                myFirstVertex, myLastVertex);
        m_faceElements->insertIntoIndexMatrixNoMainDiagonal(indexList,
                m_nodes->getGlobalDegreesOfFreedom(),
                m_nodes->getGlobalDegreesOfFreedom(),
                myFirstVertex, myLastVertex);
        m_points->insertIntoIndexMatrixNoMainDiagonal(indexList,
                m_nodes->getGlobalDegreesOfFreedom(),
                m_nodes->getGlobalDegreesOfFreedom(),
                myFirstVertex, myLastVertex);

        // create the local matrix pattern
        Paso_Pattern *pattern = createPatternFromIndexMatrix(
                indexList, 0, myNumVertices, 0, globalNumVertices, 0);

#ifdef USE_PARMETIS
        if (mpiSize > 1 && Check_Inputs_For_Parmetis(mpiSize, myRank, distribution, &(m_mpiInfo->comm)) > 0) {
            int wgtflag = 0;
            int numflag = 0; /* pattern->ptr is C style: starting from 0 instead of 1 */
            int ncon = 1;
            int edgecut;
            int options[2];
            float *tpwgts = TMPMEMALLOC(ncon * mpiSize, float);
            float *ubvec = TMPMEMALLOC(ncon, float);
            for (dim_t i = 0; i < ncon * mpiSize; i++)
                tpwgts[i] = 1.0 / (float)mpiSize;
            for (dim_t i = 0; i < ncon; i++)
                ubvec[i] = 1.05;
            options[0] = 3;
            options[1] = 15;
            ParMETIS_V3_PartGeomKway(&distribution[0], pattern->ptr,
                    pattern->index, NULL, NULL, &wgtflag, &numflag,
                    &dim, xyz, &ncon, &mpiSize, tpwgts, ubvec, options,
                    &edgecut, partition, // new CPU ownership of elements
                    &m_mpiInfo->comm);
            TMPMEMFREE(ubvec);
            TMPMEMFREE(tpwgts);
        } else {
            for (dim_t i = 0; i < myNumVertices; ++i)
                partition[i] = 0; // CPU 0 owns it
        }
#else
        for (dim_t i = 0; i < myNumVertices; ++i)
            partition[i] = myRank; // CPU myRank owns it
#endif

        Paso_Pattern_free(pattern);

        // create a new distribution and labeling of the DOF
        RankVector newDistribution(mpiSize+1, 0);
#pragma omp parallel
        {
            RankVector loc_partition_count(mpiSize, 0);
#pragma omp for
            for (dim_t i = 0; i < myNumVertices; ++i)
                loc_partition_count[partition[i]]++;
#pragma omp critical
            {
                for (dim_t i = 0; i < mpiSize; ++i)
                    newDistribution[i] += loc_partition_count[i];
            }
        }

        // recvbuf will be the concatenation of each CPU's contribution
        // to newDistribution
        dim_t *recvbuf = TMPMEMALLOC(mpiSize * mpiSize, dim_t);

#ifdef ESYS_MPI
        MPI_Allgather(&newDistribution[0], mpiSize, MPI_INT, recvbuf, mpiSize, MPI_INT, m_mpiInfo->comm);
#else
        for (dim_t i = 0; i < mpiSize; ++i)
            recvbuf[i] = newDistribution[i];
#endif
        newDistribution[0] = 0;
        IndexVector newGlobalDOFid(len);
        for (Esys_MPI_rank rank = 0; rank < mpiSize; rank++) {
            int c = 0;
            for (dim_t i = 0; i < myRank; ++i)
                c += recvbuf[rank + mpiSize * i];
            for (dim_t i = 0; i < myNumVertices; ++i) {
                if (rank == partition[i]) {
                    newGlobalDOFid[i] = newDistribution[rank] + c;
                    c++;
                }
            }
            for (dim_t i = myRank + 1; i < mpiSize; ++i)
                c += recvbuf[rank + mpiSize * i];
            newDistribution[rank + 1] = newDistribution[rank] + c;
        }
        TMPMEMFREE(recvbuf);

        // now the overlap needs to be created by sending the partition
        // around
        m_nodes->resetGlobalDegreesOfFreedom(newGlobalDOFid, distribution);
        for (dim_t i = 0; i < mpiSize+1; ++i)
            distribution[i] = newDistribution[i];
    }
    TMPMEMFREE(partition);
    TMPMEMFREE(xyz);
}

void RipleyDomain::optimizeDOFLabeling(const IndexVector &distribution)
{
    Esys_MPI_rank myRank = m_mpiInfo->rank;
    index_t myFirstVertex = distribution[myRank];
    index_t myLastVertex = distribution[myRank + 1];
    dim_t myNumVertices = myLastVertex - myFirstVertex;
    dim_t len = 0;
    for (dim_t p = 0; p < m_mpiInfo->size; ++p)
        len = MAX(len, distribution[p + 1] - distribution[p]);

    IndexMatrix indexList(myNumVertices);

    // create the adjacency structure xadj and adjncy
    // insert contributions from element matrices into columns of indexList
    m_elements->insertIntoIndexMatrixNoMainDiagonal(indexList,
            m_nodes->getGlobalDegreesOfFreedom(),
            m_nodes->getGlobalDegreesOfFreedom(),
            myFirstVertex, myLastVertex);
    m_faceElements->insertIntoIndexMatrixNoMainDiagonal(indexList,
            m_nodes->getGlobalDegreesOfFreedom(),
            m_nodes->getGlobalDegreesOfFreedom(),
            myFirstVertex, myLastVertex);
    m_points->insertIntoIndexMatrixNoMainDiagonal(indexList,
            m_nodes->getGlobalDegreesOfFreedom(),
            m_nodes->getGlobalDegreesOfFreedom(),
            myFirstVertex, myLastVertex);

    // create the local matrix pattern
    Paso_Pattern *pattern = createPatternFromIndexMatrix(indexList, 0,
            myNumVertices, myFirstVertex, myLastVertex, -myFirstVertex);

    IndexVector newGlobalDOFid(len);
    Paso_Pattern_reduceBandwidth(pattern, &newGlobalDOFid[0]);
    Paso_Pattern_free(pattern);

    Esys_MPIInfo_noError(m_mpiInfo);
    checkPasoError();

    /* shift new labeling to create a global id */
#pragma omp parallel for
    for (dim_t i = 0; i < myNumVertices; ++i)
        newGlobalDOFid[i] += myFirstVertex;

    /* distribute new labeling to other processors */
    m_nodes->resetGlobalDegreesOfFreedom(newGlobalDOFid, distribution);
#if 0
    for (i = 0; i < m_nodes->getNumNodes(); ++i)
        printf("%d ", m_nodes->getGlobalDegreesOfFreedom()[i]);
    printf("\n");
#endif
}

void RipleyDomain::resolveNodeIds()
{
    // find the minimum and maximum node ID used by elements
    IndexPair minMaxId = m_elements->getNodeRange();
    index_t minId = minMaxId.first;
    index_t maxId = minMaxId.second;
    minMaxId = m_faceElements->getNodeRange();
    minId = MIN(minId, minMaxId.first);
    maxId = MAX(maxId, minMaxId.second);
    minMaxId = m_points->getNodeRange();
    minId = MIN(minId, minMaxId.first);
    maxId = MAX(maxId, minMaxId.second);

#ifdef Ripley_TRACE
    index_t globalMinId, globalMaxId;
#ifdef ESYS_MPI
    index_t idRange[2], globalIdRange[2];
    idRange[0] = -minId;
    idRange[1] = maxId;
    MPI_Allreduce(idRange, globalIdRange, 2, MPI_INT, MPI_MAX, m_mpiInfo->comm);
    globalMinId = -globalIdRange[0];
    globalMaxId = globalIdRange[1];
#else
    globalMinId = minId;
    globalMaxId = maxId;
#endif // ESYS_MPI
    cout << "Node id range used by elements is " << globalMinId << ":"
        << globalMaxId << endl;
#endif // Ripley_TRACE

    // this is only true if we have no elements at all
    if (minId > maxId) {
        maxId = -1;
        minId = 0;
    }

    // allocate mappings for new local node labeling to global node labeling
    // (newLocalToGlobalNodeLabels) and global node labeling to the new local
    // node labeling. globalToNewLocalNodeLabels[i-minId] is the new local id
    // of global node i
    dim_t len = (maxId >= minId) ? maxId - minId + 1 : 0;
    IndexVector globalToNewLocalNodeLabels(len, -1);

    // mark the nodes referred by elements in globalToNewLocalNodeLabels
    // which is currently used as a mask
    markNodes(globalToNewLocalNodeLabels, minId);

    // create a local labeling newLocalToGlobalNodeLabels of the local
    // nodes by packing the mask globalToNewLocalNodeLabels
    IndexVector newLocalToGlobalNodeLabels = packMask(globalToNewLocalNodeLabels);
    const dim_t newNumNodes = newLocalToGlobalNodeLabels.size();

    // invert the new labeling and shift the index newLocalToGlobalNodeLabels
    // to global node IDs
#pragma omp parallel for schedule(static)
    for (dim_t n = 0; n < newNumNodes; n++) {
#ifdef BOUNDS_CHECK
        if (n >= len || n < 0) {
            printf("BOUNDS_CHECK %s %d n=%d\n", __FILE__, __LINE__, n);
            exit(1);
        }
        if (newLocalToGlobalNodeLabels[n] >= len || newLocalToGlobalNodeLabels[n] < 0) {
            printf("BOUNDS_CHECK %s %d n=%d\n", __FILE__, __LINE__, n);
            exit(1);
        }
#endif
        globalToNewLocalNodeLabels[newLocalToGlobalNodeLabels[n]] = n;
        newLocalToGlobalNodeLabels[n] += minId;
    }

    // create a new node file
    m_nodes = m_nodes->gatherGlobal(newLocalToGlobalNodeLabels);

    // relabel nodes of the elements
    relabelElementNodes(globalToNewLocalNodeLabels, minId);
}

void RipleyDomain::relabelElementNodes(const IndexVector &newNode, index_t offset)
{
    m_elements->relabelNodes(newNode, offset);
    m_faceElements->relabelNodes(newNode, offset);
    m_points->relabelNodes(newNode, offset);
}

void RipleyDomain::markNodes(IndexVector &mask, index_t offset)
{
    m_elements->markNodes(mask, offset);
    m_faceElements->markNodes(mask, offset);
    m_points->markNodes(mask, offset);
}

void RipleyDomain::createMappings(const IndexVector &dofDistribution, const IndexVector &nodeDistribution)
{
    IndexVector maskReducedNodes(m_nodes->getNumNodes(), -1);
    markNodes(maskReducedNodes, 0);

    IndexVector indexReducedNodes = packMask(maskReducedNodes);
    m_nodes->createNodeFileMappings(indexReducedNodes, dofDistribution,
                                    nodeDistribution);
}

Paso_SystemMatrixPattern *RipleyDomain::makePattern(bool reduce_row_order, bool reduce_col_order) const
{
    Paso_Connector *colConnector, *rowConnector;
    NodeMapping_ptr colMap, rowMap;
    Paso_Distribution *colDistribution = NULL, *rowDistribution = NULL;

    if (reduce_col_order) {
        colMap = m_nodes->getReducedDegreesOfFreedomMapping();
        colDistribution = m_nodes->getReducedDegreesOfFreedomDistribution();
        colConnector = m_nodes->getReducedDegreesOfFreedomConnector();
    } else {
        colMap = m_nodes->getDegreesOfFreedomMapping();
        colDistribution = m_nodes->getDegreesOfFreedomDistribution();
        colConnector = m_nodes->getDegreesOfFreedomConnector();
    }

    if (reduce_row_order) {
        rowMap = m_nodes->getReducedDegreesOfFreedomMapping();
        rowDistribution = m_nodes->getReducedDegreesOfFreedomDistribution();
        rowConnector = m_nodes->getReducedDegreesOfFreedomConnector();
    } else {
        rowMap = m_nodes->getDegreesOfFreedomMapping();
        rowDistribution = m_nodes->getDegreesOfFreedomDistribution();
        rowConnector = m_nodes->getDegreesOfFreedomConnector();
    }

    //  insert contributions from element matrices into columns in indexList
    IndexMatrix indexList(rowMap->numTargets);
    m_elements->insertIntoIndexMatrix(indexList, rowMap->target, colMap->target);
    m_faceElements->insertIntoIndexMatrix(indexList, rowMap->target, colMap->target);
    m_points->insertIntoIndexMatrix(indexList, rowMap->target, colMap->target);

    // create patterns
    Paso_Pattern *mainPattern = createPatternFromIndexMatrix(indexList, 0,
            Paso_Distribution_getMyNumComponents(rowDistribution),
            0, Paso_Distribution_getMyNumComponents(colDistribution), 0);
    Paso_Pattern *colCouplePattern = createPatternFromIndexMatrix(indexList, 0,
            Paso_Distribution_getMyNumComponents(rowDistribution),
            Paso_Distribution_getMyNumComponents(colDistribution),
            colMap->numTargets,
            -Paso_Distribution_getMyNumComponents(colDistribution));
    Paso_Pattern *rowCouplePattern = createPatternFromIndexMatrix(indexList,
            Paso_Distribution_getMyNumComponents(rowDistribution),
            rowMap->numTargets, 0,
            Paso_Distribution_getMyNumComponents(colDistribution), 0);

    // if everything is in order we can create the return value
    Paso_SystemMatrixPattern *out = Paso_SystemMatrixPattern_alloc(
            MATRIX_FORMAT_DEFAULT, rowDistribution, colDistribution,
            mainPattern, colCouplePattern, rowCouplePattern,
            colConnector, rowConnector);

    // clean up
    Paso_Pattern_free(mainPattern);
    Paso_Pattern_free(colCouplePattern);
    Paso_Pattern_free(rowCouplePattern);
    Esys_MPIInfo_noError(m_mpiInfo);
    return out;
}

Paso_SystemMatrixPattern *RipleyDomain::getPattern(bool reduce_row_order, bool reduce_col_order) const
{
    return Paso_SystemMatrixPattern_getReference(makePattern(reduce_row_order, reduce_col_order));
    /*
    Paso_SystemMatrixPattern *out = NULL;

    // make sure that the requested pattern is available
    if (reduce_row_order) {
        if (reduce_col_order) {
            if (m_reducedReducedPattern == NULL)
                m_reducedReducedPattern = makePattern(reduce_row_order, reduce_col_order);
        } else {
            if (m_reducedFullPattern == NULL)
                m_reducedFullPattern = makePattern(reduce_row_order, reduce_col_order);
        }
    } else {
        if (reduce_col_order) {
            if (m_fullReducedPattern == NULL)
                m_fullReducedPattern = makePattern(reduce_row_order, reduce_col_order);
        } else {
            if (m_fullFullPattern == NULL)
                m_fullFullPattern = makePattern(reduce_row_order, reduce_col_order);
        }
    }
    if (Ripley_noError()) {
        if (reduce_row_order) {
            if (reduce_col_order) {
                out = Paso_SystemMatrixPattern_getReference(m_reducedReducedPattern);
            } else {
                out = Paso_SystemMatrixPattern_getReference(m_reducedFullPattern);
            }
        } else {
            if (reduce_col_order) {
                out = Paso_SystemMatrixPattern_getReference(m_fullReducedPattern);
            } else {
                out = Paso_SystemMatrixPattern_getReference(m_fullFullPattern);
            }
        }
    }
    return out;
    */
}

} // end of namespace ripley


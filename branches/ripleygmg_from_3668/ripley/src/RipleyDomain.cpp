
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
#include <escript/DataFactory.h>
#include <escript/FunctionSpaceFactory.h>
#include <pasowrap/SystemMatrixAdapter.h>
#include <pasowrap/TransportProblemAdapter.h>

#include <iomanip>

using namespace std;
using paso::SystemMatrixAdapter;
using paso::TransportProblemAdapter;

namespace ripley {

escript::Domain_ptr RipleyDomain::loadMesh(const string& filename)
{
    throw RipleyException("loadMesh() not implemented");
}

escript::Domain_ptr RipleyDomain::readMesh(const string& filename)
{
    throw RipleyException("readMesh() not implemented");
}

RipleyDomain::RipleyDomain(dim_t dim) :
    m_numDim(dim),
    m_status(0)
{
    m_mpiInfo = Esys_MPIInfo_alloc(MPI_COMM_WORLD);
}

RipleyDomain::~RipleyDomain()
{
    Esys_MPIInfo_free(m_mpiInfo);
}

bool RipleyDomain::isValidFunctionSpaceType(int fsType) const
{
    switch (fsType) {
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
        case Nodes:
        case ReducedNodes:
        case Elements:
        case ReducedElements:
        case FaceElements:
        case ReducedFaceElements:
        case Points:
            return true;
        default:
            break;
    }
    return false;
}

string RipleyDomain::functionSpaceTypeAsString(int fsType) const
{
    switch (fsType) {
        case DegreesOfFreedom: return "Ripley_DegreesOfFreedom";
        case ReducedDegreesOfFreedom: return "Ripley_ReducedDegreesOfFreedom";
        case Nodes: return "Ripley_Nodes";
        case ReducedNodes: return "Ripley_Reduced_Nodes";
        case Elements: return "Ripley_Elements";
        case ReducedElements: return "Ripley_Reduced_Elements";
        case FaceElements: return "Ripley_Face_Elements";
        case ReducedFaceElements: return "Ripley_Reduced_Face_Elements";
        case Points: return "Ripley_Points";
        default:
            break;
    }
    return "Invalid function space type code";
}

pair<int,int> RipleyDomain::getDataShape(int fsType) const
{
    const int ptsPerSample = (m_numDim==2 ? 4 : 8);
    switch (fsType) {
        case Nodes:
        case ReducedNodes:
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
            return pair<int,int>(1, getNumNodes());
        case Elements:
            return pair<int,int>(ptsPerSample, getNumElements());
        case FaceElements:
            return pair<int,int>(ptsPerSample/2, getNumFaceElements());
        case ReducedElements:
            return pair<int,int>(1, getNumElements());
        case ReducedFaceElements:
            return pair<int,int>(1, getNumFaceElements());
            /*
        case Points:
            return pair<int,int>(1, getNumPoints());
            */
        default:
            break;
    }

    stringstream msg;
    msg << "getDataShape(): Unsupported function space type "
        << functionSpaceTypeAsString(fsType) << " for " << getDescription();
    throw RipleyException(msg.str());
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

bool RipleyDomain::commonFunctionSpace(const vector<int>& fs, int& resultcode) const
{
   /*
    The idea is to use equivalence classes (i.e. types which can be
    interpolated back and forth):
    class 1: DOF <-> Nodes
    class 2: ReducedDOF <-> ReducedNodes
    class 3: Points
    class 4: Elements
    class 5: ReducedElements
    class 6: FaceElements
    class 7: ReducedFaceElements
    class 8: ContactElementZero <-> ContactElementOne
    class 9: ReducedContactElementZero <-> ReducedContactElementOne

    There is also a set of lines. Interpolation is possible down a line but not
    between lines.
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
        if (hasline[0]==1)
            resultcode=Points;
        else if (hasline[1]==1) {
            if (hasclass[5]==1)
                resultcode=ReducedElements;
            else
                resultcode=Elements;
        } else if (hasline[2]==1) {
            if (hasclass[7]==1)
                resultcode=ReducedFaceElements;
            else
                resultcode=FaceElements;
        } else
            throw RipleyException("Internal Ripley Error!");
    } else { // totlines==0
        if (hasclass[2]==1)
            // something from class 2
            resultcode=(hasrednodes ? ReducedNodes : ReducedDegreesOfFreedom);
        else // something from class 1
            resultcode=(hasnodes ? Nodes : DegreesOfFreedom);
    }
    return true;
}

void RipleyDomain::interpolateOnDomain(escript::Data& target,
                                       const escript::Data& in) const
{
    const RipleyDomain& inDomain=dynamic_cast<const RipleyDomain&>(*(in.getFunctionSpace().getDomain()));
    const RipleyDomain& targetDomain=dynamic_cast<const RipleyDomain&>(*(target.getFunctionSpace().getDomain()));
    if (inDomain != *this)
        throw RipleyException("Illegal domain of interpolant");
    if (targetDomain != *this)
        throw RipleyException("Illegal domain of interpolation target");

    stringstream msg;
    msg << "interpolateOnDomain() not implemented for function space "
        << functionSpaceTypeAsString(in.getFunctionSpace().getTypeCode())
        << " -> "
        << functionSpaceTypeAsString(target.getFunctionSpace().getTypeCode());

    switch (in.getFunctionSpace().getTypeCode()) {
        case Nodes:
        case ReducedNodes:
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
            switch (target.getFunctionSpace().getTypeCode()) {
                case Nodes:
                case ReducedNodes:
                case DegreesOfFreedom:
                case ReducedDegreesOfFreedom:
                    // FIXME
                    copyNodalData(target, *const_cast<escript::Data*>(&in));
                    break;

                case Elements:
                    interpolateNodesOnElements(target, *const_cast<escript::Data*>(&in), false);
                    break;

                case ReducedElements:
                    interpolateNodesOnElements(target, *const_cast<escript::Data*>(&in), true);
                    break;

                case FaceElements:
                    interpolateNodesOnFaces(target, *const_cast<escript::Data*>(&in), false);
                    break;

                case ReducedFaceElements:
                    interpolateNodesOnFaces(target, *const_cast<escript::Data*>(&in), true);
                    break;

                default:
                    throw RipleyException(msg.str());
            }
            break;
        default:
            throw RipleyException(msg.str());
    }
}

escript::Data RipleyDomain::getX() const
{
    return escript::continuousFunction(*this).getX();
}

escript::Data RipleyDomain::getNormal() const
{
    return escript::functionOnBoundary(*this).getNormal();
}

escript::Data RipleyDomain::getSize() const
{
    return escript::function(*this).getSize();
}

void RipleyDomain::setToX(escript::Data& arg) const
{
    const RipleyDomain& argDomain=dynamic_cast<const RipleyDomain&>(
            *(arg.getFunctionSpace().getDomain()));
    if (argDomain != *this)
        throw RipleyException("setToX: Illegal domain of data point locations");
    if (!arg.isExpanded())
        throw RipleyException("setToX: Expanded Data object expected");

    if (arg.getFunctionSpace().getTypeCode()==Nodes) {
        assembleCoordinates(arg);
    } else {
        // interpolate the result
        escript::Data contData=escript::Vector(0., escript::continuousFunction(*this), true);
        assembleCoordinates(contData);
        interpolateOnDomain(arg, contData);
    }
}

bool RipleyDomain::isCellOriented(int fsType) const
{
    switch(fsType) {
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
            break;
    }
    stringstream msg;
    msg << "isCellOriented(): Illegal function space type " << fsType << " on " << getDescription();
    throw RipleyException(msg.str());
}

bool RipleyDomain::canTag(int fsType) const
{
    switch(fsType) {
        case Nodes:
        case Elements:
        case ReducedElements:
        case FaceElements:
        case ReducedFaceElements:
            return true;
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
        case Points:
            return false;
        default:
            break;
    }
    stringstream msg;
    msg << "canTag(): Illegal function space type " << fsType << " on " << getDescription();
    throw RipleyException(msg.str());
}

void RipleyDomain::setTags(const int fsType, const int newTag, const escript::Data& cMask) const
{
    IndexVector* target=NULL;
    dim_t num=0;

    switch(fsType) {
        case Nodes:
            num=getNumNodes();
            target=&m_nodeTags;
            break;
        case Elements:
        case ReducedElements:
            num=getNumElements();
            target=&m_elementTags;
            break;
        case FaceElements:
        case ReducedFaceElements:
            num=getNumFaceElements();
            target=&m_faceTags;
            break;
        default: {
            stringstream msg;
            msg << "setTags(): not implemented for " << functionSpaceTypeAsString(fsType);
            throw RipleyException(msg.str());
        }
    }
    if (target->size() != num) {
        target->assign(num, -1);
    }

    escript::Data& mask=*const_cast<escript::Data*>(&cMask);
#pragma omp parallel for
    for (index_t i=0; i<num; i++) {
        if (mask.getSampleDataRO(i)[0] > 0) {
            (*target)[i]=newTag;
        }
    }
    updateTagsInUse(fsType);
}

int RipleyDomain::getTagFromSampleNo(int fsType, int sampleNo) const
{
    switch(fsType) {
        case Nodes:
            if (m_nodeTags.size() > sampleNo)
                return m_nodeTags[sampleNo];
            break;
        case Elements:
        case ReducedElements:
            if (m_elementTags.size() > sampleNo)
                return m_elementTags[sampleNo];
            break;
        case FaceElements:
        case ReducedFaceElements:
            if (m_faceTags.size() > sampleNo)
                return m_faceTags[sampleNo];
            break;
        default: {
            stringstream msg;
            msg << "getTagFromSampleNo(): not implemented for " << functionSpaceTypeAsString(fsType);
            throw RipleyException(msg.str());
        }
    }
    return -1;
}

int RipleyDomain::getNumberOfTagsInUse(int fsType) const
{
    switch(fsType) {
        case Nodes:
            return m_nodeTagsInUse.size();
        case Elements:
        case ReducedElements:
            return m_elementTagsInUse.size();
        case FaceElements:
        case ReducedFaceElements:
            return m_faceTagsInUse.size();
        default: {
            stringstream msg;
            msg << "getNumberOfTagsInUse(): not implemented for " << functionSpaceTypeAsString(fsType);
            throw RipleyException(msg.str());
        }
    }
}

const int* RipleyDomain::borrowListOfTagsInUse(int fsType) const
{
    switch(fsType) {
        case Nodes:
            return &m_nodeTagsInUse[0];
        case Elements:
        case ReducedElements:
            return &m_elementTagsInUse[0];
        case FaceElements:
        case ReducedFaceElements:
            return &m_faceTagsInUse[0];
        default: {
            stringstream msg;
            msg << "borrowListOfTagsInUse(): not implemented for " << functionSpaceTypeAsString(fsType);
            throw RipleyException(msg.str());
        }
    }
}

void RipleyDomain::Print_Mesh_Info(const bool full) const
{
    cout << "Print_Mesh_Info for " << getDescription() << " running on CPU "
        << m_mpiInfo->rank << ". MPI size: " << m_mpiInfo->size << endl;
    cout << "Number of dimensions: " << m_numDim << endl;

    // write tags
    if (m_tagMap.size() > 0) {
        cout << "Tags:" << endl;
        TagMap::const_iterator it;
        for (it=m_tagMap.begin(); it!=m_tagMap.end(); it++) {
            cout << "  " << setw(5) << it->second << " "
                << it->first << endl;
        }
    }
}

int RipleyDomain::getSystemMatrixTypeId(const int solver,
        const int preconditioner, const int package, const bool symmetry) const
{
    return SystemMatrixAdapter::getSystemMatrixTypeId(solver, preconditioner,
            package, symmetry, m_mpiInfo);
}

int RipleyDomain::getTransportTypeId(const int solver, const int preconditioner,
        const int package, const bool symmetry) const
{
    return TransportProblemAdapter::getTransportTypeId(solver, preconditioner,
            package, symmetry, m_mpiInfo);
}

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
        throw RipleyException("Domain of column function space does not match the domain of matrix generator");
    // is the function space type right?
    if (row_functionspace.getTypeCode()==ReducedDegreesOfFreedom)
        reduceRowOrder=true;
    else if (row_functionspace.getTypeCode()!=DegreesOfFreedom)
        throw RipleyException("Illegal function space type for system matrix rows");
    if (column_functionspace.getTypeCode()==ReducedDegreesOfFreedom)
        reduceColOrder=true;
    else if (column_functionspace.getTypeCode()!=DegreesOfFreedom)
        throw RipleyException("Illegal function space type for system matrix columns");

    // generate matrix
    Paso_SystemMatrixPattern* pattern=getPattern(reduceRowOrder, reduceColOrder);
    Paso_SystemMatrix* matrix = Paso_SystemMatrix_alloc(type, pattern,
            row_blocksize, column_blocksize, FALSE);
    paso::checkPasoError();
    Paso_SystemMatrixPattern_free(pattern);
    escript::ASM_ptr sma(new SystemMatrixAdapter(matrix, row_blocksize,
                row_functionspace, column_blocksize, column_functionspace));
    return sma;
}

void RipleyDomain::setNewX(const escript::Data& arg)
{
    throw RipleyException("setNewX(): Operation not supported");
}

//protected
void RipleyDomain::copyNodalData(escript::Data& out, escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
#pragma omp parallel for
    for (index_t i=0; i<in.getNumSamples(); i++) {
        const double* src = in.getSampleDataRO(i);
        copy(src, src+numComp, out.getSampleDataRW(i));
    }
}

//protected
void RipleyDomain::updateTagsInUse(int fsType) const
{
    IndexVector* tagsInUse=NULL;
    const IndexVector* tags=NULL;
    switch(fsType) {
        case Nodes:
            tags=&m_nodeTags;
            tagsInUse=&m_nodeTagsInUse;
            break;
        case Elements:
        case ReducedElements:
            tags=&m_elementTags;
            tagsInUse=&m_elementTagsInUse;
            break;
        case FaceElements:
        case ReducedFaceElements:
            tags=&m_faceTags;
            tagsInUse=&m_faceTagsInUse;
            break;
        default:
            return;
    }

    // gather global unique tag values from tags into tagsInUse
    tagsInUse->clear();
    index_t lastFoundValue = INDEX_T_MIN, minFoundValue, local_minFoundValue;

    while (true) {
        // find smallest value bigger than lastFoundValue 
        minFoundValue = INDEX_T_MAX;
#pragma omp parallel private(local_minFoundValue)
        {
            local_minFoundValue = minFoundValue;
#pragma omp for schedule(static) nowait
            for (size_t i = 0; i < tags->size(); i++) {
                const index_t v = (*tags)[i];
                if ((v > lastFoundValue) && (v < local_minFoundValue))
                    local_minFoundValue = v;
            }
#pragma omp critical
            {
                if (local_minFoundValue < minFoundValue)
                    minFoundValue = local_minFoundValue;
            }
        }
#ifdef ESYS_MPI
        local_minFoundValue = minFoundValue;
        MPI_Allreduce(&local_minFoundValue, &minFoundValue, 1, MPI_INT, MPI_MIN, m_mpiInfo->comm);
#endif
        // if we found a new value add it to the tagsInUse vector

        if (minFoundValue < INDEX_T_MAX) {
            tagsInUse->push_back(minFoundValue);
            lastFoundValue = minFoundValue;
        } else
            break;
    }
}

//
// the methods that follow have to be implemented by the subclasses
//

string RipleyDomain::getDescription() const
{
    throw RipleyException("getDescription() not implemented");
}

bool RipleyDomain::operator==(const AbstractDomain& other) const
{
    throw RipleyException("operator==() not implemented");
}

void RipleyDomain::write(const string& filename) const
{
    throw RipleyException("write() not implemented");
}

void RipleyDomain::dump(const string& filename) const
{
    throw RipleyException("dump() not implemented");
}

const int* RipleyDomain::borrowSampleReferenceIDs(int fsType) const
{
    throw RipleyException("borrowSampleReferenceIDs() not implemented");
}

bool RipleyDomain::probeInterpolationOnDomain(int fsType_source,
                                             int fsType_target) const
{
    throw RipleyException("probeInterpolationOnDomain() not implemented");
}

void RipleyDomain::interpolateACross(escript::Data& target, const escript::Data& source) const
{
    throw RipleyException("interpolateACross() not implemented");
}

bool RipleyDomain::probeInterpolationACross(int fsType_source,
        const escript::AbstractDomain&, int fsType_target) const
{
    throw RipleyException("probeInterpolationACross() not implemented");
}

void RipleyDomain::setToNormal(escript::Data& normal) const
{
    throw RipleyException("setToNormal() not implemented");
}

void RipleyDomain::setToSize(escript::Data& size) const
{
    throw RipleyException("setToSize() not implemented");
}

void RipleyDomain::setToGradient(escript::Data& grad, const escript::Data& arg) const
{
    throw RipleyException("setToGradient() not implemented");
}

void RipleyDomain::setToIntegrals(vector<double>& integrals, const escript::Data& arg) const
{
    throw RipleyException("setToIntegrals() not implemented");
}

bool RipleyDomain::ownSample(int fsType, index_t id) const
{
    throw RipleyException("ownSample() not implemented");
}

void RipleyDomain::addPDEToSystem(
        escript::AbstractSystemMatrix& mat, escript::Data& rhs,
        const escript::Data& A, const escript::Data& B, const escript::Data& C,
        const escript::Data& D, const escript::Data& X, const escript::Data& Y,
        const escript::Data& d, const escript::Data& y,
        const escript::Data& d_contact, const escript::Data& y_contact,
        const escript::Data& d_dirac,const escript::Data& y_dirac) const
{
    throw RipleyException("addPDEToSystem() not implemented");
}

void RipleyDomain::addPDEToLumpedSystem(escript::Data& mat,
        const escript::Data& D, const escript::Data& d,
        const escript::Data& d_dirac, const bool useHRZ) const
{
    throw RipleyException("addPDEToLumpedSystem() not implemented");
}

void RipleyDomain::addPDEToRHS(escript::Data& rhs, const escript::Data& X,
        const escript::Data& Y, const escript::Data& y,
        const escript::Data& y_contact, const escript::Data& y_dirac) const
{
    throw RipleyException("addPDEToRHS() not implemented");
}

void RipleyDomain::addPDEToTransportProblem(
        escript::AbstractTransportProblem& tp,
        escript::Data& source, const escript::Data& M,
        const escript::Data& A, const escript::Data& B, const escript::Data& C,
        const escript::Data& D, const escript::Data& X, const escript::Data& Y,
        const escript::Data& d, const escript::Data& y,
        const escript::Data& d_contact, const escript::Data& y_contact,
        const escript::Data& d_dirac, const escript::Data& y_dirac) const
{
    throw RipleyException("addPDEToTransportProblem() not implemented");
}

escript::ATP_ptr RipleyDomain::newTransportProblem(const bool useBackwardEuler,
        const int blocksize, const escript::FunctionSpace& functionspace,
        const int type) const
{
    throw RipleyException("newTransportProblem() not implemented");
}

Paso_SystemMatrixPattern* RipleyDomain::getPattern(bool reducedRowOrder,
                                                   bool reducedColOrder) const
{
    throw RipleyException("getPattern() not implemented");
}

dim_t RipleyDomain::getNumDataPointsGlobal() const
{
    throw RipleyException("getNumDataPointsGlobal() not implemented");
}

IndexVector RipleyDomain::getNumNodesPerDim() const
{
    throw RipleyException("getNumNodesPerDim() not implemented");
}

IndexVector RipleyDomain::getNumElementsPerDim() const
{
    throw RipleyException("getNumElementsPerDim() not implemented");
}

IndexVector RipleyDomain::getNumFacesPerBoundary() const
{
    throw RipleyException("getNumFacesPerBoundary() not implemented");
}

IndexVector RipleyDomain::getNodeDistribution() const
{
    throw RipleyException("getNodeDistribution() not implemented");
}

pair<double,double> RipleyDomain::getFirstCoordAndSpacing(dim_t dim) const
{
    throw RipleyException("getFirstCoordAndSpacing() not implemented");
}

dim_t RipleyDomain::getNumFaceElements() const
{
    throw RipleyException("getNumFaceElements() not implemented");
}

dim_t RipleyDomain::getNumElements() const
{
    throw RipleyException("getNumElements() not implemented");
}

dim_t RipleyDomain::getNumNodes() const
{
    throw RipleyException("getNumNodes() not implemented");
}

void RipleyDomain::assembleCoordinates(escript::Data& arg) const
{
    throw RipleyException("assembleCoordinates() not implemented");
}

void RipleyDomain::interpolateNodesOnElements(escript::Data& out, escript::Data& in, bool reduced) const
{
    throw RipleyException("interpolateNodesOnElements() not implemented");
}

void RipleyDomain::interpolateNodesOnFaces(escript::Data& out, escript::Data& in, bool reduced) const
{
    throw RipleyException("interpolateNodesOnFaces() not implemented");
}


} // end of namespace ripley



/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include <speckley/SpeckleyDomain.h>
#include <escript/DataFactory.h>
#include <escript/FunctionSpaceFactory.h>
#include <speckley/domainhelpers.h>

#include <iomanip>

namespace bp = boost::python;

using namespace std;

namespace speckley {

void tupleListToMap(DataMap& mapping, const bp::list& list)
{
    for (int i = 0; i < len(list); i++) {
        if (!bp::extract<bp::tuple>(list[i]).check())
            throw SpeckleyException("Passed in list contains objects"
                                    " other than tuples");
        bp::tuple t = bp::extract<bp::tuple>(list[i]);
        if (len(t) != 2 || !bp::extract<string>(t[0]).check() ||
                !bp::extract<escript::Data>(t[1]).check())
            throw SpeckleyException("The passed in list must contain tuples"
                " of the form (string, escript::data)");
        mapping[bp::extract<string>(t[0])] = bp::extract<escript::Data>(t[1]);
    }
}

SpeckleyDomain::SpeckleyDomain(dim_t dim, int order, escript::SubWorld_ptr p) :
    m_numDim(dim),
    m_status(0),
    m_order(order)
{
    if (p.get() == NULL)
        m_mpiInfo = esysUtils::makeInfo(MPI_COMM_WORLD);
    else
        m_mpiInfo = p->getMPI();

    assembler_type = DEFAULT_ASSEMBLER;
}

SpeckleyDomain::~SpeckleyDomain()
{
    // cleanup of MPI is dealt with by shared_ptr
}

bool SpeckleyDomain::operator==(const AbstractDomain& other) const
{
    const SpeckleyDomain* o=dynamic_cast<const SpeckleyDomain*>(&other);
    if (o) {
        return (m_tagMap==o->m_tagMap && m_nodeTags==o->m_nodeTags
                && m_elementTags==o->m_elementTags
                && m_faceTags==o->m_faceTags);
    }
    return false;
}

bool SpeckleyDomain::isValidFunctionSpaceType(int fsType) const
{
    switch (fsType) {
        case DegreesOfFreedom:
//        case ReducedDegreesOfFreedom:
        case Nodes:
//        case ReducedNodes:
        case Elements:
//        case ReducedElements:
//        case FaceElements:
//        case ReducedFaceElements:
        case Points:
            return true;
        default:
            break;
    }
    return false;
}

string SpeckleyDomain::functionSpaceTypeAsString(int fsType) const
{
    switch (fsType) {
        case DegreesOfFreedom: return "Speckley_DegreesOfFreedom [Solution(domain)]";
        case ReducedDegreesOfFreedom: return "Speckley_ReducedDegreesOfFreedom [ReducedSolution(domain)]";
        case Nodes: return "Speckley_Nodes [ContinuousFunction(domain)]";
        case ReducedNodes: return "Speckley_ReducedNodes [ReducedContinuousFunction(domain)]";
        case Elements: return "Speckley_Elements [Function(domain)]";
        case ReducedElements: return "Speckley_ReducedElements [ReducedFunction(domain)]";
        case FaceElements: return "Speckley_FaceElements [FunctionOnBoundary(domain)]";
        case ReducedFaceElements: return "Speckley_ReducedFaceElements [ReducedFunctionOnBoundary(domain)]";
        case Points: return "Speckley_Points [DiracDeltaFunctions(domain)]";
        default:
            break;
    }
    return "Invalid function space type code";
}

pair<int,dim_t> SpeckleyDomain::getDataShape(int fsType) const
{
    int ptsPerSample = (m_order + 1) * (m_order + 1);
    if (m_numDim==3)
        ptsPerSample *= (m_order + 1);
    switch (fsType) {
        case Nodes:
        case ReducedNodes: //FIXME: reduced
            return pair<int,dim_t>(1, getNumNodes());
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom: //FIXME: reduced
            return pair<int,dim_t>(1, getNumDOF());
        case Elements:
            return pair<int,dim_t>(ptsPerSample, getNumElements());
        case FaceElements:
            return pair<int,dim_t>(ptsPerSample/2, getNumFaceElements());
        case ReducedElements:
            return pair<int,dim_t>(1, getNumElements());
        case ReducedFaceElements:
            return pair<int,dim_t>(1, getNumFaceElements());
        case Points:
            return pair<int,dim_t>(1, m_diracPoints.size());
        default:
            break;
    }

    stringstream msg;
    msg << "getDataShape: Invalid function space type " << fsType
        << " for " << getDescription();
    throw SpeckleyException(msg.str());
}

string SpeckleyDomain::showTagNames() const
{
    stringstream ret;
    TagMap::const_iterator it;
    for (it=m_tagMap.begin(); it!=m_tagMap.end(); it++) {
        if (it!=m_tagMap.begin()) ret << ", ";
        ret << it->first;
    }
    return ret.str();
}

bool SpeckleyDomain::commonFunctionSpace(const vector<int>& fs, int& resultcode) const
{
   /*
    The idea is to use equivalence classes (i.e. types which can be
    interpolated back and forth):
    class 0: DOF <-> Nodes
    class 1: ReducedDOF <-> ReducedNodes
    class 2: Points
    class 3: Elements
    class 4: ReducedElements
    class 5: FaceElements
    class 6: ReducedFaceElements

    There is also a set of lines. Interpolation is possible down a line but not
    between lines.
    class 0 and 1 belong to all lines so aren't considered.
    line 0: class 2
    line 1: class 3,4
    line 2: class 5,6

    For classes with multiple members (eg class 1) we have vars to record if
    there is at least one instance. e.g. hasnodes is true if we have at least
    one instance of Nodes.
    */
    if (fs.empty())
        return false;
    vector<bool> hasclass(7, false);
    vector<int> hasline(3, 0);
    bool hasnodes=false;
    bool hasrednodes=false;
    for (size_t i=0; i<fs.size(); ++i) {
        switch (fs[i]) {
            case Nodes: hasnodes=true; // fall through
            case DegreesOfFreedom:
                hasclass[0]=true;
                break;
            case ReducedNodes: hasrednodes=true; // fall through
            case ReducedDegreesOfFreedom:
                hasclass[1]=true;
                break;
            case Points:
                hasline[0]=1;
                hasclass[2]=true;
                break;
            case Elements:
                hasclass[3]=true;
                hasline[1]=1;
                break;
            case ReducedElements:
                hasclass[4]=true;
                hasline[1]=1;
                break;
            case FaceElements:
                hasclass[5]=true;
                hasline[2]=1;
                break;
            case ReducedFaceElements:
                hasclass[6]=true;
                hasline[2]=1;
                break;
            default:
                return false;
        }
    }
    int numLines=hasline[0]+hasline[1]+hasline[2];

    // fail if we have more than one leaf group
    // = there are at least two branches we can't interpolate between
    if (numLines > 1)
        return false;
    else if (numLines==1) {
        // we have points
        if (hasline[0]==1)
            resultcode=Points;
        else if (hasline[1]==1) {
            if (hasclass[4])
                resultcode=ReducedElements;
            else
                resultcode=Elements;
        } else { // hasline[2]==1
            if (hasclass[6])
                resultcode=ReducedFaceElements;
            else
                resultcode=FaceElements;
        }
    } else { // numLines==0
        if (hasclass[1])
            // something from class 1
            resultcode=(hasrednodes ? ReducedNodes : ReducedDegreesOfFreedom);
        else // something from class 0
            resultcode=(hasnodes ? Nodes : DegreesOfFreedom);
    }
    return true;
}

bool SpeckleyDomain::probeInterpolationOnDomain(int fsType_source,
                                              int fsType_target) const
{
    if (!isValidFunctionSpaceType(fsType_target)) {
        stringstream msg;
        msg << "probeInterpolationOnDomain: Invalid function space type "
            << fsType_target << " for " << getDescription();
        throw SpeckleyException(msg.str());
    }

    switch (fsType_source) {
        case Nodes:
        case DegreesOfFreedom:
            return true;
        case ReducedNodes:
        case ReducedDegreesOfFreedom:
            return (fsType_target != Nodes &&
                    fsType_target != DegreesOfFreedom);
        case Elements:
        case ReducedElements:
            return (fsType_target==Elements ||
                    fsType_target==ReducedElements ||
                    fsType_target==Nodes);
        case FaceElements:
        case ReducedFaceElements:
            return (fsType_target==FaceElements ||
                    fsType_target==ReducedFaceElements);
        case Points:
            return (fsType_target==fsType_source);

        default: {
            stringstream msg;
            msg << "probeInterpolationOnDomain: Invalid function space type "
                << fsType_source << " for " << getDescription();
            throw SpeckleyException(msg.str());
        }
    }
}

signed char SpeckleyDomain::preferredInterpolationOnDomain(int fsType_source,
                                              int fsType_target) const
{
    if (!isValidFunctionSpaceType(fsType_target)) {
        stringstream msg;
        msg << "preferredInterpolationOnDomain: Invalid function space type "
            << fsType_target << " for " << getDescription();
        throw SpeckleyException(msg.str());
    }

    if (fsType_source==fsType_target) {
        return 1;
    }
    // There is a complication here in that Nodes and DegreesOfFreedom
    // can be interpolated to anything, so we need to handle the
    // reverse case for them specially

    if ((fsType_target==Nodes) || (fsType_target==DegreesOfFreedom)) {
        return -1;  // reverse interpolation
    }

    switch (fsType_source) {
        case Nodes:
        case DegreesOfFreedom:
            return 1;
        case ReducedNodes:
        case ReducedDegreesOfFreedom:
            return (fsType_target != Nodes &&
                    fsType_target != DegreesOfFreedom) ? -1 : 0;
        case Elements:
            return (fsType_target==ReducedElements) ? 1 : 0;
        case ReducedElements:
            return (fsType_target==Elements) ? -1 : 0;
        case FaceElements:
            return (fsType_target==ReducedFaceElements) ? 1 : 0;
        case ReducedFaceElements:
            return (fsType_target==FaceElements) ? -1 : 0;
        case Points:
            return false;  // other case caught by the if above
        default: {
            stringstream msg;
            msg << "probeInterpolationOnDomain: Invalid function space type "
                << fsType_source << " for " << getDescription();
            throw SpeckleyException(msg.str());
        }
    }
}

void SpeckleyDomain::interpolateOnDomain(escript::Data& target,
                                       const escript::Data& in) const
{
    const SpeckleyDomain& inDomain=dynamic_cast<const SpeckleyDomain&>(*(in.getFunctionSpace().getDomain()));
    const SpeckleyDomain& targetDomain=dynamic_cast<const SpeckleyDomain&>(*(target.getFunctionSpace().getDomain()));
    if (inDomain != *this)
        throw SpeckleyException("Illegal domain of interpolant");
    if (targetDomain != *this)
        throw SpeckleyException("Illegal domain of interpolation target");

    stringstream msg;
    msg << "interpolateOnDomain() not implemented for function space "
        << functionSpaceTypeAsString(in.getFunctionSpace().getTypeCode())
        << " -> "
        << functionSpaceTypeAsString(target.getFunctionSpace().getTypeCode());

    const int inFS = in.getFunctionSpace().getTypeCode();
    const int outFS = target.getFunctionSpace().getTypeCode();

    // simplest case: 1:1 copy
    if (inFS==outFS) {
        copyData(target, in);
    // not allowed: reduced nodes/DOF->non-reduced nodes/DOF
    } else if ((inFS==ReducedNodes || inFS==ReducedDegreesOfFreedom)
            && (outFS==Nodes || outFS==DegreesOfFreedom)) {
        throw SpeckleyException("interpolateOnDomain: Cannot interpolate "
                              "reduced data to non-reduced data.");
    } else if ((inFS==Elements && outFS==ReducedElements)
            || (inFS==FaceElements && outFS==ReducedFaceElements)) {
        if (in.actsExpanded())
            averageData(target, in);
        else
            copyData(target, in);
    } else if (inFS==Elements && outFS==Nodes) {
            interpolateElementsOnNodes(target, in, false);
    } else if ((inFS==ReducedElements && outFS==Elements)
            || (inFS==ReducedFaceElements && outFS==FaceElements)) {
        multiplyData(target, in);
    } else {
        switch (inFS) {
            case Nodes:
            case ReducedNodes:
                switch (outFS) {
                    case DegreesOfFreedom:
                    case ReducedDegreesOfFreedom:
                        if (getMPISize()==1)
                            copyData(target, in);
                        else
                            throw SpeckleyException("interpolateOnDomain(): MPI not supported");
                        break;

                    case Nodes:
                    case ReducedNodes:
                        copyData(target, in);
                        break;
                    case Elements:
                        interpolateNodesOnElements(target, in, false);
                        break;

                    case ReducedElements:
                        interpolateNodesOnElements(target, in, true);
                        break;

                    case FaceElements:
                        interpolateNodesOnFaces(target, in, false);
                        break;

                    case ReducedFaceElements:
                        interpolateNodesOnFaces(target, in, true);
                        break;
                    case Points:
                        {
                            const dim_t numComp = in.getDataPointSize();
                            const int nDirac = m_diracPoints.size();
                            target.requireWrite();
#pragma omp parallel for
                            for (int i = 0; i < nDirac; i++) {
                                const double* src = in.getSampleDataRO(m_diracPoints[i].node);
                                copy(src, src+numComp, target.getSampleDataRW(i));
                            }
                        }
                        break;
                    default:
                        throw SpeckleyException(msg.str());
                }
                break;

            case DegreesOfFreedom:
            case ReducedDegreesOfFreedom:
                switch (outFS) {
                    case Nodes:
                    case ReducedNodes:
                        if (getMPISize()==1)
                            copyData(target, in);
                        else
                            throw SpeckleyException("interpolateOnDomain(): MPI not supported");
                        break;

                    case DegreesOfFreedom:
                    case ReducedDegreesOfFreedom:
                        copyData(target, in);
                        break;

                    case Elements:
                    case ReducedElements:
                        if (getMPISize()==1) {
                            interpolateNodesOnElements(target, in, outFS==ReducedElements);
                        } else {
                            escript::Data contIn(in, (inFS==DegreesOfFreedom ?
                                        escript::continuousFunction(*this) : escript::reducedContinuousFunction(*this)));
                            interpolateNodesOnElements(target, contIn, outFS==ReducedElements);
                        }
                        break;
                    case FaceElements:
                    case ReducedFaceElements:
                        if (getMPISize()==1) {
                            interpolateNodesOnFaces(target, in, outFS==ReducedFaceElements);
                        } else {
                            escript::Data contIn(in, (inFS==DegreesOfFreedom ?
                                     escript::continuousFunction(*this) : escript::reducedContinuousFunction(*this)));
                            interpolateNodesOnFaces(target, contIn, outFS==ReducedFaceElements);
                        }
                        break;

                    default:
                        throw SpeckleyException(msg.str());
                }
                break;
            default:
                throw SpeckleyException(msg.str());
        }
    }
}

escript::Data SpeckleyDomain::getX() const
{
    return escript::continuousFunction(*this).getX();
}

escript::Data SpeckleyDomain::getNormal() const
{
    return escript::functionOnBoundary(*this).getNormal();
}

escript::Data SpeckleyDomain::getSize() const
{
    return escript::function(*this).getSize();
}

void SpeckleyDomain::setToX(escript::Data& arg) const
{
    const SpeckleyDomain& argDomain=dynamic_cast<const SpeckleyDomain&>(
            *(arg.getFunctionSpace().getDomain()));
    if (argDomain != *this)
        throw SpeckleyException("setToX: Illegal domain of data point locations");
    if (!arg.isExpanded())
        throw SpeckleyException("setToX: Expanded Data object expected");

    if (arg.getFunctionSpace().getTypeCode()==Nodes) {
        assembleCoordinates(arg);
    } else {
        // interpolate the result
        escript::Data contData=escript::Vector(0., escript::continuousFunction(*this), true);
        assembleCoordinates(contData);
        interpolateOnDomain(arg, contData);
    }
}

void SpeckleyDomain::setToGradient(escript::Data& grad, const escript::Data& arg) const
{
    const SpeckleyDomain& argDomain=dynamic_cast<const SpeckleyDomain&>(
            *(arg.getFunctionSpace().getDomain()));
    if (argDomain != *this)
        throw SpeckleyException("setToGradient: Illegal domain of gradient argument");
    const SpeckleyDomain& gradDomain=dynamic_cast<const SpeckleyDomain&>(
            *(grad.getFunctionSpace().getDomain()));
    if (gradDomain != *this)
        throw SpeckleyException("setToGradient: Illegal domain of gradient");

    switch (grad.getFunctionSpace().getTypeCode()) {
        case Elements:
        case ReducedElements:
        case FaceElements:
        case ReducedFaceElements:
            break;
        default: {
            stringstream msg;
            msg << "setToGradient: not supported for "
                << functionSpaceTypeAsString(grad.getFunctionSpace().getTypeCode());
            throw SpeckleyException(msg.str());
        }
    }

    switch (arg.getFunctionSpace().getTypeCode()) {
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
        case Nodes:
        case ReducedNodes:
        case Elements:
        case ReducedElements:
            break;
        default: {
            throw SpeckleyException("setToGradient: only supported for nodal input data");
        }
    }

    if (getMPISize() > 1) {
        if (arg.getFunctionSpace().getTypeCode()==DegreesOfFreedom) {
            escript::Data contArg(arg, escript::continuousFunction(*this));
            assembleGradient(grad, contArg);
        } else if (arg.getFunctionSpace().getTypeCode()==ReducedDegreesOfFreedom) {
            escript::Data contArg(arg, escript::reducedContinuousFunction(*this));
            assembleGradient(grad, contArg);
        } else {
            assembleGradient(grad, arg);
        }
    } else {
        assembleGradient(grad, arg);
    }
}

void SpeckleyDomain::setToIntegrals(vector<double>& integrals, const escript::Data& arg) const
{
    const SpeckleyDomain& argDomain=dynamic_cast<const SpeckleyDomain&>(
            *(arg.getFunctionSpace().getDomain()));
    if (argDomain != *this)
        throw SpeckleyException("setToIntegrals: illegal domain of integration kernel");

    switch (arg.getFunctionSpace().getTypeCode()) {
        case Nodes:
        case ReducedNodes:
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
            {
                escript::Data funcArg(arg, escript::function(*this));
                assembleIntegrate(integrals, funcArg);
            }
            break;
        case Elements:
        case ReducedElements:
        case FaceElements:
        case ReducedFaceElements:
            assembleIntegrate(integrals, arg);
            break;
        default: {
            stringstream msg;
            msg << "setToIntegrals: not supported for "
                << functionSpaceTypeAsString(arg.getFunctionSpace().getTypeCode());
            throw SpeckleyException(msg.str());
        }
    }

}

bool SpeckleyDomain::isCellOriented(int fsType) const
{
    switch(fsType) {
        case Nodes:
        case ReducedNodes:
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
    msg << "isCellOriented: invalid function space type " << fsType
        << " on " << getDescription();
    throw SpeckleyException(msg.str());
}

bool SpeckleyDomain::canTag(int fsType) const
{
    switch(fsType) {
        case Nodes:
        case Elements:
        case ReducedElements:
        case FaceElements:
        case Points:
        case ReducedFaceElements:
            return true;
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
        case ReducedNodes:
            return false;
        default:
            break;
    }
    stringstream msg;
    msg << "canTag: invalid function space type " << fsType << " on "
        << getDescription();
    throw SpeckleyException(msg.str());
}

void SpeckleyDomain::setTags(const int fsType, const int newTag, const escript::Data& mask) const
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
            msg << "setTags: invalid function space type " << fsType;
            throw SpeckleyException(msg.str());
        }
    }
    if (target->size() != num) {
        target->assign(num, -1);
    }

#pragma omp parallel for
    for (index_t i=0; i<num; i++) {
        if (mask.getSampleDataRO(i)[0] > 0) {
            (*target)[i]=newTag;
        }
    }
    updateTagsInUse(fsType);
}

int SpeckleyDomain::getTagFromSampleNo(int fsType, int sampleNo) const
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
        case Points:
            if (m_diracPoints.size() > sampleNo)
                return m_diracPoints[sampleNo].tag;
            break;
        case FaceElements:
        case ReducedFaceElements:
            if (m_faceTags.size() > sampleNo)
                return m_faceTags[sampleNo];
            break;
        default: {
            stringstream msg;
            msg << "getTagFromSampleNo: invalid function space type " << fsType;
            throw SpeckleyException(msg.str());
        }
    }
    return -1;
}

int SpeckleyDomain::getNumberOfTagsInUse(int fsType) const
{
    switch(fsType) {
        case Points:
            return m_diracPointNodeIDs.size();
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
            msg << "getNumberOfTagsInUse: invalid function space type "
                << fsType;
            throw SpeckleyException(msg.str());
        }
    }
}

const int* SpeckleyDomain::borrowListOfTagsInUse(int fsType) const
{
    switch(fsType) {
        case Points:
            return &m_diracPointNodeIDs[0];
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
            msg << "borrowListOfTagsInUse: invalid function space type "
                << fsType;
            throw SpeckleyException(msg.str());
        }
    }
}

void SpeckleyDomain::Print_Mesh_Info(bool full) const
{
    cout << "Print_Mesh_Info for " << getDescription() << " running on CPU "
        << m_mpiInfo->rank << ". MPI size: " << m_mpiInfo->size << endl;
    cout << "Number of dimensions: " << m_numDim << endl;
    cout << "Number of elements per rank: " << getNumElements() << endl;
    // write tags
    if (m_tagMap.size() > 0) {
        cout << "Tags:" << endl;
        TagMap::const_iterator it;
        for (it=m_tagMap.begin(); it!=m_tagMap.end(); it++) {
            cout << "  " << setw(5) << it->second << " " << it->first << endl;
        }
    }
}

int SpeckleyDomain::getSystemMatrixTypeId(int solver, int preconditioner,
                                        int package, bool symmetry) const
{
	throw SpeckleyException("System matrices not supported by Speckley");
}

int SpeckleyDomain::getTransportTypeId(int solver, int preconditioner,
                                     int package, bool symmetry) const
{
    throw SpeckleyException("Transport problems not supported by Speckley");
}

escript::ASM_ptr SpeckleyDomain::newSystemMatrix(int row_blocksize,
        const escript::FunctionSpace& row_functionspace, int column_blocksize,
        const escript::FunctionSpace& column_functionspace, int type) const
{
	throw SpeckleyException("Speckley domains do not support system matrices");
}

void SpeckleyDomain::addToSystem(escript::AbstractSystemMatrix& mat,
                               escript::Data& rhs, const DataMap& coefs,
                               Assembler_ptr assembler) const
{
   throw SpeckleyException("Speckley domains do not support system matrices");
}

void SpeckleyDomain::addToSystemFromPython(escript::AbstractSystemMatrix& mat,
                                         escript::Data& rhs,
                                         const bp::list& data,
                                         Assembler_ptr assembler) const
{
    DataMap mapping;
    tupleListToMap(mapping, data);
    addToSystem(mat, rhs, mapping, assembler);
}

Assembler_ptr SpeckleyDomain::createAssemblerFromPython(const string type,
                                                const bp::list& options) const
{
    DataMap mapping;
    tupleListToMap(mapping, options);
    return createAssembler(type, mapping);
}

void SpeckleyDomain::addToRHSFromPython(escript::Data& rhs,
            const bp::list& data, Assembler_ptr assembler) const
{
    DataMap mapping;
    tupleListToMap(mapping, data);
    rhs.expand(); //yes, this is absolutely neccessary
    addToRHS(rhs, mapping, assembler);
}

void SpeckleyDomain::addToRHS(escript::Data& rhs, const DataMap& coefs,
                            Assembler_ptr assembler) const
{
    if (isNotEmpty("y_contact", coefs))
        throw SpeckleyException(
                    "addPDEToRHS: Speckley does not support contact elements");

    if (rhs.isEmpty()) {
        if (isNotEmpty("X", coefs) || isNotEmpty("Y", coefs))
            throw SpeckleyException(
                    "addPDEToRHS: right hand side coefficients are provided "
                    "but no right hand side vector given");
        else
            return;
    }

    assemblePDE(NULL, rhs, coefs, assembler);
    assemblePDEBoundary(NULL, rhs, coefs, assembler);
    assemblePDEDirac(NULL, rhs, coefs, assembler);
}

escript::ATP_ptr SpeckleyDomain::newTransportProblem(int blocksize,
        const escript::FunctionSpace& functionspace, int type) const
{
    throw SpeckleyException("Speckley domains do not support transport problems");
}

void SpeckleyDomain::addPDEToTransportProblemFromPython(
        escript::AbstractTransportProblem& tp, escript::Data& source,
        const bp::list& data, Assembler_ptr assembler) const
{
    DataMap mapping;
    tupleListToMap(mapping, data);
    addPDEToTransportProblem(tp, source, mapping, assembler);
}

void SpeckleyDomain::addPDEToTransportProblem(
        escript::AbstractTransportProblem& tp,
        escript::Data& source, const DataMap& coefs,
        Assembler_ptr assembler) const
{
    throw SpeckleyException("Speckley domains do not support transport problems");
}

void SpeckleyDomain::setNewX(const escript::Data& arg)
{
    throw SpeckleyException("setNewX(): operation not supported");
}

//protected
void SpeckleyDomain::copyData(escript::Data& out, const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    const dim_t numSamples = in.getNumSamples();
    out.requireWrite();
#pragma omp parallel for
    for (index_t i=0; i<numSamples; i++) {
        const double* src = in.getSampleDataRO(i);
        copy(src, src+numComp, out.getSampleDataRW(i));
    }
}

//protected
void SpeckleyDomain::averageData(escript::Data& out, const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    const dim_t dpp = in.getNumDataPointsPerSample();
    const dim_t numSamples = in.getNumSamples();
    out.requireWrite();
#pragma omp parallel for
    for (index_t i=0; i<numSamples; i++) {
        const double* src = in.getSampleDataRO(i);
        double* dest = out.getSampleDataRW(i);
        for (index_t c=0; c<numComp; c++) {
            double res=0.;
            for (index_t q=0; q<dpp; q++)
                res+=src[c+q*numComp];
            *dest++ = res/dpp;
        }
    }
}

//protected
void SpeckleyDomain::multiplyData(escript::Data& out, const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    const dim_t dpp = out.getNumDataPointsPerSample();
    const dim_t numSamples = in.getNumSamples();
    out.requireWrite();
#pragma omp parallel for
    for (index_t i=0; i<numSamples; i++) {
        const double* src = in.getSampleDataRO(i);
        double* dest = out.getSampleDataRW(i);
        for (index_t c=0; c<numComp; c++) {
            for (index_t q=0; q<dpp; q++)
                dest[c+q*numComp] = src[c];
        }
    }
}

//protected
void SpeckleyDomain::updateTagsInUse(int fsType) const
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
        case Points:
            throw SpeckleyException("updateTagsInUse for Speckley dirac points "
                    "not supported");
        default:
            return;
    }

    // gather global unique tag values from tags into tagsInUse
    tagsInUse->clear();
    index_t lastFoundValue = INDEX_T_MIN, minFoundValue, local_minFoundValue;
    const long numTags = tags->size();

    while (true) {
        // find smallest value bigger than lastFoundValue
        minFoundValue = INDEX_T_MAX;
#pragma omp parallel private(local_minFoundValue)
        {
            local_minFoundValue = minFoundValue;
            long i;     // should be size_t but omp mutter mutter
#pragma omp for schedule(static) private(i) nowait
            for (i = 0; i < numTags; i++) {
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

//protected
void SpeckleyDomain::addToSystemMatrix(escript::AbstractSystemMatrix* mat,
                                     const IndexVector& nodes, dim_t numEq,
                                     const DoubleVector& array) const
{
    throw SpeckleyException("addToSystemMatrix not yet implemented");
}

//private
void SpeckleyDomain::assemblePDE(escript::AbstractSystemMatrix* mat,
                               escript::Data& rhs, const DataMap& coefs,
                               Assembler_ptr assembler) const
{
    if (rhs.isEmpty() && isNotEmpty("X", coefs) && isNotEmpty("Y", coefs))
        throw SpeckleyException("assemblePDE: right hand side coefficients are "
                    "provided but no right hand side vector given");

    vector<int> fsTypes;
    assembler->collateFunctionSpaceTypes(fsTypes, coefs);

    if (fsTypes.empty()) {
        return;
    }

    int fs=fsTypes[0];
    if (fs != Elements)
        throw SpeckleyException("assemblePDE: illegal function space type for coefficients");

    for (vector<int>::const_iterator it=fsTypes.begin()+1; it!=fsTypes.end(); it++) {
        if (*it != fs) {
            throw SpeckleyException("assemblePDE: coefficient function spaces don't match");
        }
    }

    int numEq, numComp;
    if (!mat) {
        if (rhs.isEmpty()) {
            numEq = numComp = 1;
        } else {
            numEq = numComp = rhs.getDataPointSize();
        }
    } else {
        if (!rhs.isEmpty() && rhs.getDataPointSize() != mat->getRowBlockSize())
            throw SpeckleyException("assemblePDE: matrix row block size and number of components of right hand side don't match");
        numEq = mat->getRowBlockSize();
        numComp = mat->getColumnBlockSize();
    }

    if (numEq != numComp)
        throw SpeckleyException("assemblePDE: number of equations and number of solutions don't match");

    //TODO: check shape and num samples of coeffs

    if (numEq==1) {
        if (fs==ReducedElements) {
            assembler->assemblePDESingleReduced(mat, rhs, coefs);
        } else {
            assembler->assemblePDESingle(mat, rhs, coefs);
        }
    } else {
        if (fs==ReducedElements) {
            assembler->assemblePDESystemReduced(mat, rhs, coefs);
        } else {
            assembler->assemblePDESystem(mat, rhs, coefs);
        }
    }
}

//private
void SpeckleyDomain::assemblePDEBoundary(escript::AbstractSystemMatrix* mat,
        escript::Data& rhs, const DataMap& coefs, Assembler_ptr assembler) const
{
    if (rhs.isEmpty() && isNotEmpty("y", coefs))
        throw SpeckleyException("assemblePDEBoundary: y provided but no right hand side vector given");

    int fs=-1;
    if (isNotEmpty("d", coefs))
        fs=coefs.find("d")->second.getFunctionSpace().getTypeCode();
    if (isNotEmpty("y", coefs)) {
        DataMap::const_iterator iy = coefs.find("y");
        if (fs == -1)
            fs = iy->second.getFunctionSpace().getTypeCode();
        else if (fs != iy->second.getFunctionSpace().getTypeCode())
            throw SpeckleyException("assemblePDEBoundary: coefficient function spaces don't match");
    }
    if (fs==-1) {
        return;
    }

    if (fs != FaceElements && fs != ReducedFaceElements)
        throw SpeckleyException("assemblePDEBoundary: illegal function space type for coefficients");

    int numEq = 1, numComp = 1;
    if (!mat) {
        if (!rhs.isEmpty()) {
            numEq = numComp = rhs.getDataPointSize();
        }
    } else {
        if (!rhs.isEmpty() && rhs.getDataPointSize() != mat->getRowBlockSize())
            throw SpeckleyException("assemblePDEBoundary: matrix row block size and number of components of right hand side don't match");
        numEq = mat->getRowBlockSize();
        numComp = mat->getColumnBlockSize();
    }

    if (numEq != numComp)
        throw SpeckleyException("assemblePDEBoundary: number of equations and number of solutions don't match");

    //TODO: check shape and num samples of coeffs

    if (numEq==1) {
        if (fs==ReducedFaceElements)
            assembler->assemblePDEBoundarySingleReduced(mat, rhs, coefs);
        else
            assembler->assemblePDEBoundarySingle(mat, rhs, coefs);
    } else {
        if (fs==ReducedFaceElements)
            assembler->assemblePDEBoundarySystemReduced(mat, rhs, coefs);
        else
            assembler->assemblePDEBoundarySystem(mat, rhs, coefs);
    }
}

void SpeckleyDomain::assemblePDEDirac(escript::AbstractSystemMatrix* mat,
                                    escript::Data& rhs, const DataMap& coefs,
                                    Assembler_ptr assembler) const
{
    bool yNotEmpty = isNotEmpty("y_dirac", coefs);
    bool dNotEmpty = isNotEmpty("d_dirac", coefs);
    if (!(yNotEmpty || dNotEmpty)) {
        return;
    }
    escript::Data d = unpackData("d_dirac", coefs);
    escript::Data y = unpackData("y_dirac", coefs);
    int nEq = 1, nComp = 1;
    if (!mat) {
        if (!rhs.isEmpty()) {
            nEq = nComp = rhs.getDataPointSize();
        }
    } else {
        if (!rhs.isEmpty() && rhs.getDataPointSize()!=mat->getRowBlockSize())
            throw SpeckleyException("assemblePDEDirac: matrix row block size "
                    "and number of components of right hand side don't match");
        nEq = mat->getRowBlockSize();
        nComp = mat->getColumnBlockSize();
    }

    for (int i = 0; i < m_diracPoints.size(); i++) { //only for this rank
        const IndexVector rowIndex(1, getDofOfNode(m_diracPoints[i].node));
        if (yNotEmpty) {
            const double *EM_F = y.getSampleDataRO(i);
            double *F_p = rhs.getSampleDataRW(0);
            if (rowIndex[0] < getNumDOF()) {
                for (index_t eq = 0; eq < nEq; eq++) {
                    F_p[INDEX2(eq, rowIndex[0], nEq)] += EM_F[INDEX2(eq,i,nEq)];
                }
            }
        }
        if (dNotEmpty) {
            throw SpeckleyException("Rectangle::assemblePDEDirac can't cope with d not being empty right now");
//            const double *EM_S = d.getSampleDataRO(i);
//            std::vector<double> contents(EM_S,
//                        EM_S+nEq*nEq*nComp*rowIndex.size()); //TODO: this will break with MPI
//            addToSystemMatrix(mat, rowIndex, nEq, rowIndex, nComp, contents);
        }
    }
}

bool SpeckleyDomain::probeInterpolationACross(int fsType_source,
        const escript::AbstractDomain&, int fsType_target) const
{
    return false;
}

void SpeckleyDomain::interpolateACross(escript::Data& target, const escript::Data& source) const
{
    throw SpeckleyException("interpolateACross() not supported");
}

// Expecting ("gaussian", radius, sigma)
bool SpeckleyDomain::supportsFilter(const bp::tuple& t) const
{
    if (len(t) == 0) { // so we can handle unfiltered randoms
        return true;
    }
    return false;
}

void SpeckleyDomain::addPoints(const vector<double>& coords,
        const vector<int>& tags)
{
    for (int i = 0; i < tags.size(); i++) {
        dim_t node = findNode(&coords[i * m_numDim]);
        if (node >= 0) {
            m_diracPointNodeIDs.push_back(borrowSampleReferenceIDs(Nodes)[node]);
            DiracPoint dp;
            dp.node = node; //local
            dp.tag = tags[i];
            m_diracPoints.push_back(dp);
        }
    }
}

} // end of namespace speckley


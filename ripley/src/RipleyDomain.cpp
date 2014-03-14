
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

bool isNotEmpty(string target, map<string, escript::Data> mapping) {
    map<string, escript::Data>::iterator i = mapping.find(target);
    return i != mapping.end() && !i->second.isEmpty();
}

void tupleListToMap(map<string, escript::Data>& mapping,
        boost::python::list& list) {
    using boost::python::tuple;
    using boost::python::extract;
    for (int i = 0; i < len(list); i++) {
        if (!extract<tuple>(list[i]).check())
            throw RipleyException("Passed in list contains objects"
                                    " other than tuples");
        tuple t = extract<tuple>(list[i]);
        if (len(t) != 2 || !extract<std::string>(t[0]).check() ||
                !extract<escript::Data>(t[1]).check())
            throw RipleyException("The passed in list must contain tuples"
                " of the form (string, escript::data)");
        mapping[extract<std::string>(t[0])] = extract<escript::Data>(t[1]);
    }
}

RipleyDomain::RipleyDomain(dim_t dim, escript::SubWorld_ptr p) :
    m_numDim(dim),
    m_status(0)
{
    if (p.get()==0)	
    {
	m_mpiInfo = makeInfo(MPI_COMM_WORLD);
    }
    else
    {
	m_mpiInfo = p->getMPI();
    }
}

RipleyDomain::~RipleyDomain()
{
    // cleanup of MPI is dealt with by shared_ptr
}

bool RipleyDomain::operator==(const AbstractDomain& other) const
{
    const RipleyDomain* o=dynamic_cast<const RipleyDomain*>(&other);
    if (o) {
        return (m_tagMap==o->m_tagMap && m_nodeTags==o->m_nodeTags
                && m_elementTags==o->m_elementTags
                && m_faceTags==o->m_faceTags);
    }
    return false;
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
        case ReducedNodes: return "Ripley_ReducedNodes";
        case Elements: return "Ripley_Elements";
        case ReducedElements: return "Ripley_ReducedElements";
        case FaceElements: return "Ripley_FaceElements";
        case ReducedFaceElements: return "Ripley_ReducedFaceElements";
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
        case ReducedNodes: //FIXME: reduced
            return pair<int,int>(1, getNumNodes());
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom: //FIXME: reduced
            return pair<int,int>(1, getNumDOF());
        case Elements:
            return pair<int,int>(ptsPerSample, getNumElements());
        case FaceElements:
            return pair<int,int>(ptsPerSample/2, getNumFaceElements());
        case ReducedElements:
            return pair<int,int>(1, getNumElements());
        case ReducedFaceElements:
            return pair<int,int>(1, getNumFaceElements());
        case Points:
            return pair<int,int>(1, m_diracPoints.size());
        default:
            break;
    }

    stringstream msg;
    msg << "getDataShape: Invalid function space type " << fsType
        << " for " << getDescription();
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

bool RipleyDomain::probeInterpolationOnDomain(int fsType_source,
                                              int fsType_target) const
{
    if (!isValidFunctionSpaceType(fsType_target)) {
        stringstream msg;
        msg << "probeInterpolationOnDomain: Invalid function space type "
            << fsType_target << " for " << getDescription();
        throw RipleyException(msg.str());
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
                    fsType_target==ReducedElements);
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
            throw RipleyException(msg.str());
        }
    }
}

signed char RipleyDomain::preferredInterpolationOnDomain(int fsType_source,
                                              int fsType_target) const
{
    if (!isValidFunctionSpaceType(fsType_target)) {
        stringstream msg;
        msg << "preferredInterpolationOnDomain: Invalid function space type "
            << fsType_target << " for " << getDescription();
        throw RipleyException(msg.str());
    }

    if (fsType_source==fsType_target) {
        return 1;
    }
    // There is a complication here in that Nodes and DegreesOfFreedom
    // can be interpolated to anything, so we need to handle the
    // reverse case for them specially

    if ((fsType_target==Nodes) || (fsType_target==DegreesOfFreedom))
    {
        return -1;    // reverse interpolation
    }

    switch (fsType_source) {
        case Nodes:
        case DegreesOfFreedom:
            return 1;
        case ReducedNodes:
        case ReducedDegreesOfFreedom:
            return (fsType_target != Nodes &&
                    fsType_target != DegreesOfFreedom)?-1:0;
        case Elements:
	    return (fsType_target==ReducedElements)?1:0;
        case ReducedElements:
	    return (fsType_target==Elements)?-1:0;	  
        case FaceElements:
	    return (fsType_target==ReducedFaceElements)?1:0;
        case ReducedFaceElements:
            return (fsType_target==FaceElements)?-1:0;
	case Points: 
	    return false;	// other case caught by the if above

        default: {
            stringstream msg;
            msg << "probeInterpolationOnDomain: Invalid function space type "
                << fsType_source << " for " << getDescription();
            throw RipleyException(msg.str());
        }
    }
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

    const int inFS = in.getFunctionSpace().getTypeCode();
    const int outFS = target.getFunctionSpace().getTypeCode();

    // simplest case: 1:1 copy
    if (inFS==outFS) {
        copyData(target, in);
    // not allowed: reduced nodes/dof->non-reduced nodes/dof
    } else if ((inFS==ReducedNodes || inFS==ReducedDegreesOfFreedom)
            && (outFS==Nodes || outFS==DegreesOfFreedom)) {
        throw RipleyException("interpolateOnDomain: Cannot interpolate reduced data to non-reduced data.");
    } else if ((inFS==Elements && outFS==ReducedElements)
            || (inFS==FaceElements && outFS==ReducedFaceElements)) {
        if (in.actsExpanded())
            averageData(target, in);
        else
            copyData(target, in);
    } else if ((inFS==ReducedElements && outFS==Elements)
            || (inFS==ReducedFaceElements && outFS==FaceElements)) {
        multiplyData(target, in);
    } else {
        switch (inFS) {
            case Nodes:
            case ReducedNodes: //FIXME: reduced
                switch (outFS) {
                    case DegreesOfFreedom:
                    case ReducedDegreesOfFreedom: //FIXME: reduced
                        if (getMPISize()==1)
                            copyData(target, in);
                        else
                            nodesToDOF(target, in);
                        break;

                    case Nodes:
                    case ReducedNodes: //FIXME: reduced
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
                            target.requireWrite();
                        #pragma omp parallel for
                            for (int i = 0; i < m_diracPoints.size(); i++) {
                                const double* src = in.getSampleDataRO(m_diracPoints[i].node);
                                copy(src, src+numComp, target.getSampleDataRW(i));
                            }
                        }
                        break;
                    default:
                        throw RipleyException(msg.str());
                }
                break;

            case DegreesOfFreedom:
            case ReducedDegreesOfFreedom: //FIXME: reduced
                switch (outFS) {
                    case Nodes:
                    case ReducedNodes: //FIXME: reduced
                        if (getMPISize()==1)
                            copyData(target, in);
                        else
                            dofToNodes(target, in);
                        break;

                    case DegreesOfFreedom:
                    case ReducedDegreesOfFreedom: //FIXME: reduced
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
                        throw RipleyException(msg.str());
                }
                break;
            default:
                throw RipleyException(msg.str());
        }
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

void RipleyDomain::setToGradient(escript::Data& grad, const escript::Data& arg) const
{
    const RipleyDomain& argDomain=dynamic_cast<const RipleyDomain&>(
            *(arg.getFunctionSpace().getDomain()));
    if (argDomain != *this)
        throw RipleyException("setToGradient: Illegal domain of gradient argument");
    const RipleyDomain& gradDomain=dynamic_cast<const RipleyDomain&>(
            *(grad.getFunctionSpace().getDomain()));
    if (gradDomain != *this)
        throw RipleyException("setToGradient: Illegal domain of gradient");

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
            throw RipleyException(msg.str());
        }
    }

    switch (arg.getFunctionSpace().getTypeCode()) {
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
        case Nodes:
        case ReducedNodes:
            break;
        default: {
            throw RipleyException("setToGradient: only supported for nodal input data");
        }
    }

    if (getMPISize()>1) {
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

void RipleyDomain::setToIntegrals(vector<double>& integrals, const escript::Data& arg) const
{
    const RipleyDomain& argDomain=dynamic_cast<const RipleyDomain&>(
            *(arg.getFunctionSpace().getDomain()));
    if (argDomain != *this)
        throw RipleyException("setToIntegrals: illegal domain of integration kernel");

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
            throw RipleyException(msg.str());
        }
    }

}

bool RipleyDomain::isCellOriented(int fsType) const
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
    throw RipleyException(msg.str());
}

bool RipleyDomain::canTag(int fsType) const
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
    throw RipleyException(msg.str());
}

void RipleyDomain::setTags(const int fsType, const int newTag, const escript::Data& mask) const
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
            throw RipleyException(msg.str());
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
            msg << "getNumberOfTagsInUse: invalid function space type "
                << fsType;
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
            msg << "borrowListOfTagsInUse: invalid function space type "
                << fsType;
            throw RipleyException(msg.str());
        }
    }
}

void RipleyDomain::Print_Mesh_Info(const bool full) const
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
    if (row_domain != *this)
        throw RipleyException("newSystemMatrix: domain of row function space does not match the domain of matrix generator");
    const RipleyDomain& col_domain=dynamic_cast<const RipleyDomain&>(*(column_functionspace.getDomain()));
    if (col_domain != *this)
        throw RipleyException("newSystemMatrix: domain of column function space does not match the domain of matrix generator");
    // is the function space type right?
    if (row_functionspace.getTypeCode()==ReducedDegreesOfFreedom)
        reduceRowOrder=true;
    else if (row_functionspace.getTypeCode()!=DegreesOfFreedom)
        throw RipleyException("newSystemMatrix: illegal function space type for system matrix rows");
    if (column_functionspace.getTypeCode()==ReducedDegreesOfFreedom)
        reduceColOrder=true;
    else if (column_functionspace.getTypeCode()!=DegreesOfFreedom)
        throw RipleyException("newSystemMatrix: illegal function space type for system matrix columns");

    // generate matrix
    Paso_SystemMatrixPattern* pattern=getPattern(reduceRowOrder, reduceColOrder);
    Paso_SystemMatrix* matrix = Paso_SystemMatrix_alloc(type, pattern,
            row_blocksize, column_blocksize, FALSE);
    paso::checkPasoError();
    escript::ASM_ptr sma(new SystemMatrixAdapter(matrix, row_blocksize,
                row_functionspace, column_blocksize, column_functionspace));
    return sma;
}

void RipleyDomain::addToSystem(
        escript::AbstractSystemMatrix& mat, escript::Data& rhs,
        std::map<std::string, escript::Data> coefs) const
{
    if (isNotEmpty("d_contact", coefs) || isNotEmpty("y_contact", coefs))
        throw RipleyException(
                    "addToSystem: Ripley does not support contact elements");

    paso::SystemMatrixAdapter* sma = 
                    dynamic_cast<paso::SystemMatrixAdapter*>(&mat);
    if (!sma)
        throw RipleyException(
                    "addToSystem: Ripley only accepts Paso system matrices");

    Paso_SystemMatrix* S = sma->getPaso_SystemMatrix();
    assemblePDE(S, rhs, coefs);
    assemblePDEBoundary(S, rhs, coefs);
    assemblePDEDirac(S, rhs, coefs);
}

void RipleyDomain::addToSystemFromPython(escript::AbstractSystemMatrix& mat,
        escript::Data& rhs, boost::python::list data) const
{
    std::map<std::string, escript::Data> mapping;
    tupleListToMap(mapping, data);
    addToSystem(mat, rhs, mapping);
}

void RipleyDomain::addPDEToSystem(
        escript::AbstractSystemMatrix& mat, escript::Data& rhs,
        const escript::Data& A, const escript::Data& B, const escript::Data& C,
        const escript::Data& D, const escript::Data& X, const escript::Data& Y,
        const escript::Data& d, const escript::Data& y,
        const escript::Data& d_contact, const escript::Data& y_contact,
        const escript::Data& d_dirac,const escript::Data& y_dirac) const
{
    if (!d_contact.isEmpty() || !y_contact.isEmpty())
        throw RipleyException("addPDEToSystem: Ripley does not support contact elements");

    paso::SystemMatrixAdapter* sma=dynamic_cast<paso::SystemMatrixAdapter*>(&mat);
    if (!sma)
        throw RipleyException("addPDEToSystem: Ripley only accepts Paso system matrices");

    Paso_SystemMatrix* S = sma->getPaso_SystemMatrix();

    std::map<std::string, escript::Data> coefs;
    coefs["A"] = A; coefs["B"] = B; coefs["C"] = C; coefs["D"] = D;
    coefs["X"] = X; coefs["Y"] = Y;
    assemblePDE(S, rhs, coefs);

    std::map<std::string, escript::Data> boundary;
    boundary["d"] = d;
    boundary["y"] = y;
    assemblePDEBoundary(S, rhs, boundary);

    map<string, escript::Data> dirac;
    dirac["d_dirac"] = d_dirac;
    dirac["y_dirac"] = y_dirac;
    assemblePDEDirac(S, rhs, dirac);
}

void RipleyDomain::addPDEToRHS(escript::Data& rhs, const escript::Data& X,
        const escript::Data& Y, const escript::Data& y,
        const escript::Data& y_contact, const escript::Data& y_dirac) const
{
    if (!y_contact.isEmpty())
        throw RipleyException("addPDEToRHS: Ripley does not support contact elements");

    if (rhs.isEmpty()) {
        if (!X.isEmpty() || !Y.isEmpty())
            throw RipleyException("addPDEToRHS: right hand side coefficients are provided but no right hand side vector given");
        else
            return;
    }
    {
        std::map<std::string, escript::Data> coefs;
        coefs["X"] = X; coefs["Y"] = Y;
        assemblePDE(NULL, rhs, coefs);
    }
    std::map<std::string, escript::Data> coefs;
    coefs["y"] = y;
    assemblePDEBoundary(NULL, rhs, coefs);
    map<string, escript::Data> dirac;
    dirac["y_dirac"] = y_dirac;
    assemblePDEDirac(NULL, rhs, dirac);
}

void RipleyDomain::setAssemblerFromPython(std::string type,
                                         boost::python::list options) {
    std::map<std::string, escript::Data> mapping;
    tupleListToMap(mapping, options);
    setAssembler(type, mapping);
}

void RipleyDomain::addToRHSFromPython(escript::Data& rhs,
                                         boost::python::list data) const
{
    std::map<std::string, escript::Data> mapping;
    tupleListToMap(mapping, data);
    addToRHS(rhs, mapping);
}

void RipleyDomain::addToRHS(escript::Data& rhs,
        std::map<std::string, escript::Data> coefs) const
{
    if (isNotEmpty("y_contact", coefs))
        throw RipleyException(
                    "addPDEToRHS: Ripley does not support contact elements");

    if (rhs.isEmpty()) {
        if (isNotEmpty("X", coefs) || isNotEmpty("Y", coefs))
            throw RipleyException(
                    "addPDEToRHS: right hand side coefficients are provided "
                    "but no right hand side vector given");
        else
            return;
    }

    assemblePDE(NULL, rhs, coefs);
    assemblePDEBoundary(NULL, rhs, coefs);
    assemblePDEDirac(NULL, rhs, coefs);
}

escript::ATP_ptr RipleyDomain::newTransportProblem(const int blocksize,
        const escript::FunctionSpace& functionspace, const int type) const
{
    bool reduceOrder=false;
    // is the domain right?
    const RipleyDomain& domain=dynamic_cast<const RipleyDomain&>(*(functionspace.getDomain()));
    if (domain != *this)
        throw RipleyException("newTransportProblem: domain of function space does not match the domain of transport problem generator");
    // is the function space type right?
    if (functionspace.getTypeCode()==ReducedDegreesOfFreedom)
        reduceOrder=true;
    else if (functionspace.getTypeCode()!=DegreesOfFreedom)
        throw RipleyException("newTransportProblem: illegal function space type for transport problem");

    // generate matrix
    Paso_SystemMatrixPattern* pattern=getPattern(reduceOrder, reduceOrder);
    Paso_TransportProblem* tp = Paso_TransportProblem_alloc(pattern, blocksize);
    paso::checkPasoError();
    escript::ATP_ptr atp(new TransportProblemAdapter(tp, blocksize, functionspace));
    return atp;
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
    if (!d_contact.isEmpty() || !y_contact.isEmpty())
        throw RipleyException("addPDEToTransportProblem: Ripley does not support contact elements");

    paso::TransportProblemAdapter* tpa=dynamic_cast<paso::TransportProblemAdapter*>(&tp);
    if (!tpa)
        throw RipleyException("addPDEToTransportProblem: Ripley only accepts Paso transport problems");

    Paso_TransportProblem* ptp = tpa->getPaso_TransportProblem();
    std::map<std::string, escript::Data> coefs;
    coefs["D"] = M;
    assemblePDE(ptp->mass_matrix, source, coefs);
    coefs["A"] = A; coefs["B"] = B; coefs["C"] = C; coefs["D"] = D;
    coefs["X"] = X; coefs["Y"] = Y; coefs["d"] = d; coefs["y"] = y;
    coefs["d_dirac"] = d_dirac; coefs["y_dirac"] = y_dirac;
    assemblePDE(ptp->transport_matrix, source, coefs);
    assemblePDEBoundary(ptp->transport_matrix, source, coefs);
    assemblePDEDirac(ptp->transport_matrix, source, coefs);
}

void RipleyDomain::setNewX(const escript::Data& arg)
{
    throw RipleyException("setNewX(): operation not supported");
}

//protected
void RipleyDomain::copyData(escript::Data& out, const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    out.requireWrite();
#pragma omp parallel for
    for (index_t i=0; i<in.getNumSamples(); i++) {
        const double* src = in.getSampleDataRO(i);
        copy(src, src+numComp, out.getSampleDataRW(i));
    }
}

//protected
void RipleyDomain::averageData(escript::Data& out, const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    const dim_t dpp = in.getNumDataPointsPerSample();
    out.requireWrite();
#pragma omp parallel for
    for (index_t i=0; i<in.getNumSamples(); i++) {
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
void RipleyDomain::multiplyData(escript::Data& out, const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    const dim_t dpp = out.getNumDataPointsPerSample();
    out.requireWrite();
#pragma omp parallel for
    for (index_t i=0; i<in.getNumSamples(); i++) {
        const double* src = in.getSampleDataRO(i);
        double* dest = out.getSampleDataRW(i);
        for (index_t c=0; c<numComp; c++) {
            for (index_t q=0; q<dpp; q++)
                dest[c+q*numComp] = src[c];
        }
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
        case Points:
            throw RipleyException("updateTagsInUse for Ripley dirac points "
                    "not supported");
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
	    long i;	// should be size_t but omp mutter mutter
#pragma omp for schedule(static) private(i) nowait
            for (i = 0; i < tags->size(); i++) {
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
Paso_Pattern* RipleyDomain::createPasoPattern(const IndexVector& ptr,
        const IndexVector& index, const dim_t M, const dim_t N) const
{
    // paso will manage the memory
    index_t* indexC = new  index_t[index.size()];
    index_t* ptrC = new  index_t[ptr.size()];
    copy(index.begin(), index.end(), indexC);
    copy(ptr.begin(), ptr.end(), ptrC);
    return Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT, M, N, ptrC, indexC);
}

//protected
Paso_Pattern* RipleyDomain::createMainPattern() const
{
    IndexVector ptr(1,0);
    IndexVector index;

    for (index_t i=0; i<getNumDOF(); i++) {
        // add the DOF itself
        index.push_back(i);
        const dim_t num=insertNeighbourNodes(index, i);
        ptr.push_back(ptr.back()+num+1);
    }

    return createPasoPattern(ptr, index, ptr.size()-1, ptr.size()-1);
}

//protected
void RipleyDomain::createCouplePatterns(const vector<IndexVector>& colIndices,
        const vector<IndexVector>& rowIndices,
        const dim_t N, Paso_Pattern** colPattern, Paso_Pattern** rowPattern) const
{
    IndexVector ptr(1,0);
    IndexVector index;
    for (index_t i=0; i<getNumDOF(); i++) {
        index.insert(index.end(), colIndices[i].begin(), colIndices[i].end());
        ptr.push_back(ptr.back()+colIndices[i].size());
    }

    const dim_t M=ptr.size()-1;
    *colPattern=createPasoPattern(ptr, index, M, N);

    IndexVector rowPtr(1,0);
    IndexVector rowIndex;
    for (index_t i=0; i<rowIndices.size(); i++) {
        rowIndex.insert(rowIndex.end(), rowIndices[i].begin(), rowIndices[i].end());
        rowPtr.push_back(rowPtr.back()+rowIndices[i].size());
    }

    // M and N are now reversed
    *rowPattern=createPasoPattern(rowPtr, rowIndex, N, M);
}

//protected
void RipleyDomain::addToSystemMatrix(Paso_SystemMatrix* mat, 
       const IndexVector& nodes_Eq, dim_t num_Eq, const IndexVector& nodes_Sol,
       dim_t num_Sol, const vector<double>& array) const
{
    if (mat->type & MATRIX_FORMAT_TRILINOS_CRS)
        throw RipleyException("addToSystemMatrix: TRILINOS_CRS not supported");

    const dim_t numMyCols = mat->pattern->mainPattern->numInput;
    const dim_t numMyRows = mat->pattern->mainPattern->numOutput;
    const dim_t numSubblocks_Eq = num_Eq / mat->row_block_size;
    const dim_t numSubblocks_Sol = num_Sol / mat->col_block_size;

    const index_t* mainBlock_ptr = mat->mainBlock->pattern->ptr;
    const index_t* mainBlock_index = mat->mainBlock->pattern->index;
    double* mainBlock_val = mat->mainBlock->val;
    const index_t* col_coupleBlock_ptr = mat->col_coupleBlock->pattern->ptr;
    const index_t* col_coupleBlock_index = mat->col_coupleBlock->pattern->index;
    double* col_coupleBlock_val = mat->col_coupleBlock->val;
    const index_t* row_coupleBlock_ptr = mat->row_coupleBlock->pattern->ptr;
    const index_t* row_coupleBlock_index = mat->row_coupleBlock->pattern->index;
    double* row_coupleBlock_val = mat->row_coupleBlock->val;
    index_t offset=(mat->type & MATRIX_FORMAT_OFFSET1 ? 1:0);

#define UPDATE_BLOCK(VAL) do {\
    for (dim_t ic=0; ic<mat->col_block_size; ++ic) {\
        const dim_t i_Sol=ic+mat->col_block_size*l_col;\
        for (dim_t ir=0; ir<mat->row_block_size; ++ir) {\
            const dim_t i_Eq=ir+mat->row_block_size*l_row;\
            VAL[k*mat->block_size+ir+mat->row_block_size*ic]\
                += array[INDEX4(i_Eq, i_Sol, k_Eq, k_Sol, num_Eq, num_Sol, nodes_Eq.size())];\
        }\
    }\
} while(0)

    if (mat->type & MATRIX_FORMAT_CSC) {
        for (dim_t k_Sol = 0; k_Sol < nodes_Sol.size(); ++k_Sol) {
            // down columns of array
            for (dim_t l_col = 0; l_col < numSubblocks_Sol; ++l_col) {
                const dim_t i_col = nodes_Sol[k_Sol]*numSubblocks_Sol+l_col;
                if (i_col < numMyCols) {
                    for (dim_t k_Eq = 0; k_Eq < nodes_Eq.size(); ++k_Eq) {
                        for (dim_t l_row = 0; l_row < numSubblocks_Eq; ++l_row) {
                            const dim_t i_row = nodes_Eq[k_Eq]*numSubblocks_Eq+l_row+offset;
                            if (i_row < numMyRows+offset) {
                                for (dim_t k = mainBlock_ptr[i_col]-offset; k < mainBlock_ptr[i_col+1]-offset; ++k) {
                                    if (mainBlock_index[k] == i_row) {
                                        UPDATE_BLOCK(mainBlock_val);
                                        break;
                                    }
                                }
                            } else {
                                for (dim_t k = col_coupleBlock_ptr[i_col]-offset; k < col_coupleBlock_ptr[i_col+1]-offset; ++k) {
                                    if (row_coupleBlock_index[k] == i_row - numMyRows) {
                                        UPDATE_BLOCK(row_coupleBlock_val);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    for (dim_t k_Eq = 0; k_Eq < nodes_Eq.size(); ++k_Eq) {
                        // across rows of array
                        for (dim_t l_row=0; l_row<numSubblocks_Eq; ++l_row) {
                            const dim_t i_row = nodes_Eq[k_Eq]*numSubblocks_Eq+l_row+offset;
                            if (i_row < numMyRows+offset) {
                                for (dim_t k = col_coupleBlock_ptr[i_col-numMyCols]-offset;
                                     k < col_coupleBlock_ptr[i_col-numMyCols+1]-offset; ++k)
                                {
                                    if (col_coupleBlock_index[k] == i_row) {
                                        UPDATE_BLOCK(col_coupleBlock_val);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        for (dim_t k_Eq = 0; k_Eq < nodes_Eq.size(); ++k_Eq) {
            // down columns of array
            for (dim_t l_row = 0; l_row < numSubblocks_Eq; ++l_row) {
                const dim_t i_row = nodes_Eq[k_Eq]*numSubblocks_Eq+l_row;
                // only look at the matrix rows stored on this processor
                if (i_row < numMyRows) {
                    for (dim_t k_Sol = 0; k_Sol < nodes_Sol.size(); ++k_Sol) {
                        for (dim_t l_col = 0; l_col < numSubblocks_Sol; ++l_col) {
                            const dim_t i_col = nodes_Sol[k_Sol]*numSubblocks_Sol+l_col+offset;
                            if (i_col < numMyCols+offset) {
                                for (dim_t k = mainBlock_ptr[i_row]-offset; k < mainBlock_ptr[i_row+1]-offset; ++k) {
                                    if (mainBlock_index[k] == i_col) {
                                        UPDATE_BLOCK(mainBlock_val);
                                        break;
                                    }
                                }
                            } else {
                                for (dim_t k = col_coupleBlock_ptr[i_row]-offset; k < col_coupleBlock_ptr[i_row+1]-offset; ++k) {
                                    if (col_coupleBlock_index[k] == i_col-numMyCols) {
                                        UPDATE_BLOCK(col_coupleBlock_val);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    for (dim_t k_Sol = 0; k_Sol < nodes_Sol.size(); ++k_Sol) {
                        // across rows of array
                        for (dim_t l_col=0; l_col<numSubblocks_Sol; ++l_col) {
                            const dim_t i_col = nodes_Sol[k_Sol]*numSubblocks_Sol+l_col+offset;
                            if (i_col < numMyCols+offset) {
                                for (dim_t k = row_coupleBlock_ptr[i_row-numMyRows]-offset;
                                     k < row_coupleBlock_ptr[i_row-numMyRows+1]-offset; ++k)
                                {
                                    if (row_coupleBlock_index[k] == i_col) {
                                        UPDATE_BLOCK(row_coupleBlock_val);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#undef UPDATE_BLOCK
}

//private
void RipleyDomain::assemblePDE(Paso_SystemMatrix* mat, escript::Data& rhs,
        std::map<std::string, escript::Data> coefs) const
{
    if (rhs.isEmpty() && isNotEmpty("X", coefs) && isNotEmpty("Y", coefs))
        throw RipleyException("assemblePDE: right hand side coefficients are "
                    "provided but no right hand side vector given");

    vector<int> fsTypes;
    if (isNotEmpty("A", coefs))
        fsTypes.push_back(coefs["A"].getFunctionSpace().getTypeCode());
    if (isNotEmpty("B", coefs))
        fsTypes.push_back(coefs["B"].getFunctionSpace().getTypeCode());
    if (isNotEmpty("C", coefs))
        fsTypes.push_back(coefs["C"].getFunctionSpace().getTypeCode());
    if (isNotEmpty("D", coefs))
        fsTypes.push_back(coefs["D"].getFunctionSpace().getTypeCode());
    if (isNotEmpty("X", coefs)) //may not exist for some custom assemblers
        fsTypes.push_back(coefs["X"].getFunctionSpace().getTypeCode());
    else if (isNotEmpty("du", coefs)) //used for (some) custom assemblers
        fsTypes.push_back(coefs["du"].getFunctionSpace().getTypeCode());
    if (isNotEmpty("Y", coefs))
        fsTypes.push_back(coefs["Y"].getFunctionSpace().getTypeCode());
    
    if (fsTypes.empty()) {
        return;
    }

    int fs=fsTypes[0];
    if (fs != Elements && fs != ReducedElements)
        throw RipleyException("assemblePDE: illegal function space type for coefficients");

    for (vector<int>::const_iterator it=fsTypes.begin()+1; it!=fsTypes.end(); it++) {
        if (*it != fs) {
            throw RipleyException("assemblePDE: coefficient function spaces don't match");
        }
    }

    int numEq, numComp;
    if (!mat) {
        if (rhs.isEmpty()) {
            numEq=numComp=1;
        } else {
            numEq=numComp=rhs.getDataPointSize();
        }
    } else {
        if (!rhs.isEmpty() && rhs.getDataPointSize()!=mat->logical_row_block_size)
            throw RipleyException("assemblePDE: matrix row block size and number of components of right hand side don't match");
        numEq = mat->logical_row_block_size;
        numComp = mat->logical_col_block_size;
    }

    if (numEq != numComp)
        throw RipleyException("assemblePDE: number of equations and number of solutions don't match");

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
void RipleyDomain::assemblePDEBoundary(Paso_SystemMatrix* mat,
      escript::Data& rhs, std::map<std::string, escript::Data> coefs) const
{
    std::map<std::string, escript::Data>::iterator iy = coefs.find("y"),
                                                   id = coefs.find("d");
    if (rhs.isEmpty() && isNotEmpty("y", coefs))
        throw RipleyException("assemblePDEBoundary: y provided but no right hand side vector given");

    int fs=-1;
    if (isNotEmpty("d", coefs))
        fs=id->second.getFunctionSpace().getTypeCode();
    if (isNotEmpty("y", coefs)) {
        if (fs == -1)
            fs = iy->second.getFunctionSpace().getTypeCode();
        else if (fs != iy->second.getFunctionSpace().getTypeCode())
            throw RipleyException("assemblePDEBoundary: coefficient function spaces don't match");
    }
    if (fs==-1) {
        return;
    }
    
    if (fs != FaceElements && fs != ReducedFaceElements)
        throw RipleyException("assemblePDEBoundary: illegal function space type for coefficients");

    int numEq, numComp;
    if (!mat) {
        if (rhs.isEmpty()) {
            numEq=numComp=1;
        } else {
            numEq=numComp=rhs.getDataPointSize();
        }
    } else {
        if (!rhs.isEmpty() && rhs.getDataPointSize()!=mat->logical_row_block_size)
            throw RipleyException("assemblePDEBoundary: matrix row block size and number of components of right hand side don't match");
        numEq = mat->logical_row_block_size;
        numComp = mat->logical_col_block_size;
    }

    if (numEq != numComp)
        throw RipleyException("assemblePDEBoundary: number of equations and number of solutions don't match");

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

void RipleyDomain::assemblePDEDirac(Paso_SystemMatrix* mat,
        escript::Data& rhs, std::map<std::string, escript::Data> coefs) const
{
    bool yNotEmpty = isNotEmpty("y_dirac", coefs),
         dNotEmpty = isNotEmpty("d_dirac", coefs);
    escript::Data d = dNotEmpty ? coefs["d_dirac"] : escript::Data(),
                  y = yNotEmpty ? coefs["y_dirac"] : escript::Data();
    if (!(yNotEmpty || dNotEmpty)) {
        return;
    }
    int nEq, nComp;
    if (!mat) {
        if (rhs.isEmpty()) {
            nEq=nComp=1;
        } else {
            nEq=nComp=rhs.getDataPointSize();
        }
    } else {
        if (!rhs.isEmpty() && rhs.getDataPointSize()!=mat->logical_row_block_size)
            throw RipleyException("assemblePDEDirac: matrix row block size "
                    "and number of components of right hand side don't match");
        nEq = mat->logical_row_block_size;
        nComp = mat->logical_col_block_size;
    }
    for (int i = 0; i < m_diracPoints.size(); i++) { //only for this rank
        IndexVector rowIndex;
        rowIndex.push_back(getDofOfNode(m_diracPoints[i].node));
        if (yNotEmpty) {
            const double *EM_F = y.getSampleDataRO(i);
            double *F_p = rhs.getSampleDataRW(0);
            for (index_t eq = 0; eq < nEq; eq++) {
                F_p[INDEX2(eq, rowIndex[0], nEq)] += EM_F[INDEX2(eq,i,nEq)];
            }
        }
        if (dNotEmpty) {
            const double *EM_S = d.getSampleDataRO(i);
            std::vector<double> contents(EM_S,
                        EM_S+mat->row_block_size*nEq*nComp*rowIndex.size());
            addToSystemMatrix(mat, rowIndex, nEq, rowIndex, nComp, contents);
        }
    }
}

bool RipleyDomain::probeInterpolationACross(int fsType_source,
        const escript::AbstractDomain&, int fsType_target) const
{
    //TODO
    return false;
}

void RipleyDomain::interpolateACross(escript::Data& target, const escript::Data& source) const
{
    throw RipleyException("interpolateACross() not supported");
}

// Expecting ("gaussian", radius, sigma)
bool RipleyDomain::supportsFilter(const boost::python::tuple& t) const
{
    if (len(t)==0) {	// so we can handle unfiltered randoms
        return true;
    }
    if (len(t)!=3) {
        return false;
    }
    boost::python::extract<string> ex(t[0]);
    if (!ex.check() || (ex()!="gaussian"))
    {
        return false;
    }
    if (! boost::python::extract<unsigned int>(t[1]).check())
    {
        return false;
    }
    return boost::python::extract<double>(t[2]).check();
}

void RipleyDomain::addPoints(int numPoints, const double* points_ptr,
                     const int* tags_ptr)
{
    for (int i = 0; i < numPoints; i++) {
        int node = findNode(&points_ptr[i * m_numDim]);
        if (node >= 0) {
            m_diracPointNodeIDs.push_back(borrowSampleReferenceIDs(Nodes)[node]);
            DiracPoint dp;
            dp.node = node; //local
            dp.tag = tags_ptr[i];
            m_diracPoints.push_back(dp);
        }
    }
}

} // end of namespace ripley


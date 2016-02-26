
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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
#include <ripley/domainhelpers.h>

#include <escript/DataFactory.h>
#include <escript/FunctionSpaceFactory.h>
#include <escript/SolverOptions.h>

#ifdef USE_CUDA
#include <ripley/RipleySystemMatrix.h>
#endif

#include <paso/SystemMatrix.h>
#include <paso/Transport.h>

#include <iomanip>

namespace bp = boost::python;

using namespace std;
using escript::ValueError;
using escript::NotImplementedError;
using paso::TransportProblem;

namespace ripley {

void tupleListToMap(DataMap& mapping, const bp::list& list)
{
    for (int i = 0; i < len(list); i++) {
        if (!bp::extract<bp::tuple>(list[i]).check())
            throw ValueError("Passed in list contains objects"
                                      " other than tuples");
        bp::tuple t = bp::extract<bp::tuple>(list[i]);
        if (len(t) != 2 || !bp::extract<string>(t[0]).check() ||
                !bp::extract<escript::Data>(t[1]).check())
            throw ValueError("The passed in list must contain tuples"
                " of the form (string, escript::Data)");
        mapping[bp::extract<string>(t[0])] = bp::extract<escript::Data>(t[1]);
    }
}

RipleyDomain::RipleyDomain(dim_t dim, escript::SubWorld_ptr p) :
    m_numDim(dim),
    m_status(0)
{
    if (p.get() == NULL)
        m_mpiInfo = esysUtils::makeInfo(MPI_COMM_WORLD);
    else
        m_mpiInfo = p->getMPI();

    assembler_type = DEFAULT_ASSEMBLER;
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
        case DegreesOfFreedom:
            return "Ripley_DegreesOfFreedom [Solution(domain)]";
        case ReducedDegreesOfFreedom:
            return "Ripley_ReducedDegreesOfFreedom [ReducedSolution(domain)]";
        case Nodes:
            return "Ripley_Nodes [ContinuousFunction(domain)]";
        case ReducedNodes:
            return "Ripley_ReducedNodes [ReducedContinuousFunction(domain)]";
        case Elements:
            return "Ripley_Elements [Function(domain)]";
        case ReducedElements:
            return "Ripley_ReducedElements [ReducedFunction(domain)]";
        case FaceElements:
            return "Ripley_FaceElements [FunctionOnBoundary(domain)]";
        case ReducedFaceElements:
            return "Ripley_ReducedFaceElements [ReducedFunctionOnBoundary(domain)]";
        case Points:
            return "Ripley_Points [DiracDeltaFunctions(domain)]";
        default:
            break;
    }
    return "Invalid function space type code";
}

pair<int,dim_t> RipleyDomain::getDataShape(int fsType) const
{
    const int ptsPerSample = (m_numDim==2 ? 4 : 8);
    switch (fsType) {
        case Nodes:
        case ReducedNodes:
            return pair<int,dim_t>(1, getNumNodes());
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
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
    throw ValueError(msg.str());
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
        throw ValueError(msg.str());
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
            throw ValueError(msg.str());
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
        throw ValueError(msg.str());
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
            throw ValueError(msg.str());
        }
    }
}

void RipleyDomain::interpolateOnDomain(escript::Data& target,
                                       const escript::Data& in) const
{
    const RipleyDomain& inDomain=dynamic_cast<const RipleyDomain&>(*(in.getFunctionSpace().getDomain()));
    const RipleyDomain& targetDomain=dynamic_cast<const RipleyDomain&>(*(target.getFunctionSpace().getDomain()));
    if (inDomain != *this)
        throw ValueError("Illegal domain of interpolant");
    if (targetDomain != *this)
        throw ValueError("Illegal domain of interpolation target");

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
        throw ValueError("interpolateOnDomain: Cannot interpolate "
                              "reduced data to non-reduced data.");
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
            case ReducedNodes:
                switch (outFS) {
                    case DegreesOfFreedom:
                    case ReducedDegreesOfFreedom:
                        if (getMPISize()==1)
                            copyData(target, in);
                        else
                            nodesToDOF(target, in);
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
                        throw NotImplementedError(msg.str());
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
                            dofToNodes(target, in);
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
                        throw NotImplementedError(msg.str());
                }
                break;
            case Points:
                switch(outFS) {
                    case Nodes:
                        {
                            const dim_t numComp = in.getDataPointSize();
                            const int nDirac = m_diracPoints.size();
                            target.requireWrite();
#pragma omp parallel for
                            for (int i = 0; i < nDirac; i++) {
                                const double* src = in.getSampleDataRO(i);
                                copy(src, src+numComp, target.getSampleDataRW(m_diracPoints[i].node));
                            }
                        }
                
                }
                break;
            default:
                throw NotImplementedError(msg.str());
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
        throw ValueError("setToX: Illegal domain of data point locations");
    if (!arg.isExpanded())
        throw ValueError("setToX: Expanded Data object expected");

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
        throw ValueError("setToGradient: Illegal domain of gradient argument");
    const RipleyDomain& gradDomain=dynamic_cast<const RipleyDomain&>(
            *(grad.getFunctionSpace().getDomain()));
    if (gradDomain != *this)
        throw ValueError("setToGradient: Illegal domain of gradient");

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
            throw ValueError(msg.str());
        }
    }

    switch (arg.getFunctionSpace().getTypeCode()) {
        case DegreesOfFreedom:
        case ReducedDegreesOfFreedom:
        case Nodes:
        case ReducedNodes:
            break;
        default: {
            throw ValueError("setToGradient: only supported for nodal input data");
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

void RipleyDomain::setToIntegrals(vector<double>& integrals, const escript::Data& arg) const
{
    const RipleyDomain& argDomain=dynamic_cast<const RipleyDomain&>(
            *(arg.getFunctionSpace().getDomain()));
    if (argDomain != *this)
        throw ValueError("setToIntegrals: illegal domain of integration kernel");

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
            throw ValueError(msg.str());
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
    throw ValueError(msg.str());
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
    throw ValueError(msg.str());
}

void RipleyDomain::setTags(int fsType, int newTag, const escript::Data& mask) const
{
    vector<int>* target=NULL;
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
            throw ValueError(msg.str());
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

int RipleyDomain::getTagFromSampleNo(int fsType, dim_t sampleNo) const
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
            throw ValueError(msg.str());
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
            throw ValueError(msg.str());
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
            throw ValueError(msg.str());
        }
    }
}

void RipleyDomain::Print_Mesh_Info(bool full) const
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

int RipleyDomain::getSystemMatrixTypeId(const bp::object& options) const
{
    const escript::SolverBuddy& sb = bp::extract<escript::SolverBuddy>(options);
    int package = sb.getPackage();

    // use CUSP for single rank and supported solvers+preconditioners if CUDA
    // is available, PASO otherwise
    if (package == escript::SO_DEFAULT) {
#ifdef USE_CUDA
        if (m_mpiInfo->size == 1) {
            switch (sb.getSolverMethod()) {
                case escript::SO_DEFAULT:
                case escript::SO_METHOD_BICGSTAB:
                case escript::SO_METHOD_CGLS:
                case escript::SO_METHOD_GMRES:
                case escript::SO_METHOD_LSQR:
                case escript::SO_METHOD_PCG:
                case escript::SO_METHOD_PRES20:
                    package = escript::SO_PACKAGE_CUSP;
                    break;
                default:
                    package = escript::SO_PACKAGE_PASO;
            }
            if (package == escript::SO_PACKAGE_CUSP) {
                if (sb.getPreconditioner() != escript::SO_PRECONDITIONER_NONE &&
                        sb.getPreconditioner() != escript::SO_PRECONDITIONER_JACOBI) {
                    package = escript::SO_PACKAGE_PASO;
                }
            }
        } else {
            package = escript::SO_PACKAGE_PASO;
        }
#else // USE_CUDA
        package = escript::SO_PACKAGE_PASO;
#endif
    }

    if (package == escript::SO_PACKAGE_CUSP) {
        if (m_mpiInfo->size > 1) {
            throw NotImplementedError("CUSP matrices are not supported with more than one rank");
        }
        int type = (int)SMT_CUSP;
        if (sb.isSymmetric())
            type |= (int)SMT_SYMMETRIC;
        return type;
    }

    // in all other cases we use PASO
    return (int)SMT_PASO | paso::SystemMatrix::getSystemMatrixTypeId(
            sb.getSolverMethod(), sb.getPreconditioner(), sb.getPackage(),
            sb.isSymmetric(), m_mpiInfo);
}

int RipleyDomain::getTransportTypeId(int solver, int preconditioner,
                                     int package, bool symmetry) const
{
    return TransportProblem::getTypeId(solver, preconditioner,
            package, symmetry, m_mpiInfo);
}

escript::ASM_ptr RipleyDomain::newSystemMatrix(int row_blocksize,
        const escript::FunctionSpace& row_functionspace, int column_blocksize,
        const escript::FunctionSpace& column_functionspace, int type) const
{
    bool reduceRowOrder=false;
    bool reduceColOrder=false;
    // is the domain right?
    const RipleyDomain& row_domain=dynamic_cast<const RipleyDomain&>(*(row_functionspace.getDomain()));
    if (row_domain != *this)
        throw ValueError("newSystemMatrix: domain of row function space does not match the domain of matrix generator");
    const RipleyDomain& col_domain=dynamic_cast<const RipleyDomain&>(*(column_functionspace.getDomain()));
    if (col_domain != *this)
        throw ValueError("newSystemMatrix: domain of column function space does not match the domain of matrix generator");
    // is the function space type right?
    if (row_functionspace.getTypeCode()==ReducedDegreesOfFreedom)
        reduceRowOrder=true;
    else if (row_functionspace.getTypeCode()!=DegreesOfFreedom)
        throw ValueError("newSystemMatrix: illegal function space type for system matrix rows");
    if (column_functionspace.getTypeCode()==ReducedDegreesOfFreedom)
        reduceColOrder=true;
    else if (column_functionspace.getTypeCode()!=DegreesOfFreedom)
        throw ValueError("newSystemMatrix: illegal function space type for system matrix columns");
    // are block sizes identical?
    if (row_blocksize != column_blocksize)
        throw ValueError("newSystemMatrix: row/column block sizes must be equal");
    // are function spaces equal
    if (reduceRowOrder != reduceColOrder)
        throw ValueError("newSystemMatrix: row/column function spaces must be equal");

    // generate matrix
    //if (reduceRowOrder || reduceColOrder)
    //    throw NotImplementedError("newSystemMatrix: reduced order not supported");

    if (type & (int)SMT_CUSP) {
#ifdef USE_CUDA
        const dim_t numMatrixRows = getNumDOF();
        bool symmetric = (type & (int)SMT_SYMMETRIC);
        escript::ASM_ptr sm(new SystemMatrix(m_mpiInfo, row_blocksize,
                    row_functionspace, numMatrixRows,
                    getDiagonalIndices(symmetric), symmetric));
        return sm;
#else
        throw RipleyException("newSystemMatrix: ripley was compiled without CUDA support so CUSP solvers & matrices are not available.");
#endif
    } else if (type & (int)SMT_PASO) {
        paso::SystemMatrixPattern_ptr pattern(getPasoMatrixPattern(
                                            reduceRowOrder, reduceColOrder));
        type -= (int)SMT_PASO;
        escript::ASM_ptr sm(new paso::SystemMatrix(type, pattern,
                row_blocksize, column_blocksize, false, row_functionspace,
                column_functionspace));
        return sm;
    } else {
        throw RipleyException("newSystemMatrix: unknown matrix type ID");
    }
}

void RipleyDomain::addToSystem(escript::AbstractSystemMatrix& mat,
                               escript::Data& rhs, const DataMap& coefs,
                               Assembler_ptr assembler) const
{
    if (isNotEmpty("d_contact", coefs) || isNotEmpty("y_contact", coefs))
        throw ValueError(
                    "addToSystem: Ripley does not support contact elements");

    assemblePDE(&mat, rhs, coefs, assembler);
    assemblePDEBoundary(&mat, rhs, coefs, assembler);
    assemblePDEDirac(&mat, rhs, coefs, assembler);
}

void RipleyDomain::addToSystemFromPython(escript::AbstractSystemMatrix& mat,
                                         escript::Data& rhs,
                                         const bp::list& data,
                                         Assembler_ptr assembler) const
{
    DataMap mapping;
    tupleListToMap(mapping, data);
    addToSystem(mat, rhs, mapping, assembler);
}

Assembler_ptr RipleyDomain::createAssemblerFromPython(const string type,
                                                const bp::list& options) const
{
    DataMap mapping;
    tupleListToMap(mapping, options);
    return createAssembler(type, mapping);
}

void RipleyDomain::addToRHSFromPython(escript::Data& rhs, const bp::list& data,
                                      Assembler_ptr assembler) const
{
    DataMap mapping;
    tupleListToMap(mapping, data);
    addToRHS(rhs, mapping, assembler);
}

void RipleyDomain::addToRHS(escript::Data& rhs, const DataMap& coefs,
                            Assembler_ptr assembler) const
{
    if (isNotEmpty("y_contact", coefs))
        throw ValueError(
                    "addPDEToRHS: Ripley does not support contact elements");

    if (rhs.isEmpty()) {
        if ((isNotEmpty("X", coefs) && isNotEmpty("du", coefs))
                || isNotEmpty("Y", coefs))
            throw ValueError(
                    "addPDEToRHS: right hand side coefficients are provided "
                    "but no right hand side vector given");
        else
            return;
    }

    assemblePDE(NULL, rhs, coefs, assembler);
    assemblePDEBoundary(NULL, rhs, coefs, assembler);
    assemblePDEDirac(NULL, rhs, coefs, assembler);
}

escript::ATP_ptr RipleyDomain::newTransportProblem(int blocksize,
                  const escript::FunctionSpace& functionspace, int type) const
{
    bool reduceOrder=false;
    // is the domain right?
    const RipleyDomain& domain=dynamic_cast<const RipleyDomain&>(*(functionspace.getDomain()));
    if (domain != *this)
        throw ValueError("newTransportProblem: domain of function space does not match the domain of transport problem generator");
    // is the function space type right?
    if (functionspace.getTypeCode()==ReducedDegreesOfFreedom)
        reduceOrder=true;
    else if (functionspace.getTypeCode()!=DegreesOfFreedom)
        throw ValueError("newTransportProblem: illegal function space type for transport problem");

    // generate matrix
    paso::SystemMatrixPattern_ptr pattern(getPasoMatrixPattern(reduceOrder,
                                                               reduceOrder));
    escript::ATP_ptr tp(new paso::TransportProblem(pattern, blocksize,
                                                   functionspace));
    return tp;
}

void RipleyDomain::addPDEToTransportProblemFromPython(
                escript::AbstractTransportProblem& tp, escript::Data& source,
                const bp::list& data, Assembler_ptr assembler) const
{
    DataMap mapping;
    tupleListToMap(mapping, data);
    addPDEToTransportProblem(tp, source, mapping, assembler);
}

void RipleyDomain::addPDEToTransportProblem(
                escript::AbstractTransportProblem& tp, escript::Data& source,
                const DataMap& coefs, Assembler_ptr assembler) const
{
    if (isNotEmpty("d_contact", coefs) || isNotEmpty("y_contact", coefs))
        throw ValueError("addPDEToTransportProblem: Ripley does not support contact elements");

    TransportProblem* ptp=dynamic_cast<TransportProblem*>(&tp);
    if (!ptp)
        throw ValueError("addPDEToTransportProblem: Ripley only accepts Paso transport problems");

    escript::ASM_ptr mm(ptp->borrowMassMatrix());
    escript::ASM_ptr tm(ptp->borrowTransportMatrix());

    assemblePDE(mm.get(), source, coefs, assembler);
    assemblePDE(tm.get(), source, coefs, assembler);
    assemblePDEBoundary(tm.get(), source, coefs, assembler);
    assemblePDEDirac(tm.get(), source, coefs, assembler);
}

void RipleyDomain::addPDEToTransportProblem(
                     escript::AbstractTransportProblem& tp, escript::Data& source,
                     const escript::Data& M,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C,const  escript::Data& D,
                     const  escript::Data& X,const  escript::Data& Y,
                     const escript::Data& d, const escript::Data& y,
                     const escript::Data& d_contact,const escript::Data& y_contact,
                     const escript::Data& d_dirac,const escript::Data& y_dirac) const
{
    throw RipleyException("Programmer error: incorrect version of addPDEToTransportProblem called");

}


void RipleyDomain::setNewX(const escript::Data& arg)
{
    throw NotImplementedError("setNewX(): operation not supported");
}

//protected
void RipleyDomain::copyData(escript::Data& out, const escript::Data& in) const
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
void RipleyDomain::averageData(escript::Data& out, const escript::Data& in) const
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
void RipleyDomain::multiplyData(escript::Data& out, const escript::Data& in) const
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
void RipleyDomain::updateTagsInUse(int fsType) const
{
    vector<int>* tagsInUse=NULL;
    const vector<int>* tags=NULL;
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
            throw NotImplementedError("updateTagsInUse for Ripley dirac points"
                                      " not supported");
        default:
            return;
    }

    // gather global unique tag values from tags into tagsInUse
    tagsInUse->clear();
    int lastFoundValue = numeric_limits<int>::min();
    int minFoundValue, local_minFoundValue;
    const int numTags = tags->size();

    while (true) {
        // find smallest value bigger than lastFoundValue
        minFoundValue = numeric_limits<int>::max();
#pragma omp parallel private(local_minFoundValue)
        {
            local_minFoundValue = minFoundValue;
#pragma omp for schedule(static) nowait
            for (int i = 0; i < numTags; i++) {
                const int v = (*tags)[i];
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
        if (minFoundValue < numeric_limits<int>::max()) {
            tagsInUse->push_back(minFoundValue);
            lastFoundValue = minFoundValue;
        } else
            break;
    }
}

//protected
paso::Pattern_ptr RipleyDomain::createPasoPattern(
                            const vector<IndexVector>& indices, dim_t N) const
{
    // paso will manage the memory
    const dim_t M = indices.size();
    index_t* ptr = new index_t[M+1];
    ptr[0] = 0;
    for (index_t i=0; i < M; i++) {
        ptr[i+1] = ptr[i]+indices[i].size();
    }

    index_t* index = new index_t[ptr[M]];

#pragma omp parallel for
    for (index_t i=0; i < M; i++) {
        copy(indices[i].begin(), indices[i].end(), &index[ptr[i]]);
    }

    return paso::Pattern_ptr(new paso::Pattern(MATRIX_FORMAT_DEFAULT, M, N, ptr, index));
}

//protected
void RipleyDomain::addToSystemMatrix(escript::AbstractSystemMatrix* mat,
                                     const IndexVector& nodes, dim_t numEq,
                                     const DoubleVector& array) const
{
    paso::SystemMatrix* sm = dynamic_cast<paso::SystemMatrix*>(mat);
    if (sm) {
        addToPasoMatrix(sm, nodes, numEq, array);
    } else {
#ifdef USE_CUDA
        SystemMatrix* sm = dynamic_cast<SystemMatrix*>(mat);
        if (sm) {
            sm->add(nodes, array);
        } else {
            throw RipleyException("addToSystemMatrix: unknown system matrix type");
        }
#else
        throw RipleyException("addToSystemMatrix: unknown system matrix type");
#endif
    }
}

//private
void RipleyDomain::addToPasoMatrix(paso::SystemMatrix* mat,
                                   const IndexVector& nodes, dim_t numEq,
                                   const vector<double>& array) const
{
    const dim_t numMyCols = mat->pattern->mainPattern->numInput;
    const dim_t numMyRows = mat->pattern->mainPattern->numOutput;
    const dim_t numSubblocks_Eq = numEq / mat->row_block_size;
    const dim_t numSubblocks_Sol = numEq / mat->col_block_size;

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
                += array[INDEX4(i_Eq, i_Sol, k_Eq, k_Sol, numEq, numEq, nodes.size())];\
        }\
    }\
} while(0)

    if (mat->type & MATRIX_FORMAT_CSC) {
        for (dim_t k_Sol = 0; k_Sol < nodes.size(); ++k_Sol) {
            // down columns of array
            for (dim_t l_col = 0; l_col < numSubblocks_Sol; ++l_col) {
                const dim_t i_col = nodes[k_Sol]*numSubblocks_Sol+l_col;
                if (i_col < numMyCols) {
                    for (dim_t k_Eq = 0; k_Eq < nodes.size(); ++k_Eq) {
                        for (dim_t l_row = 0; l_row < numSubblocks_Eq; ++l_row) {
                            const dim_t i_row = nodes[k_Eq]*numSubblocks_Eq+l_row+offset;
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
                    for (dim_t k_Eq = 0; k_Eq < nodes.size(); ++k_Eq) {
                        // across rows of array
                        for (dim_t l_row=0; l_row<numSubblocks_Eq; ++l_row) {
                            const dim_t i_row = nodes[k_Eq]*numSubblocks_Eq+l_row+offset;
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
        for (dim_t k_Eq = 0; k_Eq < nodes.size(); ++k_Eq) {
            // down columns of array
            for (dim_t l_row = 0; l_row < numSubblocks_Eq; ++l_row) {
                const dim_t i_row = nodes[k_Eq]*numSubblocks_Eq+l_row;
                // only look at the matrix rows stored on this processor
                if (i_row < numMyRows) {
                    for (dim_t k_Sol = 0; k_Sol < nodes.size(); ++k_Sol) {
                        for (dim_t l_col = 0; l_col < numSubblocks_Sol; ++l_col) {
                            const dim_t i_col = nodes[k_Sol]*numSubblocks_Sol+l_col+offset;
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
                    for (dim_t k_Sol = 0; k_Sol < nodes.size(); ++k_Sol) {
                        // across rows of array
                        for (dim_t l_col=0; l_col<numSubblocks_Sol; ++l_col) {
                            const dim_t i_col = nodes[k_Sol]*numSubblocks_Sol+l_col+offset;
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
void RipleyDomain::assemblePDE(escript::AbstractSystemMatrix* mat,
                               escript::Data& rhs, const DataMap& coefs,
                               Assembler_ptr assembler) const
{
    if (rhs.isEmpty() && (isNotEmpty("X", coefs) || isNotEmpty("du", coefs))
                && isNotEmpty("Y", coefs))
        throw ValueError("assemblePDE: right hand side coefficients are "
                         "provided but no right hand side vector given");

    vector<int> fsTypes;
    assembler->collateFunctionSpaceTypes(fsTypes, coefs);

    if (fsTypes.empty()) {
        return;
    }
    
    int fs=fsTypes[0];
    if (fs != Elements && fs != ReducedElements)
        throw ValueError("assemblePDE: illegal function space type for coefficients");

    for (vector<int>::const_iterator it=fsTypes.begin()+1; it!=fsTypes.end(); it++) {
        if (*it != fs) {
            throw ValueError("assemblePDE: coefficient function spaces don't match");
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
        if (!rhs.isEmpty() && rhs.getDataPointSize() != mat->getRowBlockSize())
            throw ValueError("assemblePDE: matrix row block size and number of components of right hand side don't match");
        numEq = mat->getRowBlockSize();
        numComp = mat->getColumnBlockSize();
    }

    if (numEq != numComp)
        throw ValueError("assemblePDE: number of equations and number of solutions don't match");

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
void RipleyDomain::assemblePDEBoundary(escript::AbstractSystemMatrix* mat,
        escript::Data& rhs, const DataMap& coefs, Assembler_ptr assembler) const
{
    if (rhs.isEmpty() && isNotEmpty("y", coefs))
        throw ValueError("assemblePDEBoundary: y provided but no right hand side vector given");

    int fs=-1;
    if (isNotEmpty("d", coefs))
        fs=coefs.find("d")->second.getFunctionSpace().getTypeCode();
    if (isNotEmpty("y", coefs)) {
        DataMap::const_iterator iy = coefs.find("y");
        if (fs == -1)
            fs = iy->second.getFunctionSpace().getTypeCode();
        else if (fs != iy->second.getFunctionSpace().getTypeCode())
            throw ValueError("assemblePDEBoundary: coefficient function spaces don't match");
    }
    if (fs==-1) {
        return;
    }

    if (fs != FaceElements && fs != ReducedFaceElements)
        throw ValueError("assemblePDEBoundary: illegal function space type for coefficients");

    int numEq, numComp;
    if (!mat) {
        if (rhs.isEmpty()) {
            numEq=numComp=1;
        } else {
            numEq=numComp=rhs.getDataPointSize();
        }
    } else {
        if (!rhs.isEmpty() && rhs.getDataPointSize() != mat->getRowBlockSize())
            throw ValueError("assemblePDEBoundary: matrix row block size and number of components of right hand side don't match");
        numEq = mat->getRowBlockSize();
        numComp = mat->getColumnBlockSize();
    }

    if (numEq != numComp)
        throw ValueError("assemblePDEBoundary: number of equations and number of solutions don't match");

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

void RipleyDomain::assemblePDEDirac(escript::AbstractSystemMatrix* mat,
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
    int nEq, nComp;
    if (!mat) {
        if (rhs.isEmpty()) {
            nEq=nComp=1;
        } else {
            nEq=nComp=rhs.getDataPointSize();
        }
    } else {
        if (!rhs.isEmpty() && rhs.getDataPointSize() != mat->getRowBlockSize())
            throw ValueError("assemblePDEDirac: matrix row block size "
                    "and number of components of right hand side don't match");
        nEq = mat->getRowBlockSize();
        nComp = mat->getColumnBlockSize();
    }

    rhs.requireWrite();
    for (int i = 0; i < m_diracPoints.size(); i++) { //only for this rank
        const index_t dof = getDofOfNode(m_diracPoints[i].node);
        if (yNotEmpty) {
            const double *EM_F = y.getSampleDataRO(i);
            double *F_p = rhs.getSampleDataRW(0);
            if (dof < getNumDOF()) {
                for (index_t eq = 0; eq < nEq; eq++) {
                    F_p[INDEX2(eq, dof, nEq)] += EM_F[eq];
                }
            }
        }
        if (dNotEmpty) {
            const IndexVector rowIndex(1, dof);
            const double *EM_S = d.getSampleDataRO(i);
            vector<double> contents(EM_S, EM_S+nEq*nEq*nComp);
            addToSystemMatrix(mat, rowIndex, nEq, contents);
        }
    }
}

bool RipleyDomain::probeInterpolationAcross(int fsType_source,
                      const escript::AbstractDomain&, int fsType_target) const
{
    return false;
}

void RipleyDomain::interpolateAcross(escript::Data& target,
                                     const escript::Data& source) const
{
    throw NotImplementedError("interpolateAcross() not supported");
}

// Expecting ("gaussian", radius, sigma)
bool RipleyDomain::supportsFilter(const bp::tuple& t) const
{
    if (len(t) == 0) { // so we can handle unfiltered randoms
        return true;
    }
    if (len(t) != 3) {
        return false;
    }
    bp::extract<string> ex(t[0]);
    if (!ex.check() || (ex() != "gaussian")) {
        return false;
    }
    if (! bp::extract<unsigned int>(t[1]).check()) {
        return false;
    }
    return bp::extract<double>(t[2]).check();
}

void RipleyDomain::addPoints(const vector<double>& coords,
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

} // end of namespace ripley


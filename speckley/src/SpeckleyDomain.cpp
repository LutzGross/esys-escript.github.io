
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

#include <speckley/SpeckleyDomain.h>
#include <speckley/domainhelpers.h>

#include <escript/Data.h>
#include <escript/DataTypes.h>
#include <escript/DataFactory.h>
#include <escript/FunctionSpaceFactory.h>
#include <escript/index.h>

#include <iomanip>
#include <iostream>

namespace speckley {

namespace bp = boost::python;
using namespace std;
using escript::ValueError;

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

SpeckleyDomain::SpeckleyDomain(dim_t dim, int order, escript::JMPI jmpi) :
    escript::AbstractContinuousDomain(jmpi),
    m_numDim(dim),
    m_status(0),
    m_order(order)
{
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
                && m_elementTags==o->m_elementTags);
    }
    return false;
}

bool SpeckleyDomain::isValidFunctionSpaceType(int fsType) const
{
    switch (fsType) {
        case DegreesOfFreedom:
        case Nodes:
        case Elements:
        case ReducedElements:
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
            return pair<int,dim_t>(1, getNumNodes());
        case DegreesOfFreedom:
            return pair<int,dim_t>(1, getNumDOF());
        case Elements:
            return pair<int,dim_t>(ptsPerSample, getNumElements());
        case ReducedElements:
            return pair<int,dim_t>(1, getNumElements());
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

    There is also a set of lines. Interpolation is possible down a line but not
    between lines.
    class 0 and 1 belong to all lines so aren't considered.
    line 0: class 2
    line 1: class 3,4

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
            default:
                return false;
        }
    }
    int numLines=hasline[0]+hasline[1];

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
            return (fsType_target==Elements || fsType_target==ReducedElements || fsType_target==Nodes);
        case Points:
            return (fsType_target==fsType_source);
        case ReducedElements:
            return (fsType_target==Elements || fsType_target==Nodes);

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
            return (fsType_target==ReducedElements) ? -1 : 0;
        case ReducedElements:
            return (fsType_target==Elements) ? 1 : 0;
        case Points:
            return 0;  // other case caught by the if above
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
    else if (targetDomain != *this)
        throw SpeckleyException("Illegal domain of interpolation target");
    else if (in.isComplex() != target.isComplex())
        throw ValueError("interpolateOnDomain: complexity of input and output data must match.");

    stringstream msg;
    msg << "interpolateOnDomain() not implemented for function space "
        << functionSpaceTypeAsString(in.getFunctionSpace().getTypeCode())
        << " -> "
        << functionSpaceTypeAsString(target.getFunctionSpace().getTypeCode());

    const int inFS = in.getFunctionSpace().getTypeCode();
    const int outFS = target.getFunctionSpace().getTypeCode();

    target.requireWrite();
    // simplest case: 1:1 copy
    if (inFS==outFS) {
        if (in.isComplex())
            copyData<cplx_t>(target, in);
        else
            copyData<real_t>(target, in);
    } else if ((inFS==Elements || inFS==ReducedElements) && outFS==Nodes) {
        interpolateElementsOnNodes(target, in);
    } else if (inFS==ReducedElements && outFS==Elements) {
        if (in.isComplex())
            multiplyData<cplx_t>(target, in);
        else
            multiplyData<real_t>(target, in);
    } else if (inFS==Elements && outFS==ReducedElements) {
        reduceElements(target, in);
    } else {
        switch (inFS) {
            case Nodes:
                switch (outFS) {
                    case DegreesOfFreedom:
                    case Nodes:
                        if (in.isComplex())
                            copyData<cplx_t>(target, in);
                        else
                            copyData<real_t>(target, in);
                        break;
                    case Elements:
                        interpolateNodesOnElements(target, in, false);
                        break;
                    case ReducedElements:
                        interpolateNodesOnElements(target, in, true);
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
                switch (outFS) {
                    case Nodes:
                    case DegreesOfFreedom:
                        if (in.isComplex())
                            copyData<cplx_t>(target, in);
                        else
                            copyData<real_t>(target, in);
                        break;

                    case Elements:
                    case ReducedElements:
                        if (getMPISize()==1) {
                            interpolateNodesOnElements(target, in, outFS==ReducedElements);
                        } else {
                            escript::Data contIn(in, escript::continuousFunction(*this));
                            interpolateNodesOnElements(target, contIn, outFS==ReducedElements);
                        }
                        break;

                    default:
                        throw SpeckleyException(msg.str());
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
                throw SpeckleyException(msg.str());
        }
    }
}

escript::Data SpeckleyDomain::getX() const
{
    return escript::continuousFunction(*this).getX();
}

#ifdef ESYS_HAVE_BOOST_NUMPY
boost::python::numpy::ndarray SpeckleyDomain::getNumpyX() const
{
    return continuousFunction(*this).getNumpyX();
}

boost::python::numpy::ndarray SpeckleyDomain::getConnectivityInfo() const
{
    throw SpeckleyException("This feature is currently not supported by Speckley.");
}
#endif

escript::Data SpeckleyDomain::getNormal() const
{
    throw SpeckleyException("Speckley doesn't support getNormal");
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
        case Nodes:
        case Elements:
        case ReducedElements:
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
        case Nodes:
        case Elements:
            break;
        default: {
            throw SpeckleyException("setToGradient: only supported for nodal input data");
        }
    }

    if (grad.isComplex() != arg.isComplex())
        throw SpeckleyException("setToGradient: complexity of input and output must match");

    if (getMPISize() > 1) {
        if (arg.getFunctionSpace().getTypeCode()==DegreesOfFreedom) {
            escript::Data contArg(arg, escript::continuousFunction(*this));
            assembleGradient(grad, contArg);
        } else {
            assembleGradient(grad, arg);
        }
    } else {
        assembleGradient(grad, arg);
    }
}

void SpeckleyDomain::setToIntegrals(vector<real_t>& integrals,
                                    const escript::Data& arg) const
{
    setToIntegralsWorker<real_t>(integrals, arg);
}

void SpeckleyDomain::setToIntegrals(vector<cplx_t>& integrals,
                                    const escript::Data& arg) const
{
    setToIntegralsWorker<cplx_t>(integrals, arg);
}

template<typename Scalar>
void SpeckleyDomain::setToIntegralsWorker(vector<Scalar>& integrals,
                                          const escript::Data& arg) const
{
    const SpeckleyDomain& argDomain=dynamic_cast<const SpeckleyDomain&>(
            *(arg.getFunctionSpace().getDomain()));
    if (argDomain != *this)
        throw SpeckleyException("setToIntegrals: illegal domain of integration kernel");

    switch (arg.getFunctionSpace().getTypeCode()) {
        case Nodes:
        case DegreesOfFreedom:
            {
                escript::Data funcArg(arg, escript::function(*this));
                assembleIntegrate(integrals, funcArg);
            }
            break;
        case Points:
        case Elements:
        case ReducedElements:
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
        case DegreesOfFreedom:
            return false;
        case Elements:
        case ReducedElements:
        case Points:
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
        case Points:
            return true;
        case DegreesOfFreedom:
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
    vector<int>* target=NULL;
    dim_t num=0;

    switch(fsType) {
        case Nodes:
            num=getNumNodes();
            target=&m_nodeTags;
            break;
        case Elements:
            num=getNumElements();
            target=&m_elementTags;
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

int SpeckleyDomain::getTagFromSampleNo(int fsType, index_t sampleNo) const
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
        case Nodes:
            return m_nodeTagsInUse.size();
        case Elements:
        case ReducedElements:
            return m_elementTagsInUse.size();
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
        case Nodes:
            return &m_nodeTagsInUse[0];
        case Elements:
        case ReducedElements:
            return &m_elementTagsInUse[0];
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

int SpeckleyDomain::getSystemMatrixTypeId(const boost::python::object& options) const
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
        if (isNotEmpty("X", coefs) || isNotEmpty("du", coefs)
                    || isNotEmpty("Y", coefs))
            throw SpeckleyException(
                    "addPDEToRHS: right hand side coefficients are provided "
                    "but no right hand side vector given");
        else
            return;
    }

    assemblePDE(NULL, rhs, coefs, assembler);
    assemblePDEBoundary(NULL, rhs, coefs, assembler);
    assemblePDEDiracWrap(NULL, rhs, coefs, assembler);
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
    throw SpeckleyException("Speckley domains do not support transport problems");
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
template<typename Scalar>
void SpeckleyDomain::copyData(escript::Data& out, const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    const dim_t numSamples = in.getNumSamples();
    const Scalar zero = static_cast<Scalar>(0);
    out.requireWrite();
#pragma omp parallel for
    for (index_t i = 0; i < numSamples; i++) {
        const Scalar* src = in.getSampleDataRO(i, zero);
        copy(src, src+numComp, out.getSampleDataRW(i, zero));
    }
}

//protected
template<typename Scalar>
void SpeckleyDomain::multiplyData(escript::Data& out, const escript::Data& in) const
{
    const dim_t numComp = in.getDataPointSize();
    const dim_t dpp = out.getNumDataPointsPerSample();
    const dim_t numSamples = in.getNumSamples();
    const Scalar zero = static_cast<Scalar>(0);
    out.requireWrite();
#pragma omp parallel for
    for (index_t i=0; i<numSamples; i++) {
        const Scalar* src = in.getSampleDataRO(i, zero);
        Scalar* dest = out.getSampleDataRW(i, zero);
        for (index_t c=0; c<numComp; c++) {
            for (index_t q=0; q<dpp; q++)
                dest[c+q*numComp] = src[c];
        }
    }
}

//protected
void SpeckleyDomain::updateTagsInUse(int fsType) const
{
    std::vector<int>* tagsInUse=NULL;
    const std::vector<int>* tags=NULL;
    switch(fsType) {
        case Nodes:
            tags=&m_nodeTags;
            tagsInUse=&m_nodeTagsInUse;
            break;
        case Elements:
            tags=&m_elementTags;
            tagsInUse=&m_elementTagsInUse;
            break;
        case Points:
            throw SpeckleyException("updateTagsInUse for Speckley dirac points "
                    "not supported");
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
    if (rhs.isEmpty() && (isNotEmpty("X", coefs) || isNotEmpty("du", coefs))
            && isNotEmpty("Y", coefs))
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
    escript::Data temp(0., rhs.getDataPointShape(), rhs.getFunctionSpace(),
            rhs.actsExpanded());
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
        assembler->assemblePDESingle(mat, temp, coefs);
    } else {
        assembler->assemblePDESystem(mat, temp, coefs);
    }
#ifdef ESYS_MPI
    balanceNeighbours(temp, false); //summation without averaging
#endif
    rhs += temp; //only now add to rhs because otherwise rhs distorts
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
        assembler->assemblePDEBoundarySingle(mat, rhs, coefs);
    } else {
        assembler->assemblePDEBoundarySystem(mat, rhs, coefs);
    }
}

void SpeckleyDomain::assemblePDEDiracWrap(escript::AbstractSystemMatrix* mat,
                                    escript::Data& rhs, const DataMap& coefs,
                                    Assembler_ptr assembler) const
{
    bool complexpde = isComplexCoef("d_dirac", coefs) || isComplexCoef("D",coefs)
                    || isComplexCoef("y_dirac", coefs) || isComplexCoef("Y",coefs);
    if(complexpde)
        assembleComplexPDEDirac(mat,rhs,coefs,assembler);
    else
        assemblePDEDirac(mat,rhs,coefs,assembler);
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

    rhs.requireWrite();
    for (int i = 0; i < m_diracPoints.size(); i++) { //only for this rank
        const IndexVector rowIndex(1, m_diracPoints[i].node);
        if (yNotEmpty) {
            const double *EM_F = y.getSampleDataRO(i);
            double *F_p = rhs.getSampleDataRW(0);
            for (index_t eq = 0; eq < nEq; eq++) {
                F_p[INDEX2(eq, rowIndex[0], nEq)] += EM_F[INDEX2(eq,i,nEq)];
            }
        }
        if (dNotEmpty) {
            throw SpeckleyException("Rectangle::assemblePDEDirac currently doesn't support d");
//            const double *EM_S = d.getSampleDataRO(i);
//            std::vector<double> contents(EM_S,
//                        EM_S+nEq*nEq*nComp*rowIndex.size()); //TODO: this will break with MPI
//            addToSystemMatrix(mat, rowIndex, nEq, rowIndex, nComp, contents);
        }
    }
}

void SpeckleyDomain::assembleComplexPDEDirac(escript::AbstractSystemMatrix* mat,
                                    escript::Data& rhs, const DataMap& coefs,
                                    Assembler_ptr assembler) const
{
    bool yNotEmpty = isNotEmpty("y_dirac", coefs);
    bool dNotEmpty = isNotEmpty("d_dirac", coefs);
    if (!(yNotEmpty || dNotEmpty)) {
        return;
    }
    escript::Data d = unpackData("d_dirac", coefs);
    escript::Data yy = unpackData("y_dirac", coefs);

    // shallow copy
    // escript::Data d = Data(dd);
    escript::Data y = escript::Data(yy);

    // complicate things
    // if(!d.isEmpty())
    //     d.complicate();
    if(!y.isEmpty())
        y.complicate();
    if(!rhs.isEmpty())
        rhs.complicate();

    escript::DataTypes::cplx_t cdummy;

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

    rhs.requireWrite();
    for (int i = 0; i < m_diracPoints.size(); i++) { //only for this rank
        const IndexVector rowIndex(1, m_diracPoints[i].node);
        if (yNotEmpty) {
            const std::complex<double> *EM_F = y.getSampleDataRO(i, cdummy);
            std::complex<double> *F_p = rhs.getSampleDataRW(0, cdummy);
            for (index_t eq = 0; eq < nEq; eq++) {
                F_p[INDEX2(eq, rowIndex[0], nEq)] += EM_F[INDEX2(eq,i,nEq)];
            }
        }
        if (dNotEmpty) {
            throw SpeckleyException("Rectangle::assemblePDEDirac currently doesn't support d");
        }
    }
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
        } else if (m_mpiInfo->size == 1) {
            throw SpeckleyException("Dirac point unmapped, implementation problem in Speckley::addPoints");
        }
    }
}

} // end of namespace speckley

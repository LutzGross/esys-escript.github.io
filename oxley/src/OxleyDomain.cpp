/*****************************************************************************
*
* Copyright (c) 2003-2019 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include <string>

#include <oxley/OxleyDomain.h>
#include <oxley/OxleyData.h>

#include <escript/Data.h>
#include <escript/DataFactory.h>
#include <escript/EsysMPI.h>
#include <escript/FunctionSpaceFactory.h>
// #include <escript/index.h>
#include <escript/SolverOptions.h>

#ifdef ESYS_HAVE_PASO
#include <paso/SystemMatrix.h>
#endif
#ifdef ESYS_HAVE_TRILINOS
#include <trilinoswrap/CrsMatrixWrapper.h>
#include <trilinoswrap/TrilinosMatrixAdapter.h>
#include <trilinoswrap/types.h>
#include <Teuchos_Comm.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#include <PyTrilinos_Tpetra_Util.hpp>
#endif

namespace bp = boost::python;

using namespace std;

namespace oxley {

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

    /**
       \brief
       Constructor
    */
    OxleyDomain::OxleyDomain(dim_t dim, int order)
    {

        // These two statements configure the level of verbosity used by p4est
        sc_set_log_defaults(NULL, NULL, LOG_LEVEL);
        // p4est_init(NULL, LOG_LEVEL);

        // Register the converter used by boost
        // boost::python::to_python_converter<Teuchos::RCP<crs_matrix_type>,convert_Teuchos_RCP_to_Python_tuple>();
    }

    /**
       \brief
       Destructor
    */
    OxleyDomain::~OxleyDomain(){

    }

    bool OxleyDomain::isValidFunctionSpaceType(int fsType) const
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

    std::string OxleyDomain::functionSpaceTypeAsString(int fsType) const
    {
        switch (fsType) {
        case DegreesOfFreedom:
            return "Oxley_DegreesOfFreedom [Solution(domain)]";
        case ReducedDegreesOfFreedom:
            return "Oxley_ReducedDegreesOfFreedom [ReducedSolution(domain)]";
        case Nodes:
            return "Oxley_Nodes [ContinuousFunction(domain)]";
        case ReducedNodes:
            return "Oxley_ReducedNodes [ReducedContinuousFunction(domain)]";
        case Elements:
            return "Oxley_Elements [Function(domain)]";
        case ReducedElements:
            return "Oxley_ReducedElements [ReducedFunction(domain)]";
        case FaceElements:
            return "Oxley_FaceElements [FunctionOnBoundary(domain)]";
        case ReducedFaceElements:
            return "Oxley_ReducedFaceElements [ReducedFunctionOnBoundary(domain)]";
        case Points:
            return "Oxley_Points [DiracDeltaFunctions(domain)]";
        default:
            break;
    }
    return "Invalid function space type code";
    }

    bool OxleyDomain::operator==(const AbstractDomain& other) const
    {
        throw OxleyException("currently not implemented");
    }

    int OxleyDomain::getTagFromSampleNo(int fsType, index_t sampleNo) const
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

    void OxleyDomain::interpolateOnDomain(escript::Data& target, const escript::Data& in) const
    {

        const OxleyDomain& inDomain=dynamic_cast<const OxleyDomain&>(*(in.getFunctionSpace().getDomain()));
        const OxleyDomain& targetDomain=dynamic_cast<const OxleyDomain&>(*(target.getFunctionSpace().getDomain()));
        if (inDomain != *this)
            throw ValueError("Illegal domain of interpolant");
        if (targetDomain != *this)
            throw ValueError("Illegal domain of interpolation target");
        if (target.isComplex() != in.isComplex())
            throw ValueError("Complexity of input and output must match");

        const int inFS = in.getFunctionSpace().getTypeCode();
        const int outFS = target.getFunctionSpace().getTypeCode();

        stringstream msg;
        msg << "interpolateOnDomain() not implemented for function space "
            << functionSpaceTypeAsString(inFS) 
            << " -> "
            << functionSpaceTypeAsString(outFS);

        // simplest case: 1:1 copy
        if (inFS==outFS) {
            if (in.isComplex())
                copyData<cplx_t>(target, in);
            else
                copyData<real_t>(target, in);
        // not allowed: reduced nodes/DOF->non-reduced nodes/DOF
        } else if ((inFS==ReducedNodes || inFS==ReducedDegreesOfFreedom)
                && (outFS==Nodes || outFS==DegreesOfFreedom)) {
            throw ValueError("interpolateOnDomain: Cannot interpolate reduced data to non-reduced data.");
        } else if ((inFS==Elements && outFS==ReducedElements)
                || (inFS==FaceElements && outFS==ReducedFaceElements)) {
            if (in.actsExpanded()) {
                if (in.isComplex())
                    averageData<cplx_t>(target, in);
                else
                    averageData<real_t>(target, in);
            } else {
                if (in.isComplex())
                    copyData<cplx_t>(target, in);
                else
                    copyData<real_t>(target, in);
            }
        } else if ((inFS==ReducedElements && outFS==Elements)
                || (inFS==ReducedFaceElements && outFS==FaceElements)) {
            if (in.isComplex())
                multiplyData<cplx_t>(target, in);
            else
                multiplyData<real_t>(target, in);
        } else {
            switch (inFS) {
                case Nodes:
                case ReducedNodes:
                    switch (outFS) {
                        case DegreesOfFreedom:
                        case ReducedDegreesOfFreedom:
                            if (getMPISize()==1) {
                                if (in.isComplex())
                                    copyData<cplx_t>(target, in);
                                else
                                    copyData<real_t>(target, in);
                            } else {
                                if (in.isComplex())
                                    throw escript::NotImplementedError("nodesToDOF not implemented for complex Data");
                                else
                                    nodesToDOF(target, in);
                            }
                            break;

                        case Nodes:
                        case ReducedNodes:
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
                            throw escript::NotImplementedError(msg.str());
                    }
                    break;

                case DegreesOfFreedom:
                case ReducedDegreesOfFreedom:
                    switch (outFS) {
                        case Nodes:
                        case ReducedNodes:
                            if (getMPISize()==1)
                                if (in.isComplex())
                                    copyData<cplx_t>(target, in);
                                else
                                    copyData<real_t>(target, in);
                            else
                                if (in.isComplex())
                                    dofToNodes<cplx_t>(target, in);
                                else
                                    dofToNodes<real_t>(target, in);
                            break;

                        case DegreesOfFreedom:
                        case ReducedDegreesOfFreedom:
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
                            throw escript::NotImplementedError(msg.str());
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
                    throw escript::NotImplementedError(msg.str());
            }
        }
    }

    //protected
    template<typename Scalar>
    void OxleyDomain::dofToNodes(escript::Data& out, const escript::Data& in) const
    {
        // expand data object if necessary
        const_cast<escript::Data*>(&in)->expand();
        const dim_t numComp = in.getDataPointSize();
        const dim_t numNodes = getNumNodes();
        const Scalar zero = static_cast<Scalar>(0);
        out.requireWrite();
    #ifdef ESYS_HAVE_TRILINOS
        using namespace esys_trilinos;

        const_TrilinosGraph_ptr graph(getTrilinosGraph());
        Teuchos::RCP<const MapType> colMap;
        Teuchos::RCP<const MapType> rowMap;
        MapType colPointMap;
        MapType rowPointMap;

        if (numComp > 1) {
            colPointMap = BlockVectorType<Scalar>::makePointMap(
                                            *graph->getColMap(), numComp);
            rowPointMap = BlockVectorType<Scalar>::makePointMap(
                                            *graph->getRowMap(), numComp);
            colMap = Teuchos::rcpFromRef(colPointMap);
            rowMap = Teuchos::rcpFromRef(rowPointMap);
        } else {
            colMap = graph->getColMap();
            rowMap = graph->getRowMap();
        }

        const ImportType importer(rowMap, colMap);
        const Teuchos::ArrayView<const Scalar> localIn(
                      in.getSampleDataRO(0, zero), in.getNumDataPoints()*numComp);
        Teuchos::RCP<VectorType<Scalar> > lclData = rcp(new VectorType<Scalar>(
                                            rowMap, localIn, localIn.size(), 1));
        Teuchos::RCP<VectorType<Scalar> > gblData = rcp(new VectorType<Scalar>(
                                            colMap, 1));
        gblData->doImport(*lclData, importer, Tpetra::INSERT);
        Teuchos::ArrayRCP<const Scalar> gblArray(gblData->getData(0));

    #pragma omp parallel for
        for (index_t i = 0; i < numNodes; i++) {
            const Scalar* src = &gblArray[getDofOfNode(i) * numComp];
            copy(src, src+numComp, out.getSampleDataRW(i, zero));
        }
    #elif defined(ESYS_HAVE_PASO)
        paso::Coupler_ptr<Scalar> coupler(new paso::Coupler<Scalar>(m_connector,
                                                             numComp, m_mpiInfo));
        coupler->startCollect(in.getSampleDataRO(0, zero));
        const dim_t numDOF = getNumDOF();
        const Scalar* buffer = coupler->finishCollect();

    #pragma omp parallel for
        for (index_t i = 0; i < numNodes; i++) {
            const index_t dof = getDofOfNode(i);
            const Scalar* src = (dof < numDOF ? in.getSampleDataRO(dof, zero)
                                              : &buffer[(dof - numDOF) * numComp]);
            copy(src, src+numComp, out.getSampleDataRW(i, zero));
        }
    #endif // ESYS_HAVE_PASO
    }


    //protected
    template<typename Scalar>
    void OxleyDomain::copyData(escript::Data& out, const escript::Data& in) const
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
    void OxleyDomain::averageData(escript::Data& out, const escript::Data& in) const
    {
        const dim_t numComp = in.getDataPointSize();
        const dim_t dpp = in.getNumDataPointsPerSample();
        const dim_t numSamples = in.getNumSamples();
        const Scalar zero = static_cast<Scalar>(0);
        out.requireWrite();
#pragma omp parallel for
        for (index_t i = 0; i < numSamples; i++) {
            const Scalar* src = in.getSampleDataRO(i, zero);
            Scalar* dest = out.getSampleDataRW(i, zero);
            for (index_t c = 0; c < numComp; c++) {
                Scalar res = zero;
                for (index_t q = 0; q < dpp; q++)
                    res += src[c+q*numComp];
                *dest++ = res / static_cast<real_t>(dpp);
            }
        }
    }

    //protected
    template<typename Scalar>
    void OxleyDomain::multiplyData(escript::Data& out, const escript::Data& in) const
    {
        const dim_t numComp = in.getDataPointSize();
        const dim_t dpp = out.getNumDataPointsPerSample();
        const dim_t numSamples = in.getNumSamples();
        const Scalar zero = static_cast<Scalar>(0);
        out.requireWrite();
#pragma omp parallel for
        for (index_t i = 0; i < numSamples; i++) {
            const Scalar* src = in.getSampleDataRO(i, zero);
            Scalar* dest = out.getSampleDataRW(i, zero);
            for (index_t c = 0; c < numComp; c++) {
                for (index_t q = 0; q < dpp; q++)
                    dest[c+q*numComp] = src[c];
            }
        }
    }

    bool OxleyDomain::probeInterpolationOnDomain(int fsType_source, int fsType_target) const
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

    signed char OxleyDomain::preferredInterpolationOnDomain(int fsType_source, int fsType_target) const
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

    bool OxleyDomain::commonFunctionSpace(const vector<int>& fs, int& resultcode) const
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

    escript::Data OxleyDomain::getX() const
    {
        throw OxleyException("programming error"); 
    }

    std::string OxleyDomain::getDescription() const
    {
        throw OxleyException("programming error");
    }

    escript::Data OxleyDomain::getNormal() const
    {
        return escript::functionOnBoundary(*this).getNormal();
    }

    escript::Data OxleyDomain::getSize() const
    {
        return escript::function(*this).getSize();
    }

    void OxleyDomain::setToX(escript::Data& arg) const
    {
        const OxleyDomain& argDomain=dynamic_cast<const OxleyDomain&>( *(arg.getFunctionSpace().getDomain()));
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

    void OxleyDomain::setToGradient(escript::Data& grad, const escript::Data& arg) const
    {
        const OxleyDomain& argDomain=dynamic_cast<const OxleyDomain&>( *(arg.getFunctionSpace().getDomain()));
        if (argDomain != *this)
            throw ValueError("setToGradient: Illegal domain of gradient argument");
        const OxleyDomain& gradDomain=dynamic_cast<const OxleyDomain&>(
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

    void OxleyDomain::setToIntegrals(vector<real_t>& integrals, const escript::Data& arg) const
    {
        setToIntegralsWorker<real_t>(integrals, arg);
    }

    void OxleyDomain::setToIntegrals(vector<cplx_t>& integrals, const escript::Data& arg) const
    {
        setToIntegralsWorker<cplx_t>(integrals, arg);
    }

    template<typename Scalar>
    void OxleyDomain::setToIntegralsWorker(std::vector<Scalar>& integrals,
                                            const escript::Data& arg) const
    {
        const OxleyDomain& argDomain = dynamic_cast<const OxleyDomain&>(
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
            case Points:
                {
                    assembleIntegrate(integrals, arg);
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

    bool OxleyDomain::isCellOriented(int fsType) const
    {
        throw OxleyException("currently not implemented"); // This is temporary
    }

    int OxleyDomain::getNumberOfTagsInUse(int fsType) const
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

    const int* OxleyDomain::borrowListOfTagsInUse(int fsType) const
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

    bool OxleyDomain::canTag(int fsType) const
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

    void OxleyDomain::updateTagsInUse(int fsType) const
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
                throw escript::NotImplementedError("updateTagsInUse for Oxley dirac points"
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

    void OxleyDomain::writeToVTK(std::string filename, bool writeTagInfo) const
    {
        throw OxleyException("unknown error");
    }

    void OxleyDomain::refineMesh(std::string algorithm)
    {
        throw OxleyException("unknown error");
    }

    void OxleyDomain::refineBoundary(std::string boundary, double dx)
    {
        throw OxleyException("unknown error");
    }

    void OxleyDomain::refineRegion(double x0, double x1, double y0, double y1)
    {
        throw OxleyException("unknown error");   
    }

    void OxleyDomain::refinePoint(double x0, double y0)
    {
        throw OxleyException("unknown error");   
    }

    void OxleyDomain::refineCircle(double x0, double y0, double r)
    {
        throw OxleyException("unknown error");   
    }

    void OxleyDomain::setTags(int fsType, int newTag, const escript::Data& mask) const
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
                // throw OxleyException("not implemented yet"); //TODO
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

    pair<int,dim_t> OxleyDomain::getDataShape(int fsType) const
    {
        const int ptsPerSample = (m_numDim == 2 ? 4 : 8);
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
        msg << "getDataShape: Invalid function space type " << fsType << " for " << getDescription();
        throw ValueError(msg.str());
    }

    std::string OxleyDomain::showTagNames() const
    {
        stringstream ret;
        TagMap::const_iterator it;
        for (it=m_tagMap.begin(); it!=m_tagMap.end(); it++) {
            if (it!=m_tagMap.begin()) ret << ", ";
            ret << it->first;
        }
        return ret.str();
    }

    //protected
    template<>
    void OxleyDomain::addToSystemMatrix<real_t>(escript::AbstractSystemMatrix* mat,
                                             const IndexVector& nodes, dim_t numEq,
                                             const DoubleVector& array) const
    {
    #ifdef ESYS_HAVE_PASO
        paso::SystemMatrix* psm = dynamic_cast<paso::SystemMatrix*>(mat);
        if (psm) {
            addToPasoMatrix(psm, nodes, numEq, array);
            return;
        }
    #endif
    #ifdef ESYS_HAVE_TRILINOS
        esys_trilinos::TrilinosMatrixAdapter* tm = dynamic_cast<esys_trilinos::TrilinosMatrixAdapter*>(mat);
        if (tm) {
            tm->add(nodes, array);
            return;
        }
    #endif
        throw OxleyException("addToSystemMatrix: unknown system matrix type");
    }

    //protected
    template<>
    void OxleyDomain::addToSystemMatrix<cplx_t>(escript::AbstractSystemMatrix* mat,
                                             const IndexVector& nodes, dim_t numEq,
                                             const vector<cplx_t>& array) const
    {
    #ifdef ESYS_HAVE_MUMPS
        paso::SystemMatrix* psm = dynamic_cast<paso::SystemMatrix*>(mat);
        if (psm) {
            addToPasoMatrix(psm, nodes, numEq, array);
            return;
        }
    #endif
    #ifdef ESYS_HAVE_TRILINOS
        esys_trilinos::TrilinosMatrixAdapter* tm = dynamic_cast<esys_trilinos::TrilinosMatrixAdapter*>(mat);
        if (tm) {
            tm->add(nodes, array);
            return;
        }
    #endif
        throw OxleyException("addToSystemMatrix: require Trilinos or MUMPS matrices for "
                              "complex-valued assembly!");
    }

    int OxleyDomain::getSystemMatrixTypeId(const bp::object& options) const
    {
        const escript::SolverBuddy& sb = bp::extract<escript::SolverBuddy>(options);
        int package = sb.getPackage();
        escript::SolverOptions method = sb.getSolverMethod();
#ifdef ESYS_HAVE_TRILINOS
        bool isDirect = escript::isDirectSolver(method);
#endif

        // the configuration of oxley should have taken care that we have either
        // paso or trilinos so here's how we prioritize
#if defined(ESYS_HAVE_PASO) && defined(ESYS_HAVE_TRILINOS)
        // we have Paso & Trilinos so use Trilinos for parallel direct solvers
        if (package == escript::SO_DEFAULT) {
            if ((method == escript::SO_METHOD_DIRECT && getMPISize() > 1)
                    || isDirect
                    || sb.isComplex()) {
                package = escript::SO_PACKAGE_TRILINOS;
            }
        }
#endif
#ifdef ESYS_HAVE_PASO
        if (package == escript::SO_DEFAULT)
            package = escript::SO_PACKAGE_PASO;
#endif
#ifdef ESYS_HAVE_TRILINOS
        if (package == escript::SO_DEFAULT)
            package = escript::SO_PACKAGE_TRILINOS;
#endif
        if (package == escript::SO_PACKAGE_TRILINOS) {
#ifdef ESYS_HAVE_TRILINOS
            int type = (int) SMT_TRILINOS;
            if (sb.isComplex())
                type |= (int) SMT_COMPLEX;
            // This is required because MueLu (AMG) and Amesos2 (direct) do not
            // support block matrices at this point. Remove if they ever do...
            if (sb.getPreconditioner() == escript::SO_PRECONDITIONER_AMG ||
                    sb.getPreconditioner() == escript::SO_PRECONDITIONER_ILUT ||
                    isDirect) {
                type |= (int) SMT_UNROLL;
            }
            return type;
#else
            throw OxleyException("Trilinos requested but not built with Trilinos.");
#endif
        }
#ifdef ESYS_HAVE_PASO
        if (sb.isComplex()) {
            throw escript::NotImplementedError("Paso does not support complex-valued matrices");
        }
        // in all other cases we use PASO
        return (int)SMT_PASO | paso::SystemMatrix::getSystemMatrixTypeId(method, sb.getPreconditioner(), sb.getPackage(), sb.isSymmetric(), m_mpiInfo);
#else
        throw OxleyException("Unable to find a working solver library!");
#endif
    }


    int OxleyDomain::getTransportTypeId(int solver, int preconditioner, int package, bool symmetry) const
    {
#ifdef ESYS_USE_PASO
        return paso::TransportProblem::getTypeId(solver, preconditioner, package, symmetry, m_mpiInfo);
#else
        throw OxleyException("Transport solvers require Paso but Oxley was not compiled with Paso!");
#endif
    }

    
    escript::ASM_ptr OxleyDomain::newSystemMatrix(int row_blocksize,
        const escript::FunctionSpace& row_functionspace, int column_blocksize,
        const escript::FunctionSpace& column_functionspace, int type) const
    {
        bool reduceRowOrder=false;
        bool reduceColOrder=false;
        // is the domain right?
        const OxleyDomain& row_domain=dynamic_cast<const OxleyDomain&>(*(row_functionspace.getDomain()));
        if (row_domain != *this)
            throw ValueError("newSystemMatrix: domain of row function space does not match the domain of matrix generator");
        const OxleyDomain& col_domain=dynamic_cast<const OxleyDomain&>(*(column_functionspace.getDomain()));
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
            throw OxleyException("eScript does not support CUDA.");
        } else if (type & (int)SMT_TRILINOS) {
    #ifdef ESYS_HAVE_TRILINOS
            esys_trilinos::const_TrilinosGraph_ptr graph(getTrilinosGraph());
            bool isComplex = (type & (int)SMT_COMPLEX);
            bool unroll = (type & (int)SMT_UNROLL);
            escript::ASM_ptr sm(new esys_trilinos::TrilinosMatrixAdapter(m_mpiInfo, row_blocksize, row_functionspace, graph, isComplex, unroll));
            return sm;
    #else
            throw OxleyException("newSystemMatrix: oxley was not compiled with Trilinos support so the Trilinos solver stack cannot be used.");
    #endif
        } else if (type & (int)SMT_PASO) {
    #ifdef ESYS_HAVE_PASO
            paso::SystemMatrixPattern_ptr pattern(getPasoMatrixPattern(reduceRowOrder, reduceColOrder));
            type -= (int)SMT_PASO;
            escript::ASM_ptr sm(new paso::SystemMatrix(type, pattern, row_blocksize, column_blocksize, false, row_functionspace, column_functionspace));
            return sm;
    #else
            throw OxleyException("newSystemMatrix: oxley was not compiled with Paso support so the Paso solver stack cannot be used.");
    #endif
        } else {
            throw OxleyException("newSystemMatrix: unknown matrix type ID");
        }
    }


#ifdef ESYS_HAVE_TRILINOS
//protected
esys_trilinos::const_TrilinosGraph_ptr OxleyDomain::createTrilinosGraph(
                                            const IndexVector& YaleRows,
                                            const IndexVector& YaleColumns) const
{
    using namespace esys_trilinos;

    const dim_t numMatrixRows = getNumDOF();

    IndexVector rowTemp(getNumDataPointsGlobal());
    if(getMPISize() == 1)
    {
    #pragma omp for
        for(long i = 0; i < getNumDataPointsGlobal(); i++)
            rowTemp[i] = i;
    }
    else
    {
        OxleyException("Not yet implemented"); //TODO
    }

    // rowMap
    // This is using the constructor on line 868 of file  Tpetra_Map_def.hpp.
    TrilinosMap_ptr rowMap(new MapType(getNumDataPointsGlobal(), rowTemp, 0, TeuchosCommFromEsysComm(m_mpiInfo->comm)));

    // colMap
    TrilinosMap_ptr colMap(new MapType(getNumDataPointsGlobal(), rowTemp, 0, TeuchosCommFromEsysComm(m_mpiInfo->comm)));
    
    // rowPtr
    const vector<IndexVector>& conns(getConnections(true));
    Teuchos::ArrayRCP<size_t> rowPtr(numMatrixRows+1);
    for (size_t i=0; i < numMatrixRows; i++) {
        rowPtr[i+1] = rowPtr[i] + conns[i].size();
    }
    Teuchos::ArrayRCP<LO> colInd(rowPtr[numMatrixRows]);

    // colInd
#pragma omp parallel for
    for (index_t i=0; i < numMatrixRows; i++) {
        copy(conns[i].begin(), conns[i].end(), &colInd[rowPtr[i]]);
    }

    // for(int i = 0; i < getNumDataPointsGlobal(); i++)
    //     std::cout << "myRows["<<i<<"]: " << rowTemp[i]<<std::endl;
    // for(int i = 0; i < getNumDataPointsGlobal(); i++)
    //     std::cout << "colMap["<<i<<"]: " << rowTemp[i]<<std::endl;
    // for(int i = 0; i < numMatrixRows+1; i++)
    //     std::cout << "rowPtr["<<i<<"]: " << rowPtr[i]<<std::endl;
    // for(int i = 0; i < rowPtr[numMatrixRows]; i++)
    //     std::cout << "colInd["<<i<<"]: " << colInd[i]<<std::endl;

    // params
    TrilinosGraph_ptr graph(new GraphType(rowMap, colMap, rowPtr, colInd));
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
    params->set("Optimize Storage", true);
    graph->fillComplete(rowMap, rowMap, params);
    return graph;
}
#endif

#ifdef ESYS_HAVE_PASO
    void OxleyDomain::addToPasoMatrix(paso::SystemMatrix* mat,
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
#endif // ESYS_HAVE_PASO


#ifdef ESYS_HAVE_PASO
    void OxleyDomain::createPasoConnector(const RankVector& neighbour, 
                                           const IndexVector& offsetInSharedSend,
                                           const IndexVector& offsetInSharedRecv,
                                           const IndexVector& sendShared,
                                           const IndexVector& recvShared)
    {
        const index_t* sendPtr = neighbour.empty() ? NULL : &sendShared[0];
        const index_t* recvPtr = neighbour.empty() ? NULL : &recvShared[0];
        paso::SharedComponents_ptr snd_shcomp(new paso::SharedComponents(getNumDOF(), neighbour, sendPtr, offsetInSharedSend));
        paso::SharedComponents_ptr rcv_shcomp(new paso::SharedComponents(getNumDOF(), neighbour, recvPtr, offsetInSharedRecv));
        m_connector.reset(new paso::Connector(snd_shcomp, rcv_shcomp));
    }
#endif

    //protected
#ifdef ESYS_HAVE_PASO
    paso::Pattern_ptr OxleyDomain::createPasoPattern(
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
#endif // ESYS_HAVE_PASO

void OxleyDomain::addToSystem(escript::AbstractSystemMatrix& mat,
                           escript::Data& rhs, const DataMap& coefs,
                           Assembler_ptr assembler) const
    {
        if (isNotEmpty("d_contact", coefs) || isNotEmpty("y_contact", coefs))
            throw ValueError(
                        "addToSystem: Oxley does not support contact elements");

        assemblePDE(&mat, rhs, coefs, assembler);
        assemblePDEBoundary(&mat, rhs, coefs, assembler);
        assemblePDEDirac(&mat, rhs, coefs, assembler);
    }

void OxleyDomain::makeZ()
    {
        if(z_needs_update)
        {
            const Tpetra::global_size_t tn = getNumNodes(); //Total number of nodes
            const Tpetra::global_size_t nh = 0.5*getNumHangingNodes(); // Number of hanging nodes
            const Tpetra::global_size_t nn = tn - nh;

            #ifdef OXLEY_PRINT_DEBUG_IZ
                std::cout << "total nodes=" << tn << ", number hanging nodes=" << nh << std::endl;
            #endif

            const esys_trilinos::GO indexBase = 0;

            auto comm = esys_trilinos::TeuchosCommFromEsysComm(m_mpiInfo->comm);
            cplx_t tmp_num(0.5,0);
            const cplx_t half = static_cast<cplx_t> (tmp_num);

            Teuchos::RCP<const Tpetra::Map<>> zRowMap    = Teuchos::rcp (new Tpetra::Map<>(nh, indexBase, comm));
            Teuchos::RCP<const Tpetra::Map<>> zColMap    = Teuchos::rcp (new Tpetra::Map<>(nn, indexBase, comm));
            Teuchos::RCP<const Tpetra::Map<>> zRangeMap  = Teuchos::rcp (new Tpetra::Map<>(nh, indexBase, comm));
            Teuchos::RCP<const Tpetra::Map<>> zDomainMap = Teuchos::rcp (new Tpetra::Map<>(nn, indexBase, comm));
            Teuchos::RCP<oxley::OxleyDomain::crs_matrix_type> z (new oxley::OxleyDomain::crs_matrix_type(zRowMap, zColMap, 4));

            for(int i = 0; i < getNumHangingNodes(); i++)
            {
                int a = hanging_faces[i].first;
                int b = hanging_faces[i].second;
                if(a>getNumNodes())
                {
                    int c = a;
                    a=b; 
                    b=c;
                }

                #ifdef OXLEY_PRINT_DEBUG_IZ
                    std::cout << "Z  element: (" << a-nn << ", " << b << ") = " << 0.5 << std::endl;
                #endif

                const esys_trilinos::GO gblRowAz = zRowMap->getGlobalElement(a-nn);
                const esys_trilinos::GO gblColBz = zColMap->getGlobalElement(b);
                z->insertGlobalValues(gblRowAz,
                                    Teuchos::tuple<esys_trilinos::GO>(gblColBz),
                                    Teuchos::tuple<cplx_t> (half));
                                    // Teuchos::tuple<scalar_type> (half));
            }

            // Tell the matrix that we are finished adding entries to it.
            z->fillComplete(zDomainMap,zRangeMap);
            pZ=&z;
            z_needs_update=false; 
        }
    }

void OxleyDomain::makeIZ()
    {    
        if(iz_needs_update)
        {
            const Tpetra::global_size_t tn = getNumNodes(); //Total number of nodes
            const Tpetra::global_size_t nh = 0.5*getNumHangingNodes(); // Number of hanging nodes
            const Tpetra::global_size_t nn = tn - nh;

            #ifdef OXLEY_PRINT_DEBUG_IZ
                std::cout << "total nodes=" << tn << ", number hanging nodes=" << nh << std::endl;
            #endif

            const esys_trilinos::GO indexBase = 0;
            auto comm = esys_trilinos::TeuchosCommFromEsysComm(m_mpiInfo->comm);

            Teuchos::RCP<const Tpetra::Map<>> izRowMap    = Teuchos::rcp (new Tpetra::Map<>(tn, indexBase, comm));
            Teuchos::RCP<const Tpetra::Map<>> izColMap    = Teuchos::rcp (new Tpetra::Map<>(nn, indexBase, comm));
            Teuchos::RCP<const Tpetra::Map<>> izRangeMap  = Teuchos::rcp (new Tpetra::Map<>(tn, indexBase, comm));
            Teuchos::RCP<const Tpetra::Map<>> izDomainMap = Teuchos::rcp (new Tpetra::Map<>(nn, indexBase, comm));
            Teuchos::RCP<oxley::OxleyDomain::crs_matrix_type> iz (new oxley::OxleyDomain::crs_matrix_type(izRowMap, izColMap, 4));
            /////////////////////////
            // Fill in iz
            /////////////////////////
            // This is I
            const cplx_t one  = static_cast<cplx_t> (1.0);
            const cplx_t half = static_cast<cplx_t> (0.5);
            for (esys_trilinos::LO lclRow = 0; lclRow < static_cast<esys_trilinos::LO>(nn); ++lclRow) 
            {
                const esys_trilinos::GO gblRow = izRowMap->getGlobalElement(lclRow);
                const esys_trilinos::GO gblCol = izColMap->getGlobalElement(lclRow);
                iz->insertGlobalValues(gblRow,
                                       Teuchos::tuple<esys_trilinos::GO>(gblCol),
                                       Teuchos::tuple<cplx_t>(one));
                #ifdef OXLEY_PRINT_DEBUG_IZ
                    std::cout << "iz element: (" << gblRow << ", " << gblRow << ") = " << 1.0 << std::endl;
                #endif
            }

            // This is Z
            for(int i = 0; i < getNumHangingNodes(); i++)
            {
                int a = hanging_faces[i].first;
                int b = hanging_faces[i].second;
                if(a>getNumNodes())
                {
                    int c = a;
                    a=b; 
                    b=c;
                }

                #ifdef OXLEY_PRINT_DEBUG_IZ
                    std::cout << "iz element: (" << a << ", " << b << ") = " << 0.5 << std::endl;
                #endif

                const esys_trilinos::GO gblRowA = izRowMap->getGlobalElement(a);
                const esys_trilinos::GO gblColB = izColMap->getGlobalElement(b);
                iz->insertGlobalValues(gblRowA,
                                       Teuchos::tuple<esys_trilinos::GO>(gblColB),
                                       Teuchos::tuple<cplx_t> (half));
            }

            // Tell the matrix that we are finished adding entries to it.
            iz->fillComplete(izDomainMap,izRangeMap);
            pIZ=&iz;
            iz_needs_update=false;
        }
    }

Teuchos::RCP<Tpetra::CrsMatrix<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> * OxleyDomain::getZ()
{
    return pZ;
}

Teuchos::RCP<Tpetra::CrsMatrix<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> * OxleyDomain::getIZ()
{
    return pIZ;
}

void OxleyDomain::finaliseA(escript::AbstractSystemMatrix& mat)
    {
        ////////////////////////////////////////////////
        if(getNumHangingNodes() > 0)
        {
            // Teuchos::RCP<Tpetra::CrsMatrix<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> iz_tmp(*pIZ);
            // if(iz_tmp->isFillComplete()==false)
            //     iz_tmp->fillComplete();

            escript::AbstractSystemMatrix * pMat = &mat;
            esys_trilinos::CrsMatrixWrapper<cplx_t> * cm = dynamic_cast<esys_trilinos::CrsMatrixWrapper<cplx_t>*>(pMat);
            
            if(cm)
            {
                // cm->IztAIz(iz_tmp);
                cm->IztAIz(*pIZ);
            }
        }
    }

void OxleyDomain::finaliseRhs(escript::Data& rhs)
    {
        #ifdef OXLEY_PRINT_DEBUG_ADDTOSYSTEM
            std::cout << "Before" << std::endl;
            rhs.print();
        #endif

        if(getNumHangingNodes() > 0)
        {
            using Teuchos::Array;
            using Teuchos::ArrayView;
            using Teuchos::ArrayRCP;
            using Teuchos::arcp;
            using Teuchos::RCP;
            using Teuchos::rcp;
            using Teuchos::tuple;
            typedef Tpetra::Map<> map_type;
            typedef Tpetra::Vector<>::scalar_type scalar_type;
            // typedef Tpetra::Vector<>::local_ordinal_type local_ordinal_type;
            typedef Tpetra::Vector<>::global_ordinal_type global_ordinal_type;
            // typedef Tpetra::CrsMatrix<> crs_matrix_type;
            const Tpetra::global_size_t tn = getNumNodes(); //Total number of nodes
            const Tpetra::global_size_t nh = 0.5*getNumHangingNodes(); // Number of hanging nodes
            const Tpetra::global_size_t nn = tn - nh;
            const global_ordinal_type indexBase = 0;
            auto comm = esys_trilinos::TeuchosCommFromEsysComm(m_mpiInfo->comm);
            

            /////////////////////////
            // Now do the multiplication
            /////////////////////////

            // RHS
            //recast rhs as a vector
            // long ndp = rhs.getNumDataPoints();
            RCP<const map_type> f_map = rcp(new map_type(nn, indexBase, comm));
            RCP<const map_type> g_map = rcp(new map_type(nh, indexBase, comm));
            //TODO replace these two vectors with a single Tpetra::MultiVector
            
            Tpetra::MultiVector<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT> f(f_map,true);
            Tpetra::MultiVector<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT> g(g_map,true);
            
            // Copy the data
            // vector f
            #pragma omp parallel for
            for(int i = 0; i < nn; i++)
            {
                const global_ordinal_type gblrow = i;
                const double *value = rhs.getSampleDataRO(i);
                f.replaceGlobalValue(gblrow,0,*value);
            }
            
            // vector g
            #pragma omp parallel for
            for(int i = nn; i < tn; i++)
            {
                const global_ordinal_type gblrow = i;
                const double *value = rhs.getSampleDataRO(i);
                g.replaceGlobalValue(gblrow,0,*value);
            }
            

            Teuchos::RCP<Tpetra::CrsMatrix<cplx_t,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>> Z(*pZ);
            Z->fillComplete();

            // multiplication using trilinos
            const scalar_type one = static_cast<scalar_type> (1.0);
            Z->apply(g,f,Teuchos::TRANS,one,one);

            auto result_view = f.getLocalViewHost();
            auto result_view_1d = Kokkos::subview(result_view, Kokkos::ALL(), 0);

            // write the new vector back into rhs
            escript::DataTypes::cplx_t dummy;
            rhs.requireWrite();
            #pragma omp parallel for
            for(int i = 0; i < nn; i++)
            {
                escript::DataTypes::cplx_t * value =  rhs.getSampleDataRW(i, dummy);
                *value=result_view_1d(i);
                #ifdef OXLEY_PRINT_DEBUG_IZ
                    std::cout << "rhs element: (" << i << ") = " << result_view_1d(i) << std::endl;
                #endif
            }
        }

        #ifdef OXLEY_PRINT_DEBUG_ADDTOSYSTEM
            std::cout << "After" << std::endl;
            rhs.print();
        #endif
    }

void OxleyDomain::addToSystemFromPython(escript::AbstractSystemMatrix& mat,
                                         escript::Data& rhs,
                                         const bp::list& data,
                                         Assembler_ptr assembler) const
    {
        DataMap mapping;
        tupleListToMap(mapping, data);
        addToSystem(mat, rhs, mapping, assembler);
    }


Assembler_ptr OxleyDomain::createAssemblerFromPython(const string type,
                                                const bp::list& options) const
    {
        DataMap mapping;
        tupleListToMap(mapping, options);
        return createAssembler(type, mapping);
    }

void OxleyDomain::addToRHS(escript::Data& rhs, const DataMap& coefs,
                            Assembler_ptr assembler) const
    {
        if (isNotEmpty("y_contact", coefs))
            throw ValueError(
                        "addPDEToRHS: Oxley does not support contact elements");

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

void OxleyDomain::addToRHSFromPython(escript::Data& rhs, const bp::list& data,
                                      Assembler_ptr assembler) const
    {
        DataMap mapping;
        tupleListToMap(mapping, data);
        addToRHS(rhs, mapping, assembler);
    }

//private
void OxleyDomain::assemblePDEBoundary(escript::AbstractSystemMatrix* mat,
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

#ifdef ESYS_HAVE_TRILINOS
        esys_trilinos::TrilinosMatrixAdapter* tm = dynamic_cast<esys_trilinos::TrilinosMatrixAdapter*>(mat);
        if (tm) {
            tm->resumeFill();
        }
#endif

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

    #ifdef ESYS_HAVE_TRILINOS
        if (tm) {
            tm->fillComplete(true);
        }
    #endif
}

//private
void OxleyDomain::assemblePDE(escript::AbstractSystemMatrix* mat,
                               escript::Data& rhs, const DataMap& coefs,
                               Assembler_ptr assembler) const
{
    if (rhs.isEmpty() && (isNotEmpty("X", coefs) 
        || isNotEmpty("du", coefs)) && isNotEmpty("Y", coefs))
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

#ifdef ESYS_HAVE_TRILINOS
    esys_trilinos::TrilinosMatrixAdapter* tm = dynamic_cast<esys_trilinos::TrilinosMatrixAdapter*>(mat);
    if (tm) {
        tm->resumeFill();
    }
#endif

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

void OxleyDomain::assemblePDEDirac(escript::AbstractSystemMatrix* mat,
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

void OxleyDomain::assemblePDEHanging(Teuchos::RCP<Tpetra::CrsMatrix<double,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>>* mat,
                                    Assembler_ptr assembler) const
{
    assembler->assemblePDEHanging(mat);
}

bool OxleyDomain::probeInterpolationAcross(int fsType_source,
                      const escript::AbstractDomain&, int fsType_target) const
{
    return false;
}

void OxleyDomain::updateSolutionInformation(escript::Data solution)
{
    throw OxleyException("programming error");
}

void OxleyDomain::updateMeshInformation()
{
    throw OxleyException("programming error");
}

escript::Data OxleyDomain::getUpdatedSolution()
{
    throw OxleyException("programming error");
}

#ifdef ESYS_HAVE_BOOST_NUMPY
boost::python::numpy::ndarray OxleyDomain::getNumpyX() const
{
    return continuousFunction(*this).getNumpyX();
}
#endif

} // end of namespace oxley


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

#include "Assemble.h"
#include "Util.h"

namespace finley {

using escript::ValueError;
using escript::NotImplementedError;

template<typename Scalar>
void Assemble_CopyNodalData(const NodeFile* nodes, escript::Data& out,
                            const escript::Data& in)
{
    if (!nodes)
        return;

    const int mpiSize = nodes->MPIInfo->size;
    const int numComps = out.getDataPointSize();
    const int in_data_type = in.getFunctionSpace().getTypeCode();
    const int out_data_type = out.getFunctionSpace().getTypeCode();

    // check out and in
    if (numComps != in.getDataPointSize()) {
        throw ValueError("Assemble_CopyNodalData: number of components of input and output Data do not match.");
    } else if (!out.actsExpanded()) {
        throw ValueError("Assemble_CopyNodalData: expanded Data object is expected for output data.");
    } else if (in.isComplex() != out.isComplex()) {
        throw ValueError("Assemble_CopyNodalData: complexity of input and output Data must match.");
    }

    // more sophisticated test needed for overlapping node/DOF counts
    if (in_data_type == FINLEY_NODES) {
        if (!in.numSamplesEqual(1, nodes->getNumNodes())) {
            throw ValueError("Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
    } else if (in_data_type == FINLEY_REDUCED_NODES) {
        if (!in.numSamplesEqual(1, nodes->getNumReducedNodes())) {
            throw ValueError("Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
    } else if (in_data_type == FINLEY_DEGREES_OF_FREEDOM) {
        if (!in.numSamplesEqual(1, nodes->getNumDegreesOfFreedom())) {
            throw ValueError("Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
        if (((out_data_type == FINLEY_NODES) || (out_data_type == FINLEY_DEGREES_OF_FREEDOM)) && !in.actsExpanded() && (mpiSize>1)) {
            throw ValueError("Assemble_CopyNodalData: FINLEY_DEGREES_OF_FREEDOM to FINLEY_NODES or FINLEY_DEGREES_OF_FREEDOM requires expanded input data on more than one processor.");
        }
    } else if (in_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (!in.numSamplesEqual(1, nodes->getNumReducedDegreesOfFreedom())) {
            throw ValueError("Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
        if ((out_data_type == FINLEY_DEGREES_OF_FREEDOM) && !in.actsExpanded() && (mpiSize>1)) {
            throw ValueError("Assemble_CopyNodalData: FINLEY_REDUCED_DEGREES_OF_FREEDOM to FINLEY_DEGREES_OF_FREEDOM requires expanded input data on more than one processor.");
        }
    } else {
        throw ValueError( "Assemble_CopyNodalData: illegal function space type for target object");
    }

    dim_t numOut = 0;
    switch (out_data_type) {
        case FINLEY_NODES:
            numOut = nodes->getNumNodes();
            break;

        case FINLEY_REDUCED_NODES:
            numOut = nodes->getNumReducedNodes();
            break;

        case FINLEY_DEGREES_OF_FREEDOM:
            numOut = nodes->getNumDegreesOfFreedom();
            break;

        case FINLEY_REDUCED_DEGREES_OF_FREEDOM:
            numOut = nodes->getNumReducedDegreesOfFreedom();
            break;

        default:
            throw ValueError("Assemble_CopyNodalData: illegal function space type for source object");
    }

    if (!out.numSamplesEqual(1, numOut)) {
        throw ValueError("Assemble_CopyNodalData: illegal number of samples of output Data object");
    }

    const Scalar zero = static_cast<Scalar>(0);
    const size_t numComps_size = numComps * sizeof(Scalar);

    /**************************** FINLEY_NODES ******************************/
    if (in_data_type == FINLEY_NODES) {
        out.requireWrite();
        if (out_data_type == FINLEY_NODES) {
#pragma omp parallel for
            for (index_t n = 0; n < numOut; n++) {
                memcpy(out.getSampleDataRW(n, zero), in.getSampleDataRO(n, zero), numComps_size);
            }
        } else if (out_data_type == FINLEY_REDUCED_NODES) {
            const IndexVector& map = nodes->borrowReducedNodesTarget();
            const dim_t mapSize = map.size();
#pragma omp parallel for
            for (index_t n = 0; n < mapSize; n++) {
                memcpy(out.getSampleDataRW(n, zero), in.getSampleDataRO(map[n], zero),
                       numComps_size);
            }
        } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
            const IndexVector& map = nodes->borrowDegreesOfFreedomTarget();
#pragma omp parallel for
            for (index_t n = 0; n < numOut; n++) {
                memcpy(out.getSampleDataRW(n, zero), in.getSampleDataRO(map[n], zero),
                       numComps_size);
            }
        } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            const std::vector<index_t>& map = nodes->borrowReducedDegreesOfFreedomTarget();
#pragma omp parallel for
            for (index_t n=0; n<numOut; n++) {
                memcpy(out.getSampleDataRW(n, zero), in.getSampleDataRO(map[n], zero),
                       numComps_size);
            }
        }

    /************************ FINLEY_REDUCED_NODES **************************/
    } else if (in_data_type == FINLEY_REDUCED_NODES) {
        if (out_data_type == FINLEY_NODES) {
            throw ValueError("Assemble_CopyNodalData: cannot copy from reduced nodes to nodes.");
        } else if (out_data_type == FINLEY_REDUCED_NODES) {
            out.requireWrite();
            const dim_t nNodes = nodes->getNumNodes();
#pragma omp parallel for
            for (index_t n=0; n < nNodes; n++) {
                memcpy(out.getSampleDataRW(n, zero), in.getSampleDataRO(n, zero), numComps_size);
            }
       } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
            throw ValueError("Assemble_CopyNodalData: cannot copy from reduced nodes to degrees of freedom.");
       } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            out.requireWrite();
            const index_t* target = nodes->borrowTargetReducedNodes();
            const IndexVector& map = nodes->borrowReducedDegreesOfFreedomTarget();
#pragma omp parallel for
            for (index_t n = 0; n < numOut; n++) {
               memcpy(out.getSampleDataRW(n, zero),
                      in.getSampleDataRO(target[map[n]], zero), numComps_size);
            }
        }
    /********************** FINLEY_DEGREES_OF_FREEDOM ***********************/
    } else if (in_data_type == FINLEY_DEGREES_OF_FREEDOM) {
        out.requireWrite();
        if (out_data_type == FINLEY_NODES) {
            const_cast<escript::Data*>(&in)->resolve();
            const index_t* target = nodes->borrowTargetDegreesOfFreedom();
#ifdef ESYS_HAVE_TRILINOS
            using namespace esys_trilinos;

            Teuchos::RCP<const MapType> colMap;
            Teuchos::RCP<const MapType> rowMap;
            MapType colPointMap;
            MapType rowPointMap;
            if (numComps > 1) {
                colPointMap = BlockVectorType<Scalar>::makePointMap(
                                             *nodes->trilinosColMap, numComps);
                rowPointMap = BlockVectorType<Scalar>::makePointMap(
                                             *nodes->trilinosRowMap, numComps);
                colMap = Teuchos::rcpFromRef(colPointMap);
                rowMap = Teuchos::rcpFromRef(rowPointMap);
            } else {
                colMap = nodes->trilinosColMap;
                rowMap = nodes->trilinosRowMap;
            }

            const ImportType importer(rowMap, colMap);
            const Teuchos::ArrayView<const Scalar> localIn(
                                               in.getSampleDataRO(0, zero),
                                               in.getNumDataPoints()*numComps);
            Teuchos::RCP<VectorType<Scalar> > lclData =
                rcp(new VectorType<Scalar>(rowMap, localIn, localIn.size(), 1));
            Teuchos::RCP<VectorType<Scalar> > gblData = rcp(new VectorType<Scalar>(colMap, 1));
            gblData->doImport(*lclData, importer, Tpetra::INSERT);
            Teuchos::ArrayRCP<const Scalar> gblArray(gblData->getData(0));
#pragma omp parallel for
            for (index_t i = 0; i < numOut; i++) {
                const Scalar* src = &gblArray[target[i] * numComps];
                std::copy(src, src+numComps, out.getSampleDataRW(i, zero));
            }
#elif defined(ESYS_HAVE_PASO)
            paso::Coupler_ptr<Scalar> coupler(new paso::Coupler<Scalar>(
                        nodes->degreesOfFreedomConnector, numComps, nodes->MPIInfo));
            coupler->startCollect(in.getSampleDataRO(0, zero));
            const Scalar* recv_buffer = coupler->finishCollect();
            const index_t upperBound = nodes->getNumDegreesOfFreedom();
#pragma omp parallel for
            for (index_t n = 0; n < numOut; n++) {
                const index_t k = target[n];
                if (k < upperBound) {
                    memcpy(out.getSampleDataRW(n, zero), in.getSampleDataRO(k, zero),
                           numComps_size);
                } else {
                    memcpy(out.getSampleDataRW(n, zero),
                           &recv_buffer[(k - upperBound) * numComps],
                           numComps_size);
                }
            }
#endif // Trilinos / Paso
        } else if (out_data_type == FINLEY_REDUCED_NODES) {
            const_cast<escript::Data*>(&in)->resolve();
            const index_t* target = nodes->borrowTargetDegreesOfFreedom();
            const IndexVector& map = nodes->borrowReducedNodesTarget();
#ifdef ESYS_HAVE_TRILINOS
            using namespace esys_trilinos;

            Teuchos::RCP<const MapType> colMap;
            Teuchos::RCP<const MapType> rowMap;
            MapType colPointMap;
            MapType rowPointMap;
            if (numComps > 1) {
                colPointMap = BlockVectorType<Scalar>::makePointMap(
                                             *nodes->trilinosColMap, numComps);
                rowPointMap = BlockVectorType<Scalar>::makePointMap(
                                             *nodes->trilinosRowMap, numComps);
                colMap = Teuchos::rcpFromRef(colPointMap);
                rowMap = Teuchos::rcpFromRef(rowPointMap);
            } else {
                colMap = nodes->trilinosColMap;
                rowMap = nodes->trilinosRowMap;
            }

            const ImportType importer(rowMap, colMap);
            const Teuchos::ArrayView<const Scalar> localIn(
                                               in.getSampleDataRO(0, zero),
                                               in.getNumDataPoints()*numComps);
            Teuchos::RCP<VectorType<Scalar> > lclData = rcp(new VectorType<Scalar>(
                                               rowMap, localIn, localIn.size(), 1));
            Teuchos::RCP<VectorType<Scalar> > gblData = rcp(new VectorType<Scalar>(colMap, 1));
            gblData->doImport(*lclData, importer, Tpetra::INSERT);
            Teuchos::ArrayRCP<const Scalar> gblArray(gblData->getData(0));
#pragma omp parallel for
            for (index_t i = 0; i < numOut; i++) {
                const Scalar* src = &gblArray[target[map[i]] * numComps];
                std::copy(src, src+numComps, out.getSampleDataRW(i, zero));
            }
#elif defined(ESYS_HAVE_PASO)
            paso::Coupler_ptr<Scalar> coupler(new paso::Coupler<Scalar>(
                        nodes->degreesOfFreedomConnector, numComps, nodes->MPIInfo));
            coupler->startCollect(in.getSampleDataRO(0, zero));
            const Scalar* recv_buffer = coupler->finishCollect();
            const index_t upperBound = nodes->getNumDegreesOfFreedom();
            const dim_t mapSize = map.size();

#pragma omp parallel for
            for (index_t n = 0; n < mapSize; n++) {
                const index_t k = target[map[n]];
                if (k < upperBound) {
                    memcpy(out.getSampleDataRW(n, zero), in.getSampleDataRO(k, zero),
                           numComps_size);
                } else {
                    memcpy(out.getSampleDataRW(n, zero),
                           &recv_buffer[(k-upperBound)*numComps],
                           numComps_size);
                }
            }
#endif // Trilinos / Paso
        } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
#pragma omp parallel for
            for (index_t n = 0; n < numOut; n++) {
                memcpy(out.getSampleDataRW(n, zero), in.getSampleDataRO(n, zero),
                       numComps_size);
            }
        } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            const index_t* target = nodes->borrowTargetDegreesOfFreedom();
            const std::vector<index_t>& map = nodes->borrowReducedDegreesOfFreedomTarget();
#pragma omp parallel for
            for (index_t n=0; n<numOut; n++) {
                memcpy(out.getSampleDataRW(n, zero),
                       in.getSampleDataRO(target[map[n]], zero), numComps_size);
            }
        }

    /****************** FINLEY_REDUCED_DEGREES_OF_FREEDOM *******************/
    } else if (in_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (out_data_type == FINLEY_NODES) {
            throw ValueError("Assemble_CopyNodalData: cannot copy from reduced degrees of freedom to nodes.");
        } else if (out_data_type == FINLEY_REDUCED_NODES) {
            const_cast<escript::Data*>(&in)->resolve();
            const index_t* target = nodes->borrowTargetReducedDegreesOfFreedom();
            const IndexVector& map = nodes->borrowReducedNodesTarget();
            out.requireWrite();
#ifdef ESYS_HAVE_TRILINOS
            using namespace esys_trilinos;

            Teuchos::RCP<const MapType> colMap;
            Teuchos::RCP<const MapType> rowMap;
            MapType colPointMap;
            MapType rowPointMap;
            if (numComps > 1) {
                colPointMap = BlockVectorType<Scalar>::makePointMap(
                                      *nodes->trilinosReducedColMap, numComps);
                rowPointMap = BlockVectorType<Scalar>::makePointMap(
                                      *nodes->trilinosReducedRowMap, numComps);
                colMap = Teuchos::rcpFromRef(colPointMap);
                rowMap = Teuchos::rcpFromRef(rowPointMap);
            } else {
                colMap = nodes->trilinosReducedColMap;
                rowMap = nodes->trilinosReducedRowMap;
            }

            const ImportType importer(rowMap, colMap);
            const Teuchos::ArrayView<const Scalar> localIn(
                                               in.getSampleDataRO(0, zero),
                                               in.getNumDataPoints()*numComps);
            Teuchos::RCP<VectorType<Scalar> > lclData = rcp(new VectorType<Scalar>(
                                               rowMap, localIn, localIn.size(), 1));
            Teuchos::RCP<VectorType<Scalar> > gblData = rcp(new VectorType<Scalar>(colMap, 1));
            gblData->doImport(*lclData, importer, Tpetra::INSERT);
            Teuchos::ArrayRCP<const Scalar> gblArray(gblData->getData(0));
#pragma omp parallel for
            for (index_t i = 0; i < numOut; i++) {
                const Scalar* src = &gblArray[target[map[i]] * numComps];
                std::copy(src, src+numComps, out.getSampleDataRW(i, zero));
            }
#elif defined(ESYS_HAVE_PASO)
            paso::Coupler_ptr<Scalar> coupler(new paso::Coupler<Scalar>(
                        nodes->reducedDegreesOfFreedomConnector, numComps, nodes->MPIInfo));
            coupler->startCollect(in.getSampleDataRO(0, zero));
            const index_t upperBound = nodes->getNumReducedDegreesOfFreedom();
            const dim_t mapSize = map.size();
            const Scalar* recv_buffer = coupler->finishCollect();
#pragma omp parallel for
            for (index_t n = 0; n < mapSize; n++) {
                const index_t k = target[map[n]];
                if (k < upperBound) {
                    memcpy(out.getSampleDataRW(n, zero), in.getSampleDataRO(k, zero),
                           numComps_size);
                } else {
                    memcpy(out.getSampleDataRW(n, zero),
                           &recv_buffer[(k - upperBound) * numComps],
                           numComps_size);
                }
            }
#endif // Trilinos / Paso
        } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
            out.requireWrite();
#pragma omp parallel for
            for (index_t n = 0; n < numOut; n++) {
                memcpy(out.getSampleDataRW(n, zero), in.getSampleDataRO(n, zero), numComps_size);
            }
        } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
            throw ValueError("Assemble_CopyNodalData: cannot copy from reduced degrees of freedom to degrees of freedom.");
        }
    } // in_data_type
}

// instantiate our two supported versions
template void Assemble_CopyNodalData<escript::DataTypes::real_t>(
           const NodeFile* nodes, escript::Data& out, const escript::Data& in);
template void Assemble_CopyNodalData<escript::DataTypes::cplx_t>(
           const NodeFile* nodes, escript::Data& out, const escript::Data& in);

} // namespace finley


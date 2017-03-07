
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
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

#include "Assemble.h"
#include "Util.h"

namespace dudley {

void Assemble_CopyNodalData(const NodeFile* nodes, escript::Data& out,
                            const escript::Data& in)
{
    if (!nodes)
        return;

    if (in.isComplex() || out.isComplex())
    {
        throw escript::ValueError("Assemble_CopyNodalData: complex arguments not supported.");
    }
    const int mpiSize = nodes->MPIInfo->size;
    const int numComps = out.getDataPointSize();
    const int in_data_type = in.getFunctionSpace().getTypeCode();
    const int out_data_type = out.getFunctionSpace().getTypeCode();

    // check out and in
    if (numComps != in.getDataPointSize()) {
        throw escript::ValueError("Assemble_CopyNodalData: number of components of input and output Data do not match.");
    } else if (!out.actsExpanded()) {
        throw escript::ValueError("Assemble_CopyNodalData: expanded Data object is expected for output data.");
    }

    // more sophisticated test needed for overlapping node/DOF counts
    if (in_data_type == DUDLEY_NODES) {
        if (!in.numSamplesEqual(1, nodes->getNumNodes())) {
            throw escript::ValueError("Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
    } else if (in_data_type == DUDLEY_DEGREES_OF_FREEDOM) {
        if (!in.numSamplesEqual(1, nodes->getNumDegreesOfFreedom())) {
            throw escript::ValueError("Assemble_CopyNodalData: illegal number of samples of input Data object");
        }
        if ((((out_data_type == DUDLEY_NODES) || (out_data_type == DUDLEY_DEGREES_OF_FREEDOM)) && !in.actsExpanded() && (mpiSize > 1))) {

            throw DudleyException("Assemble_CopyNodalData: DUDLEY_DEGREES_OF_FREEDOM to DUDLEY_NODES or DUDLEY_DEGREES_OF_FREEDOM requires expanded input data on more than one processor.");
        }
    } else {
        throw escript::ValueError( "Assemble_CopyNodalData: illegal function space type for target object");
    }

    dim_t numOut = 0;
    switch (out_data_type) {
        case DUDLEY_NODES:
            numOut = nodes->getNumNodes();
            break;

        case DUDLEY_DEGREES_OF_FREEDOM:
            numOut = nodes->getNumDegreesOfFreedom();
            break;

        default:
            throw escript::ValueError("Assemble_CopyNodalData: illegal function space type for source object");
    }

    if (!out.numSamplesEqual(1, numOut)) {
        throw escript::ValueError("Assemble_CopyNodalData: illegal number of samples of output Data object");
    }

    const size_t numComps_size = numComps * sizeof(double);

    /**************************** DUDLEY_NODES ******************************/
    if (in_data_type == DUDLEY_NODES) {
        out.requireWrite();
        if (out_data_type == DUDLEY_NODES) {
#pragma omp parallel for
            for (index_t n = 0; n < numOut; n++) {
                memcpy(out.getSampleDataRW(n, static_cast<escript::DataTypes::real_t>(0)), in.getSampleDataRO(n, static_cast<escript::DataTypes::real_t>(0)), numComps_size);
            }
        } else if (out_data_type == DUDLEY_DEGREES_OF_FREEDOM) {
            const index_t* map = nodes->borrowDegreesOfFreedomTarget();
#pragma omp parallel for
            for (index_t n = 0; n < numOut; n++) {
                memcpy(out.getSampleDataRW(n, static_cast<escript::DataTypes::real_t>(0)), in.getSampleDataRO(map[n], static_cast<escript::DataTypes::real_t>(0)),
                       numComps_size);
            }
        }
    /********************** DUDLEY_DEGREES_OF_FREEDOM ***********************/
    } else if (in_data_type == DUDLEY_DEGREES_OF_FREEDOM) {
        out.requireWrite();
        if (out_data_type == DUDLEY_NODES) {
            const_cast<escript::Data*>(&in)->resolve();
            const index_t* target = nodes->borrowTargetDegreesOfFreedom();
#ifdef ESYS_HAVE_PASO
            paso::Coupler_ptr coupler(new paso::Coupler(nodes->degreesOfFreedomConnector, numComps, nodes->MPIInfo));
            coupler->startCollect(in.getDataRO());
            const double* recv_buffer = coupler->finishCollect();
            const index_t upperBound = nodes->getNumDegreesOfFreedom();
#pragma omp parallel for
            for (index_t n = 0; n < numOut; n++) {
                const index_t k = target[n];
                if (k < upperBound) {
                    memcpy(out.getSampleDataRW(n, static_cast<escript::DataTypes::real_t>(0)), in.getSampleDataRO(k, static_cast<escript::DataTypes::real_t>(0)),
                           numComps_size);
                } else {
                    memcpy(out.getSampleDataRW(n, static_cast<escript::DataTypes::real_t>(0)),
                           &recv_buffer[(k - upperBound) * numComps],
                           numComps_size);
                }
            }
#elif defined(ESYS_HAVE_TRILINOS)
            using namespace esys_trilinos;

            const_TrilinosGraph_ptr graph(nodes->getTrilinosGraph());
            Teuchos::RCP<const MapType> colMap;
            Teuchos::RCP<const MapType> rowMap;
            MapType colPointMap;
            MapType rowPointMap;
            if (numComps > 1) {
                colPointMap = RealBlockVector::makePointMap(*graph->getColMap(),
                                                            numComps);
                rowPointMap = RealBlockVector::makePointMap(*graph->getRowMap(),
                                                            numComps);
                colMap = Teuchos::rcpFromRef(colPointMap);
                rowMap = Teuchos::rcpFromRef(rowPointMap);
            } else {
                colMap = graph->getColMap();
                rowMap = graph->getRowMap();
            }

            const ImportType importer(rowMap, colMap);
            const Teuchos::ArrayView<const real_t> localIn(
                                               in.getSampleDataRO(0, static_cast<escript::DataTypes::real_t>(0)),
                                               in.getNumDataPoints()*numComps);
            Teuchos::RCP<RealVector> lclData = rcp(new RealVector(rowMap,
                                                  localIn, localIn.size(), 1));
            Teuchos::RCP<RealVector> gblData = rcp(new RealVector(colMap, 1));
            gblData->doImport(*lclData, importer, Tpetra::INSERT);
            Teuchos::ArrayRCP<const real_t> gblArray(gblData->getData(0));
#pragma omp parallel for
            for (index_t i = 0; i < numOut; i++) {
                const real_t* src = &gblArray[target[i] * numComps];
                std::copy(src, src+numComps, out.getSampleDataRW(i));
            }
#endif
        } else if (out_data_type == DUDLEY_DEGREES_OF_FREEDOM) {
#pragma omp parallel for
            for (index_t n = 0; n < numOut; n++) {
                memcpy(out.getSampleDataRW(n, static_cast<escript::DataTypes::real_t>(0)), in.getSampleDataRO(n, static_cast<escript::DataTypes::real_t>(0)),
                       numComps_size);
            }
        }
    } // in_data_type
}

} // namespace dudley


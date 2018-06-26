
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
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

#include "DudleyDomain.h"
#include "IndexList.h"

#include <boost/scoped_array.hpp>

namespace dudley {

#ifdef ESYS_HAVE_PASO
paso::SystemMatrixPattern_ptr DudleyDomain::getPasoPattern() const
{
    // make sure that the pattern is available
    if (!pasoPattern)
        pasoPattern = makePasoPattern();

    return pasoPattern;
}

paso::SystemMatrixPattern_ptr DudleyDomain::makePasoPattern() const
{
    const dim_t myNumTargets = m_nodes->getNumDegreesOfFreedom();
    const dim_t numTargets = m_nodes->getNumDegreesOfFreedomTargets();
    const index_t* target = m_nodes->borrowTargetDegreesOfFreedom();
    boost::scoped_array<IndexList> index_list(new IndexList[numTargets]);

#pragma omp parallel
    {
        // insert contributions from element matrices into columns in indexlist
        IndexList_insertElements(index_list.get(), m_elements, target);
        IndexList_insertElements(index_list.get(), m_faceElements, target);
        IndexList_insertElements(index_list.get(), m_points, target);
    }

    // create pattern
    paso::Pattern_ptr mainPattern(paso::Pattern::fromIndexListArray(0,
              myNumTargets, index_list.get(), 0, myNumTargets, 0));
    paso::Pattern_ptr colCouplePattern(paso::Pattern::fromIndexListArray(0,
              myNumTargets, index_list.get(), myNumTargets, numTargets,
              -myNumTargets));
    paso::Pattern_ptr rowCouplePattern(paso::Pattern::fromIndexListArray(
              myNumTargets, numTargets, index_list.get(), 0, myNumTargets, 0));

    paso::Connector_ptr connector(m_nodes->degreesOfFreedomConnector);
    paso::SystemMatrixPattern_ptr out(new paso::SystemMatrixPattern(
                MATRIX_FORMAT_DEFAULT, m_nodes->dofDistribution,
                m_nodes->dofDistribution, mainPattern, colCouplePattern,
                rowCouplePattern, connector, connector));
    return out;
}
#endif // ESYS_HAVE_PASO

} // namespace dudley


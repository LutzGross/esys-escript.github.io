
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


/****************************************************************************

  Dudley: Mesh

*****************************************************************************/

#ifdef USE_TRILINOS

#include "Mesh.h"
#include "IndexList.h"

#include <boost/scoped_array.hpp>

using namespace esys_trilinos;

namespace dudley {

typedef std::vector<index_t> IndexVector;

esys_trilinos::const_TrilinosGraph_ptr createTrilinosGraph(Dudley_Mesh* mesh)
{
    throw DudleyException("createTrilinosGraph: Not implemented");
}

} // namespace dudley

#endif // USE_TRILINOS


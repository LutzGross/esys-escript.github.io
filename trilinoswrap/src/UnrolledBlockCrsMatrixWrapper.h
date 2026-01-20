
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __ESYS_TRILINOS_UNROLLEDBLOCKCRSMATRIXWRAPPER_H__
#define __ESYS_TRILINOS_UNROLLEDBLOCKCRSMATRIXWRAPPER_H__

#include <trilinoswrap/CrsMatrixWrapper.h>

namespace esys_trilinos {

template<typename ST>
class UnrolledBlockCrsMatrixWrapper : public CrsMatrixWrapper<ST>
{
public:
    typedef typename CrsMatrixWrapper<ST>::Matrix Matrix;

    /**
       \brief
       Creates a new Trilinos CRS matrix wrapper using a compatible
       fill-complete unrolled Trilinos matrix graph and given block size.
    */
    UnrolledBlockCrsMatrixWrapper(const_TrilinosGraph_ptr graph, int blocksize);

    void add(const std::vector<LO>& rowIndex, const std::vector<ST>& array);

private:
    int blockSize;
};

} // namespace esys_trilinos

#endif // __ESYS_TRILINOS_UNROLLEDBLOCKCRSMATRIXWRAPPER_H__


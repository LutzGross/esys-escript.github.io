
/*****************************************************************************
*
* Copyright (c) 2016 by The University of Queensland
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

#include "UnrolledBlockCrsMatrixWrapper.h" 

#include <escript/index.h>

namespace esys_trilinos {

template<typename ST>
UnrolledBlockCrsMatrixWrapper<ST>::UnrolledBlockCrsMatrixWrapper(
        const_TrilinosGraph_ptr graph, int blocksize) :
    CrsMatrixWrapper<ST>(graph),
    blockSize(blocksize)
{
}

template<typename ST>
void UnrolledBlockCrsMatrixWrapper<ST>::add(const std::vector<LO>& rowIdx,
                                            const std::vector<ST>& array)
{
    const size_t emSize = rowIdx.size();
    std::vector<LO> cols(emSize * blockSize);
    std::vector<ST> vals(emSize * blockSize);
    for (size_t ri = 0; ri < emSize; ri++) {
        for (int rj = 0; rj < blockSize; rj++) {
            const LO row = rowIdx[ri] * blockSize + rj;
            if (row <= this->maxLocalRow) {
                for (int ci = 0; ci < emSize; ci++) {
                    for (int cj = 0; cj < blockSize; cj++) {
                        cols[ci*blockSize + cj] = rowIdx[ci] * blockSize + cj;
                        const size_t srcIdx = INDEX4(rj, cj, ri, ci, blockSize, blockSize, emSize);
                        vals[ci*blockSize + cj] = array[srcIdx];
                    }
                }
                this->mat.sumIntoLocalValues(row, cols, vals);
            }
        }
    }
}

// instantiate
template class UnrolledBlockCrsMatrixWrapper<real_t>;
template class UnrolledBlockCrsMatrixWrapper<cplx_t>;

}  // end of namespace


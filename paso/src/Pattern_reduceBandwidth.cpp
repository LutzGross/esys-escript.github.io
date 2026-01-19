
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/


/****************************************************************************/

/*   Paso: pattern: a simple algorithm to create a new labeling */
/*   of input/output with a minimum bandwidth. The algorithm    */
/*   starts from a vertex with a minimum number of neighbours   */
/*   (connections in the pattern). From this root it searches   */
/*   recursively the neighbours. In the last level a new root   */
/*   is picked (again the vertex with the minimum number of     */
/*   neighbours) and the process is repeated. The algorithm is  */
/*   repeated until the maximum number of vertices in a level   */
/*   cannot be reduced.                                         */

/****************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Pattern.h"
#include "PasoException.h"

namespace paso {

// calculate initial bandwidth for a given labeling
dim_t Pattern::getBandwidth(index_t* label) const
{
    dim_t bandwidth = 0;

#pragma omp parallel
    {
        dim_t local_bandwidth=0;
#pragma omp for
        for (dim_t i=0; i < numOutput; ++i) {
            const index_t k = label[i];
            for (index_t iptr=ptr[i]; iptr < ptr[i+1]; ++iptr) {
                local_bandwidth = std::max(local_bandwidth,std::abs(k-label[index[iptr]]));
            }
        }
#pragma omp critical
        {
            bandwidth=std::max(local_bandwidth, bandwidth);
        }
    }
    return bandwidth;
}

// structure used to create an ordering by increasing degree
struct DegreeAndIdx
{
    dim_t deg;
    index_t idx;
};

int comparDegreeAndIdx(const void *arg1, const void *arg2)
{
    DegreeAndIdx a1=*(DegreeAndIdx*)arg1;
    DegreeAndIdx a2=*(DegreeAndIdx*)arg2;
    if (a1.deg < a2.deg) {
        return -1;
    } else if (a1.deg > a2.deg) {
        return 1;
    } else if (a1.idx < a2.idx) {
        return -1;
    } else if (a1.idx > a2.idx) {
        return 1;
    }
    return 0;
}

/*  dropTree() drops a tree in pattern from root
 *  root                 - on input the starting point of the tree.
 *  AssignedLevel        - array of length numOutput indicating the level
 *                         assigned to vertex
 *  VerticesInTree       - on output contains the vertices used in tree
 *                         (array of length numOutput)
 *  numLevels            - on output the number of levels used
 *  firstVertexInLevel   - on output firstVertexInLevel[i] points to the first
 *                         vertex of level i in VerticesInTree
 *                         firstVertexInLevel[i+1]-firstVertexInLevel[i] is the
 *                         number of vertices in level i (array of length
 *                         numOutput+1)
 *  max_LevelWidth_abort - input param which triggers early return if
 *                         LevelWidth becomes >= max_LevelWidth_abort
 */
bool dropTree(index_t root, const Pattern* pattern, index_t* AssignedLevel,
              index_t* VerticesInTree, dim_t* numLevels,
              index_t* firstVertexInLevel, dim_t max_LevelWidth_abort,
              dim_t N)
{
#pragma omp parallel for
    for (dim_t i=0; i < pattern->numInput; ++i)
        AssignedLevel[i]=-1;

    dim_t nlvls = 0;
    AssignedLevel[root] = 0;
    VerticesInTree[0] = root;
    firstVertexInLevel[0]=0;
    dim_t level_top = firstVertexInLevel[0]+1;

    while (firstVertexInLevel[nlvls] < level_top) {
        nlvls++;
        firstVertexInLevel[nlvls] = level_top;
        if (firstVertexInLevel[nlvls]-firstVertexInLevel[nlvls-1] >= max_LevelWidth_abort)
            return false;

        for (dim_t i=firstVertexInLevel[nlvls-1]; i < firstVertexInLevel[nlvls]; ++i) {
            const index_t k = VerticesInTree[i];
            for (index_t j=pattern->ptr[k]; j<pattern->ptr[k+1]; ++j) {
                const index_t itest = pattern->index[j];
                if (AssignedLevel[itest] < 0) {
#ifdef BOUNDS_CHECK
                    ESYS_ASSERT(itest >= 0 && itest < N, "BOUNDS_CHECK: itest=" << itest << ", N=" << N);
                    ESYS_ASSERT(level_top >= 0 && level_top < N, "BOUNDS_CHECK: level_top=" << level_top << ", N=" << N);
#endif
                    AssignedLevel[itest] = nlvls;
                    VerticesInTree[level_top] = itest;
                    level_top++;
                }
            }
        }
    }
    *numLevels=nlvls;
    return true;
}

/* the driver */
void Pattern::reduceBandwidth(index_t* oldToNew)
{
    if (numOutput != numInput) {
        throw PasoException("Pattern::reduceBandwidth: pattern needs to be for a square matrix.");
    } else if (numOutput == 0) {
        return;
    }

    dim_t i, N=numOutput, numLevels, deg, bandwidth;
    index_t k;

    /* printf("relabeling of %d DOFs started.\n",N); */
    DegreeAndIdx* degAndIdx = new DegreeAndIdx[N];
    index_t* oldLabel = new index_t[N];
    index_t* AssignedLevel = new index_t[N];
    index_t* VerticesInTree = new index_t[N];
    index_t* firstVertexInLevel = new index_t[N+1];

    // get the initial bandwidth
#pragma omp parallel for private(i)
    for (i=0; i<N; ++i) oldToNew[i]=i;

    dim_t initial_bandwidth = getBandwidth(oldToNew);
    // printf("initial bandwidth = %d\n",initial_bandwidth);

    #pragma omp parallel for private(i)
    for (i=0; i<N; ++i) {
        oldToNew[i] = -1;
        degAndIdx[i].idx = i;
        degAndIdx[i].deg = ptr[i+1] - ptr[i];
    }

    // create an ordering with increasing degree
    qsort(degAndIdx, (size_t)N, sizeof(DegreeAndIdx), comparDegreeAndIdx);
    index_t root = degAndIdx[0].idx;
    dim_t numLabeledVertices = 0;

    while (root >= 0) {
        dim_t max_LevelWidth = N+1;
        dim_t numVerticesInTree = 0;
        while (dropTree(root, this, AssignedLevel, VerticesInTree,
                        &numLevels, firstVertexInLevel, max_LevelWidth, N)) {
            // find new maximum level width
            max_LevelWidth=0;
#ifdef BOUNDS_CHECK
            ESYS_ASSERT(numLevels <= N, "BOUNDS_CHECK: numLevels=" << numLevels << ", N=" << N);
#endif
            for (i = 0; i < numLevels; ++i) {
                max_LevelWidth=std::max(max_LevelWidth, firstVertexInLevel[i+1]-firstVertexInLevel[i]);
            }
            // find a vertex in the last level which has minimum degree
            dim_t min_deg = N+1;
            root=-1;
            for (i = firstVertexInLevel[numLevels-1]; i<firstVertexInLevel[numLevels];++i) {
                k=VerticesInTree[i];
                deg=ptr[k+1]-ptr[k];
                if (deg<min_deg) {
                    min_deg=deg;
                    root=k;
                }
            }
            // save the vertices in the current tree
            numVerticesInTree=firstVertexInLevel[numLevels];
#ifdef BOUNDS_CHECK
            ESYS_ASSERT(numLabeledVertices+firstVertexInLevel[numLevels] <= N,
                    "BOUNDS_CHECK: numLabeledVertices=" << numLabeledVertices
                    << ", root=" << root << ", N=" << N
                    << ", first[numLevels]=" << firstVertexInLevel[numLevels]);
#endif
            for (i = 0; i < firstVertexInLevel[numLevels]; ++i) {
                oldLabel[numLabeledVertices+i]=VerticesInTree[i];
            }
        }
#ifdef BOUNDS_CHECK
            ESYS_ASSERT(numLabeledVertices+numVerticesInTree <= N,
                    "BOUNDS_CHECK: numLabeledVertices=" << numLabeledVertices
                    << ", root=" << root << ", N=" << N
                    << ", numVerticesInTree=" << numVerticesInTree);
#endif
        // now the vertices in the current tree
        for (i=0; i<numVerticesInTree; ++i) {
            oldToNew[oldLabel[numLabeledVertices+i]]=numLabeledVertices+i;
        }
        numLabeledVertices+=numVerticesInTree;
        // new search for a vertex which is not labeled yet
        root=-1;
        for (i=0; i<N; ++i) {
            if (oldToNew[degAndIdx[i].idx] < 0) {
                root=degAndIdx[i].idx;
                break;
            }
        }
    } // end of while root loop
    bandwidth=getBandwidth(oldToNew);
    // printf("bandwidth after DOF relabeling= %d\n",bandwidth);
    if (bandwidth >= initial_bandwidth) {
        //printf("initial labeling used.\n");
#pragma omp parallel for private(i)
        for (i=0; i<N; ++i) oldToNew[i]=i;
    }
    delete[] degAndIdx;
    delete[] oldLabel;
    delete[] AssignedLevel;
    delete[] VerticesInTree;
    delete[] firstVertexInLevel;
}

} // namespace paso

